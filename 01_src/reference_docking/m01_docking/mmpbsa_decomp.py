"""
MMPBSA Per-Residue Decomposition — Core Module (01g)
=======================================================
Prepares AMBER system and runs MMPBSA.py per-residue decomposition.

Fills the gap between:
  01d Footprint — vdW + ES in-vacuo per-residue (overestimates charged pairs)
  01f GB/SA     — global solvation score (no per-residue breakdown)

MMPBSA.py idecomp=1 gives: vdW + ES + GB + SA per residue, answering:
  - Are ARG598/LYS599 salt bridges real or in-vacuo artifacts?
  - Does ASP494 rescue with solvation screening?
  - Do hydrophobic residues (TRP392/TRP495) rise in relative ranking?

Two modes:
  single_point — 1-frame MMPBSA on best pose. Minutes, no error bars.
  md           — OpenMM explicit-solvent MD → strip → MMPBSA on trajectory.
                 Hours (5-10 ns), gives mean ± std per residue.

Pipeline (this module):
    1. Extract best pose from scored mol2
    2. Parametrize ligand (antechamber → GAFF2 + AM1-BCC → lib/frcmod)
    3. Build topologies with tleap (complex, receptor, ligand)
    4a. single_point: cpptraj 1-frame mdcrd
    4b. md: OpenMM solvate → minimize → heat → equilibrate → produce → strip
    5. MMPBSA.py with idecomp=1, igb=2

Analysis (01h mmpbsa_analysis.py):
    6. Parse FINAL_DECOMP_MMPBSA.dat → CSV (PDB numbering)
    7. Compare with 01d footprint (ranking delta, zone summary)

Location: 01_src/reference_docking/m01_docking/mmpbsa_decomp.py
Project: reference_docking
Module: 01g (core)
Version: 1.0 (2026-03-26)

References:
    Miller et al. J Chem Theory Comput 2012, 8(9):3314-21 (MMPBSA.py)
    Onufriev et al. Proteins 2004, 55(2):383-94 (OBC GB model, igb=2)
    Hou et al. J Chem Inf Model 2011, 51(1):69-82 (MM-GBSA best practices)
"""

import logging
import json
import re
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import parmed

logger = logging.getLogger(__name__)


# =============================================================================
# CONSTANTS — TEMPLATES
# =============================================================================

# tleap input for gas-phase topologies (complex, receptor, ligand)
TLEAP_GAS_TEMPLATE = """\
source leaprc.protein.{receptor_ff}
source leaprc.{ligand_ff}
source leaprc.water.tip3p

# Load ligand parameters
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}

# Load receptor
REC = loadpdb {receptor_pdb}

# Complex = receptor + ligand
COM = combine {{REC LIG}}

# Save gas-phase topologies
saveamberparm COM {complex_prmtop} {complex_inpcrd}
saveamberparm REC {receptor_prmtop} {receptor_inpcrd}
saveamberparm LIG {ligand_prmtop} {ligand_inpcrd}

quit
"""

# tleap input for solvated system (MD mode)
TLEAP_SOLVATED_TEMPLATE = """\
source leaprc.protein.{receptor_ff}
source leaprc.{ligand_ff}
source leaprc.water.tip3p

# Load ligand parameters
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}

# Load receptor
REC = loadpdb {receptor_pdb}

# Complex = receptor + ligand
COM = combine {{REC LIG}}

# Solvate with TIP3P
solvateoct COM TIP3PBOX {solvent_padding}

# Neutralize + add ions to target ionic strength
addionsrand COM Na+ 0
addionsrand COM Cl- 0
{extra_ions}

# Save solvated topology
saveamberparm COM {solvated_prmtop} {solvated_inpcrd}

quit
"""

# MMPBSA.py input file
MMPBSA_INPUT_TEMPLATE = """\
&general
  startframe={startframe}, endframe={endframe}, interval={interval},
  verbose=2, keep_files=0,
/
&gb
  igb={igb}, saltcon={saltcon},
/
&decomp
  idecomp={idecomp},
  dec_verbose=3,
  print_res="all"
/
"""

# cpptraj script for 1-frame trajectory (single_point mode)
CPPTRAJ_SINGLE_FRAME = """\
parm {prmtop}
trajin {inpcrd}
trajout {mdcrd} mdcrd nobox
go
quit
"""

# cpptraj script to strip water+ions from MD trajectory
CPPTRAJ_STRIP_TEMPLATE = """\
parm {solvated_prmtop}
trajin {trajectory} {first} {last} {offset}
strip :WAT,Na+,Cl-
trajout {dry_mdcrd} mdcrd nobox
go
quit
"""


# =============================================================================
# 1. POSE EXTRACTION
# =============================================================================

def extract_pose_from_mol2(
        scored_mol2: Union[str, Path],
        output_mol2: Union[str, Path],
        selection: str = "best_score",
        pose_index: int = 1,
) -> Dict[str, Any]:
    """
    Extract a single pose from a multi-pose DOCK6 scored mol2.

    Args:
        scored_mol2: Path to scored mol2 (multiple @<TRIPOS>MOLECULE blocks)
        output_mol2: Output path for single-pose mol2
        selection:   "best_score" (lowest Grid_Score) or "pose_index"
        pose_index:  1-based index (only if selection="pose_index")

    Returns:
        Dict with success, pose_index, grid_score, n_atoms
    """
    text = Path(scored_mol2).read_text()

    # DOCK6 scored mol2 has ####### metadata lines before the first
    # @<TRIPOS>MOLECULE. The Grid_Score in those header lines belongs
    # to the FIRST pose (the mol2 block that follows). We need to
    # associate those scores with the correct structural block.

    # Split into mol2 blocks (each starting with @<TRIPOS>MOLECULE)
    blocks = text.split("@<TRIPOS>MOLECULE")
    # First element is pre-header (####### lines), rest are mol2 blocks
    pre_header = blocks[0] if blocks else ""
    mol2_blocks = [b for b in blocks[1:] if b.strip()]

    if not mol2_blocks:
        return {"success": False, "error": "No poses found in scored mol2"}

    # Parse Grid_Score from each block
    # The pre-header scores belong to pose 1
    poses = []
    for i, block in enumerate(mol2_blocks):
        full_block = "@<TRIPOS>MOLECULE" + block

        # Only consider blocks that have actual atom data
        if "@<TRIPOS>ATOM" not in full_block:
            continue

        score = None
        # Check the block itself for scores
        for line in block.split("\n"):
            if "Grid_Score" in line:
                match = re.search(r"Grid_Score:\s+([-\d.]+)", line)
                if match:
                    score = float(match.group(1))
                    break

        # First pose: also check pre-header for score if not found in block
        if score is None and i == 0 and pre_header:
            for line in pre_header.split("\n"):
                if "Grid_Score" in line:
                    match = re.search(r"Grid_Score:\s+([-\d.]+)", line)
                    if match:
                        score = float(match.group(1))
                        break

        poses.append({"index": i + 1, "block": full_block, "grid_score": score})

    logger.info(f"  Found {len(poses)} poses in scored mol2")

    # Select pose
    if selection == "best_score":
        scored_poses = [p for p in poses if p["grid_score"] is not None]
        if not scored_poses:
            logger.warning("  No Grid_Score found, using first pose")
            selected = poses[0]
        else:
            selected = min(scored_poses, key=lambda p: p["grid_score"])
        logger.info(f"  Selected pose {selected['index']} "
                     f"(Grid_Score = {selected['grid_score']})")
    elif selection == "pose_index":
        if pose_index < 1 or pose_index > len(poses):
            return {"success": False,
                    "error": f"pose_index {pose_index} out of range (1-{len(poses)})"}
        selected = poses[pose_index - 1]
        logger.info(f"  Selected pose {selected['index']} by index")
    else:
        return {"success": False, "error": f"Unknown selection method: {selection}"}

    # Write single pose
    Path(output_mol2).parent.mkdir(parents=True, exist_ok=True)
    Path(output_mol2).write_text(selected["block"])

    # Count atoms
    n_atoms = 0
    in_atom = False
    for line in selected["block"].split("\n"):
        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            continue
        if line.startswith("@<TRIPOS>") and in_atom:
            break
        if in_atom and line.strip():
            n_atoms += 1

    return {
        "success": True,
        "pose_index": selected["index"],
        "grid_score": selected["grid_score"],
        "n_poses_total": len(poses),
        "n_atoms": n_atoms,
    }


# =============================================================================
# 2. LIGAND PARAMETRIZATION
# =============================================================================

def parametrize_ligand(
        ligand_mol2: Union[str, Path],
        output_dir: Union[str, Path],
        charge_method: str = "bcc",
        ligand_ff: str = "gaff2",
        timeout: int = 300,
) -> Dict[str, Any]:
    """
    Parametrize ligand with antechamber → parmchk2 for AMBER.

    Args:
        ligand_mol2:   Single-pose ligand mol2
        output_dir:    Directory for output files
        charge_method: "bcc" (AM1-BCC) or "gas" (Gasteiger)
        ligand_ff:     "gaff2" or "gaff"
        timeout:       Seconds for antechamber

    Returns:
        Dict with success, mol2_path, frcmod_path
    """
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    out_mol2 = output_dir / "ligand_amber.mol2"
    out_frcmod = output_dir / "ligand.frcmod"

    # --- Detect net charge from mol2 ---
    net_charge = _detect_charge_from_mol2(ligand_mol2)
    logger.info(f"  Ligand net charge: {net_charge}")

    # --- antechamber: mol2 → GAFF2 mol2 with AM1-BCC charges ---
    # Resolve to absolute path — antechamber runs with cwd=output_dir
    ligand_mol2_abs = str(Path(ligand_mol2).resolve())
    at_type = "gaff2" if ligand_ff == "gaff2" else "gaff"
    cmd_ante = [
        "antechamber",
        "-i", ligand_mol2_abs,
        "-fi", "mol2",
        "-o", str(out_mol2),
        "-fo", "mol2",
        "-c", charge_method,
        "-at", at_type,
        "-nc", str(net_charge),
        "-pf", "y",
        "-dr", "no",
    ]

    logger.info(f"  antechamber: {charge_method} charges, {at_type} types")
    logger.debug(f"  cmd: {' '.join(cmd_ante)}")

    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd_ante, capture_output=True, text=True,
            timeout=timeout, cwd=str(output_dir),
        )
        runtime = time.time() - t0

        if proc.returncode != 0 or not out_mol2.exists():
            logger.warning(f"  antechamber failed ({runtime:.1f}s)")
            logger.debug(f"  stderr: {proc.stderr[:500]}")

            # Fallback: try Gasteiger if BCC failed
            if charge_method == "bcc":
                logger.info("  Retrying with Gasteiger charges...")
                cmd_ante[cmd_ante.index("bcc")] = "gas"
                proc = subprocess.run(
                    cmd_ante, capture_output=True, text=True,
                    timeout=timeout, cwd=str(output_dir),
                )
                if proc.returncode != 0 or not out_mol2.exists():
                    return {"success": False,
                            "error": f"antechamber failed (bcc + gas): {proc.stderr[:200]}"}
                logger.info("  Gasteiger fallback OK")
            else:
                return {"success": False,
                        "error": f"antechamber failed: {proc.stderr[:200]}"}

    except subprocess.TimeoutExpired:
        return {"success": False,
                "error": f"antechamber timeout ({timeout}s)"}

    # --- parmchk2: generate missing parameters ---
    cmd_parmchk = [
        "parmchk2",
        "-i", str(out_mol2),
        "-f", "mol2",
        "-o", str(out_frcmod),
        "-s", at_type,
    ]

    logger.info("  parmchk2: generating frcmod")
    proc = subprocess.run(
        cmd_parmchk, capture_output=True, text=True,
        timeout=60, cwd=str(output_dir),
    )
    if proc.returncode != 0 or not out_frcmod.exists():
        return {"success": False,
                "error": f"parmchk2 failed: {proc.stderr[:200]}"}

    return {
        "success": True,
        "mol2_path": str(out_mol2),
        "frcmod_path": str(out_frcmod),
        "charge_method": charge_method if out_mol2.exists() else "gas",
        "atom_types": at_type,
    }


def _detect_charge_from_mol2(mol2_path: Union[str, Path]) -> int:
    """Detect net formal charge from mol2 partial charges."""
    text = Path(mol2_path).read_text()
    charges = []
    in_atom = False
    for line in text.split("\n"):
        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            continue
        if line.startswith("@<TRIPOS>") and in_atom:
            break
        if in_atom and line.strip():
            parts = line.split()
            if len(parts) >= 9:
                try:
                    charges.append(float(parts[8]))
                except (ValueError, IndexError):
                    pass
    if charges:
        total = sum(charges)
        return round(total)
    return 0


# =============================================================================
# 3. TLEAP TOPOLOGY BUILD
# =============================================================================

def _sanitize_pdb_for_tleap(
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
) -> Dict[str, Any]:
    """
    Rename protonation-state-dependent residues for AMBER ff14SB / tleap.

    Uses parmed to READ the PDB structurally (no column edge cases), then
    applies renames to the ORIGINAL PDB lines (no reformatting, no
    renumbering, no terminal name insertion). Best of both worlds.

    Renames applied (based on atom inventory per residue):
        PDB2PQR terminals: NASP/NGLU/NHIS/etc. → strip N/C prefix
        HIS + HD1 + HE2  →  HIP   (doubly protonated)
        HIS + HD1 only   →  HID   (δ-protonated)
        HIS + HE2 only   →  HIE   (ε-protonated)
        GLU + HE2         →  GLH   (protonated carboxylate)
        ASP + HD2         →  ASH   (protonated carboxylate)
        CYS (SG bonded to SG) → CYX (disulfide)
    """
    STD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
              "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
              "THR", "TRP", "TYR", "VAL"}
    HIS_NAMES = {"HIS", "HIE", "HID", "HIP", "HSE", "HSD", "HSP"}

    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)

    # ── Phase 1: Use parmed to analyze structure ──
    pdb = parmed.load_file(str(input_pdb))
    renames = []

    # Build rename map: (chain, resnum, orig_name) → new_name
    # Use residue index to map back to original PDB lines
    rename_by_idx = {}  # parmed residue index → new 3-letter name

    for i, residue in enumerate(pdb.residues):
        atoms = {a.name for a in residue.atoms}
        orig = residue.name

        # PDB2PQR terminal prefixes (NASP, NGLU, CALA, etc.)
        if len(orig) == 4 and orig[0] in ("N", "C") and orig[1:] in STD_AA:
            std = orig[1:]
            renames.append(f"{orig}{residue.number} → {std} (terminal)")
            rename_by_idx[i] = std
            orig = std  # use stripped name for protonation checks below

        # Histidines
        if orig in HIS_NAMES:
            hd1, he2 = "HD1" in atoms, "HE2" in atoms
            new = "HIP" if (hd1 and he2) else "HID" if hd1 else "HIE"
            if new != orig:
                renames.append(f"{orig}{residue.number} → {new}")
                rename_by_idx[i] = new

        # Protonated glutamate
        elif orig == "GLU" and "HE2" in atoms:
            renames.append(f"GLU{residue.number} → GLH")
            rename_by_idx[i] = "GLH"

        # Protonated aspartate
        elif orig == "ASP" and "HD2" in atoms:
            renames.append(f"ASP{residue.number} → ASH")
            rename_by_idx[i] = "ASH"

        # Disulfide cysteine
        elif orig == "CYS":
            sg_atoms = [a for a in residue.atoms if a.name == "SG"]
            if sg_atoms:
                sg = sg_atoms[0]
                bonded_to_sg = any(
                    (b.atom1.name == "SG" and b.atom1.residue != residue) or
                    (b.atom2.name == "SG" and b.atom2.residue != residue)
                    for b in sg.bonds
                )
                if bonded_to_sg:
                    renames.append(f"CYS{residue.number} → CYX")
                    rename_by_idx[i] = "CYX"

    # ── Identify N-terminal residue indices ──
    n_terminal_idxs = set()
    chains_seen = set()
    for i, residue in enumerate(pdb.residues):
        chain = residue.chain or "A"
        if chain not in chains_seen:
            chains_seen.add(chain)
            n_terminal_idxs.add(i)

    # ── Phase 2: Map parmed residue indices to original PDB lines ──
    # Walk the original PDB and track which residue index each ATOM line
    # belongs to, by counting residue transitions the same way parmed does.
    original_lines = Path(input_pdb).read_text().split("\n")
    out_lines = []
    res_idx = -1
    prev_res_key = None

    for line in original_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Extract residue identity from the line
            res_name_raw = line[17:21]  # up to 4 chars (handles NASP etc.)
            chain = line[21:22]
            res_num_raw = line[22:26]
            res_key = (res_name_raw, chain, res_num_raw)

            if res_key != prev_res_key:
                res_idx += 1
                prev_res_key = res_key

            if res_idx in rename_by_idx:
                new_name = rename_by_idx[res_idx]
                # Write new name into columns 17-20 (right-justified in 3 chars)
                # If original was 4-char (NASP), overwrite cols 17-20 with " XXX"
                line = line[:17] + new_name.rjust(3) + " " + line[21:]

            # N-terminal: rename backbone H → H1 (PDB2PQR writes H,
            # AMBER ff14SB expects H1 for N-terminal NH3+)
            if res_idx in n_terminal_idxs:
                atom_name = line[12:16]
                if atom_name.strip() == "H":
                    line = line[:12] + " H1 " + line[16:]

        out_lines.append(line)

    n_hfix = len(n_terminal_idxs)
    if n_hfix:
        renames.append(f"N-terminal H → H1 ({n_hfix} chain(s))")

    Path(output_pdb).write_text("\n".join(out_lines))

    if renames:
        logger.info(f"  PDB → tleap: {len(renames)} residue(s) renamed (parmed-guided):")
        for r in renames:
            logger.info(f"    {r}")
    else:
        logger.info("  PDB → tleap: no renames needed")

    return {"renames": renames}


def build_topologies(
        receptor_pdb: Union[str, Path],
        ligand_mol2: Union[str, Path],
        ligand_frcmod: Union[str, Path],
        output_dir: Union[str, Path],
        receptor_ff: str = "ff14SB",
        ligand_ff: str = "gaff2",
        solvate: bool = False,
        solvent_padding: float = 10.0,
        ionic_strength: float = 0.15,
) -> Dict[str, Any]:
    """
    Build AMBER topologies with tleap.

    Gas-phase topologies (complex, receptor, ligand) are always built —
    these are required by MMPBSA.py. Solvated topology is only built
    when mode=md.

    Args:
        receptor_pdb:    Protonated receptor PDB (from 00b)
        ligand_mol2:     GAFF2-parametrized ligand mol2
        ligand_frcmod:   Missing parameters file
        output_dir:      Output directory
        receptor_ff:     Protein force field
        ligand_ff:       Ligand force field
        solvate:         If True, also build solvated system (for MD)
        solvent_padding: TIP3P box padding in Å
        ionic_strength:  NaCl concentration in M

    Returns:
        Dict with success and paths to all topology files
    """
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Sanitize PDB for tleap (histidine naming, etc.) ---
    sanitized_pdb = output_dir / "receptor_tleap.pdb"
    san_fixes = _sanitize_pdb_for_tleap(receptor_pdb, sanitized_pdb)

    # --- Gas-phase topologies (always needed for MMPBSA) ---
    gas_input = TLEAP_GAS_TEMPLATE.format(
        receptor_ff=receptor_ff,
        ligand_ff=ligand_ff,
        ligand_frcmod=str(Path(ligand_frcmod).resolve()),
        ligand_mol2=str(Path(ligand_mol2).resolve()),
        receptor_pdb=str(sanitized_pdb),
        complex_prmtop=str(output_dir / "complex.prmtop"),
        complex_inpcrd=str(output_dir / "complex.inpcrd"),
        receptor_prmtop=str(output_dir / "receptor.prmtop"),
        receptor_inpcrd=str(output_dir / "receptor.inpcrd"),
        ligand_prmtop=str(output_dir / "ligand.prmtop"),
        ligand_inpcrd=str(output_dir / "ligand.inpcrd"),
    )

    tleap_in = output_dir / "tleap_gas.in"
    tleap_in.write_text(gas_input)

    logger.info("  tleap: building gas-phase topologies")
    proc = subprocess.run(
        ["tleap", "-f", str(tleap_in)],
        capture_output=True, text=True, timeout=120,
        cwd=str(output_dir),
    )

    tleap_log = output_dir / "tleap_gas.log"
    tleap_log.write_text(proc.stdout + "\n" + proc.stderr)

    # Validate gas-phase outputs
    gas_files = {
        "complex_prmtop": output_dir / "complex.prmtop",
        "complex_inpcrd": output_dir / "complex.inpcrd",
        "receptor_prmtop": output_dir / "receptor.prmtop",
        "receptor_inpcrd": output_dir / "receptor.inpcrd",
        "ligand_prmtop": output_dir / "ligand.prmtop",
        "ligand_inpcrd": output_dir / "ligand.inpcrd",
    }

    missing = [k for k, v in gas_files.items() if not v.exists()]
    if missing:
        log_text = proc.stdout + proc.stderr
        errors = [l for l in log_text.split("\n")
                  if "ERROR" in l.upper() or "FATAL" in l.upper()]
        error_msg = "; ".join(errors[:3]) if errors else "Unknown tleap error"
        return {"success": False,
                "error": f"tleap gas-phase failed. Missing: {missing}. {error_msg}",
                "tleap_log": str(tleap_log)}

    n_atoms_complex = _parse_tleap_atom_count(proc.stdout, "COM")
    logger.info(f"  Gas-phase complex: {n_atoms_complex} atoms")

    result = {
        "success": True,
        **{k: str(v) for k, v in gas_files.items()},
        "n_atoms_complex": n_atoms_complex,
        "tleap_gas_log": str(tleap_log),
    }

    # --- Solvated topology (only for MD mode) ---
    if solvate:
        extra_ions = f"addionsrand COM Na+ 10\naddionsrand COM Cl- 10"

        solv_input = TLEAP_SOLVATED_TEMPLATE.format(
            receptor_ff=receptor_ff,
            ligand_ff=ligand_ff,
            ligand_frcmod=str(Path(ligand_frcmod).resolve()),
            ligand_mol2=str(Path(ligand_mol2).resolve()),
            receptor_pdb=str(sanitized_pdb),
            solvent_padding=solvent_padding,
            extra_ions=extra_ions,
            solvated_prmtop=str(output_dir / "solvated.prmtop"),
            solvated_inpcrd=str(output_dir / "solvated.inpcrd"),
        )

        tleap_solv_in = output_dir / "tleap_solvated.in"
        tleap_solv_in.write_text(solv_input)

        logger.info("  tleap: building solvated topology")
        proc_solv = subprocess.run(
            ["tleap", "-f", str(tleap_solv_in)],
            capture_output=True, text=True, timeout=120,
            cwd=str(output_dir),
        )

        tleap_solv_log = output_dir / "tleap_solvated.log"
        tleap_solv_log.write_text(proc_solv.stdout + "\n" + proc_solv.stderr)

        solv_prmtop = output_dir / "solvated.prmtop"
        solv_inpcrd = output_dir / "solvated.inpcrd"

        if not solv_prmtop.exists() or not solv_inpcrd.exists():
            log_text = proc_solv.stdout + proc_solv.stderr
            errors = [l for l in log_text.split("\n")
                      if "ERROR" in l.upper() or "FATAL" in l.upper()]
            error_msg = "; ".join(errors[:3]) if errors else "Unknown tleap error"
            return {"success": False,
                    "error": f"tleap solvation failed: {error_msg}",
                    "tleap_log": str(tleap_solv_log)}

        n_atoms_solvated = _parse_tleap_atom_count(proc_solv.stdout, "COM")
        logger.info(f"  Solvated complex: {n_atoms_solvated} atoms")

        result["solvated_prmtop"] = str(solv_prmtop)
        result["solvated_inpcrd"] = str(solv_inpcrd)
        result["n_atoms_solvated"] = n_atoms_solvated
        result["tleap_solvated_log"] = str(tleap_solv_log)

    return result


def _parse_tleap_atom_count(log_text: str, unit_name: str = "COM") -> int:
    """Parse atom count from tleap log output."""
    for line in log_text.split("\n"):
        if "atoms" in line.lower() and unit_name.lower() in line.lower():
            match = re.search(r"(\d+)\s+atoms", line)
            if match:
                return int(match.group(1))
    return 0


# =============================================================================
# 4a. SINGLE-POINT TRAJECTORY
# =============================================================================

def generate_single_frame_trajectory(
        prmtop: Union[str, Path],
        inpcrd: Union[str, Path],
        output_mdcrd: Union[str, Path],
) -> Dict[str, Any]:
    """
    Convert inpcrd → 1-frame mdcrd for MMPBSA.py single-point.
    """
    output_dir = Path(output_mdcrd).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    cpptraj_script = CPPTRAJ_SINGLE_FRAME.format(
        prmtop=str(Path(prmtop).resolve()),
        inpcrd=str(Path(inpcrd).resolve()),
        mdcrd=str(Path(output_mdcrd).resolve()),
    )

    script_path = output_dir / "cpptraj_single.in"
    script_path.write_text(cpptraj_script)

    logger.info("  cpptraj: generating 1-frame trajectory")
    proc = subprocess.run(
        ["cpptraj", "-i", str(script_path)],
        capture_output=True, text=True, timeout=60,
    )

    if not Path(output_mdcrd).exists():
        return {"success": False,
                "error": f"cpptraj failed: {proc.stderr[:300]}"}

    return {"success": True, "mdcrd": str(output_mdcrd), "n_frames": 1}


# =============================================================================
# 4b. MOLECULAR DYNAMICS (OpenMM)
# =============================================================================

def run_md_openmm(
        solvated_prmtop: Union[str, Path],
        solvated_inpcrd: Union[str, Path],
        output_dir: Union[str, Path],
        production_ns: float = 5.0,
        dt: float = 0.002,
        temperature: float = 300.0,
        pressure: float = 1.0,
        save_interval_ps: float = 10.0,
        restraint_k: float = 10.0,
        restrained_steps: int = 5000,
        unrestrained_steps: int = 5000,
        heating_ps: float = 100.0,
        equilibration_ps: float = 500.0,
) -> Dict[str, Any]:
    """
    Run explicit-solvent MD with OpenMM on GPU.

    Protocol:
        1. Minimization (restrained heavy atoms → unrestrained)
        2. Heating (0→T, NVT, 100 ps)
        3. Equilibration (NPT, 500 ps)
        4. Production (NPT, configurable ns)

    Args:
        solvated_prmtop: Solvated AMBER topology
        solvated_inpcrd: Solvated AMBER coordinates
        output_dir:      Directory for trajectory + logs
        production_ns:   Production MD length in nanoseconds
        dt:              Timestep in picoseconds
        temperature:     Temperature in Kelvin
        pressure:        Pressure in bar
        save_interval_ps: Frame save interval in ps
        restraint_k:     Restraint force constant (kcal/mol/Å²)
        restrained_steps:  Steps for restrained minimization
        unrestrained_steps: Steps for unrestrained minimization
        heating_ps:      Heating time in ps
        equilibration_ps: Equilibration time in ps

    Returns:
        Dict with success, trajectory path, n_frames, timings
    """
    try:
        import openmm
        import openmm.app as app
        import openmm.unit as unit
    except ImportError:
        return {"success": False,
                "error": "OpenMM not installed. Run: conda install -c conda-forge openmm"}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    trajectory_dcd = output_dir / "production.dcd"
    state_log = output_dir / "production.log"

    t_total_start = time.time()

    logger.info("  OpenMM MD Protocol")
    logger.info(f"    Production:    {production_ns} ns")
    logger.info(f"    Timestep:      {dt} ps")
    logger.info(f"    Temperature:   {temperature} K")
    logger.info(f"    Save interval: {save_interval_ps} ps")

    # --- Load system ---
    logger.info("  Loading solvated system...")
    prmtop = app.AmberPrmtopFile(str(solvated_prmtop))
    inpcrd = app.AmberInpcrdFile(str(solvated_inpcrd))

    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=10.0 * unit.angstroms,
        constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,
    )

    n_atoms = system.getNumParticles()
    logger.info(f"    System: {n_atoms} atoms")

    # --- 1. MINIMIZATION (restrained) ---
    logger.info("  [1/5] Restrained minimization...")
    t0 = time.time()

    restraint_force = openmm.CustomExternalForce(
        "k*periodicdistance(x,y,z,x0,y0,z0)^2"
    )
    restraint_force.addGlobalParameter(
        "k", restraint_k * unit.kilocalories_per_mole / unit.angstroms ** 2
    )
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    positions = inpcrd.positions
    for i in range(n_atoms):
        mass = system.getParticleMass(i)
        if mass > 2.0 * unit.amu:
            restraint_force.addParticle(
                i, positions[i].value_in_unit(unit.nanometers)
            )

    restraint_idx = system.addForce(restraint_force)

    integrator = openmm.LangevinMiddleIntegrator(
        temperature * unit.kelvin, 1.0 / unit.picosecond, dt * unit.picoseconds
    )

    platform = openmm.Platform.getPlatformByName("CUDA")
    platform_properties = {"Precision": "mixed"}
    simulation = app.Simulation(
        prmtop.topology, system, integrator, platform, platform_properties
    )
    simulation.context.setPositions(positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    simulation.minimizeEnergy(maxIterations=restrained_steps)
    logger.info(f"    Restrained min: {time.time() - t0:.1f}s")

    # --- 2. MINIMIZATION (unrestrained) ---
    logger.info("  [2/5] Unrestrained minimization...")
    t0 = time.time()
    system.removeForce(restraint_idx)
    simulation.context.reinitialize(preserveState=True)
    simulation.minimizeEnergy(maxIterations=unrestrained_steps)
    logger.info(f"    Unrestrained min: {time.time() - t0:.1f}s")

    # --- 3. HEATING (NVT, 0→T) ---
    logger.info(f"  [3/5] Heating 0→{temperature}K ({heating_ps} ps, NVT)...")
    t0 = time.time()

    heating_steps = int(heating_ps / dt)
    n_heating_stages = 10
    steps_per_stage = heating_steps // n_heating_stages

    for stage in range(n_heating_stages):
        frac = (stage + 1) / n_heating_stages
        target_temp = frac * temperature
        integrator.setTemperature(target_temp * unit.kelvin)
        simulation.step(steps_per_stage)

    logger.info(f"    Heating: {time.time() - t0:.1f}s")

    # --- 4. EQUILIBRATION (NPT) ---
    logger.info(f"  [4/5] Equilibration ({equilibration_ps} ps, NPT)...")
    t0 = time.time()

    barostat = openmm.MonteCarloBarostat(
        pressure * unit.bar, temperature * unit.kelvin
    )
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)
    integrator.setTemperature(temperature * unit.kelvin)

    eq_steps = int(equilibration_ps / dt)
    simulation.step(eq_steps)
    logger.info(f"    Equilibration: {time.time() - t0:.1f}s")

    # --- 5. PRODUCTION (NPT) ---
    logger.info(f"  [5/5] Production ({production_ns} ns, NPT)...")
    t0 = time.time()

    prod_steps = int(production_ns * 1000 / dt)
    save_interval_steps = int(save_interval_ps / dt)

    simulation.reporters.append(
        app.DCDReporter(str(trajectory_dcd), save_interval_steps)
    )
    simulation.reporters.append(
        app.StateDataReporter(
            str(state_log), save_interval_steps * 10,
            step=True, time=True, potentialEnergy=True,
            temperature=True, speed=True,
        )
    )

    simulation.step(prod_steps)
    prod_time = time.time() - t0
    n_frames = prod_steps // save_interval_steps
    logger.info(f"    Production: {prod_time:.1f}s ({n_frames} frames)")

    total_time = time.time() - t_total_start
    logger.info(f"  Total MD time: {total_time:.0f}s ({total_time / 60:.1f} min)")

    return {
        "success": True,
        "trajectory_dcd": str(trajectory_dcd),
        "state_log": str(state_log),
        "n_frames": n_frames,
        "production_ns": production_ns,
        "total_time_sec": round(total_time, 1),
    }


def strip_solvent_from_trajectory(
        solvated_prmtop: Union[str, Path],
        trajectory: Union[str, Path],
        output_dir: Union[str, Path],
        gas_prmtop: Union[str, Path],
        mmpbsa_interval_ps: float = 100.0,
        save_interval_ps: float = 10.0,
) -> Dict[str, Any]:
    """
    Strip water + ions from MD trajectory for MMPBSA analysis.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    dry_mdcrd = output_dir / "dry_trajectory.mdcrd"
    offset = max(1, int(mmpbsa_interval_ps / save_interval_ps))

    cpptraj_script = CPPTRAJ_STRIP_TEMPLATE.format(
        solvated_prmtop=str(Path(solvated_prmtop).resolve()),
        trajectory=str(Path(trajectory).resolve()),
        first=1,
        last="last",
        offset=offset,
        dry_mdcrd=str(dry_mdcrd.resolve()),
    )

    script_path = output_dir / "cpptraj_strip.in"
    script_path.write_text(cpptraj_script)

    logger.info(f"  cpptraj: stripping solvent (every {offset} frames)")
    proc = subprocess.run(
        ["cpptraj", "-i", str(script_path)],
        capture_output=True, text=True, timeout=300,
    )

    if not dry_mdcrd.exists():
        return {"success": False,
                "error": f"cpptraj strip failed: {proc.stderr[:300]}"}

    n_frames = 0
    for line in proc.stdout.split("\n"):
        if "frames" in line.lower() and "writing" in line.lower():
            match = re.search(r"(\d+)\s+frames", line)
            if match:
                n_frames = int(match.group(1))

    logger.info(f"  Dry trajectory: {n_frames} frames for MMPBSA")

    return {
        "success": True,
        "dry_mdcrd": str(dry_mdcrd),
        "n_frames": n_frames,
        "frame_offset": offset,
    }


# =============================================================================
# 5. MMPBSA EXECUTION
# =============================================================================

def run_mmpbsa(
        complex_prmtop: Union[str, Path],
        receptor_prmtop: Union[str, Path],
        ligand_prmtop: Union[str, Path],
        trajectory: Union[str, Path],
        output_dir: Union[str, Path],
        n_frames: int = 1,
        idecomp: int = 1,
        igb: int = 2,
        saltcon: float = 0.15,
        solvated_prmtop: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Run MMPBSA.py with per-residue decomposition.

    Args:
        complex_prmtop:  Gas-phase complex topology
        receptor_prmtop: Gas-phase receptor topology
        ligand_prmtop:   Gas-phase ligand topology
        trajectory:      mdcrd trajectory (1-frame or stripped MD)
        output_dir:      Output directory
        n_frames:        Number of frames in trajectory
        idecomp:         1=per-residue, 3=pairwise
        igb:             GB model (2=OBC-II, 5=GBn)
        saltcon:         Salt concentration in M
        solvated_prmtop: Solvated topology (optional, for -sp flag)

    Returns:
        Dict with success, output files, runtime
    """
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    mmpbsa_in = output_dir / "mmpbsa.in"
    mmpbsa_input = MMPBSA_INPUT_TEMPLATE.format(
        startframe=1,
        endframe=n_frames,
        interval=1,
        igb=igb,
        saltcon=saltcon,
        idecomp=idecomp,
    )
    mmpbsa_in.write_text(mmpbsa_input)

    cmd = [
        "MMPBSA.py",
        "-i", str(mmpbsa_in),
        "-o", str(output_dir / "FINAL_RESULTS_MMPBSA.dat"),
        "-do", str(output_dir / "FINAL_DECOMP_MMPBSA.dat"),
        "-cp", str(Path(complex_prmtop).resolve()),
        "-rp", str(Path(receptor_prmtop).resolve()),
        "-lp", str(Path(ligand_prmtop).resolve()),
        "-y", str(Path(trajectory).resolve()),
    ]

    if solvated_prmtop and Path(solvated_prmtop).exists():
        cmd.extend(["-sp", str(Path(solvated_prmtop).resolve())])

    logger.info("  MMPBSA.py: running per-residue decomposition")
    logger.info(f"    idecomp={idecomp}, igb={igb}, saltcon={saltcon}")
    logger.info(f"    Frames: {n_frames}")
    logger.debug(f"    cmd: {' '.join(cmd)}")

    t0 = time.time()
    proc = subprocess.run(
        cmd, capture_output=True, text=True,
        timeout=3600,
        cwd=str(output_dir),
    )
    runtime = time.time() - t0

    (output_dir / "mmpbsa_stdout.log").write_text(proc.stdout)
    (output_dir / "mmpbsa_stderr.log").write_text(proc.stderr)

    decomp_file = output_dir / "FINAL_DECOMP_MMPBSA.dat"
    results_file = output_dir / "FINAL_RESULTS_MMPBSA.dat"

    if proc.returncode != 0:
        error_lines = [l for l in proc.stderr.split("\n")
                       if "error" in l.lower() or "fatal" in l.lower()]
        error_msg = "; ".join(error_lines[:3]) if error_lines else proc.stderr[-300:]
        return {"success": False,
                "error": f"MMPBSA.py failed ({runtime:.1f}s): {error_msg}",
                "runtime_sec": round(runtime, 1)}

    if not decomp_file.exists():
        return {"success": False,
                "error": "MMPBSA.py completed but FINAL_DECOMP_MMPBSA.dat not found",
                "runtime_sec": round(runtime, 1)}

    logger.info(f"  MMPBSA.py completed in {runtime:.1f}s")

    return {
        "success": True,
        "decomp_file": str(decomp_file),
        "results_file": str(results_file),
        "runtime_sec": round(runtime, 1),
    }


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_mmpbsa_decomp(
        scored_mol2: Union[str, Path],
        receptor_mol2: Union[str, Path],
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        mode: str = "single_point",
        pose_selection: str = "best_score",
        pose_index: int = 1,
        receptor_ff: str = "ff14SB",
        ligand_ff: str = "gaff2",
        charge_method: str = "bcc",
        mmpbsa_idecomp: int = 1,
        mmpbsa_igb: int = 2,
        mmpbsa_saltcon: float = 0.15,
        md_params: Optional[Dict] = None,
        antechamber_timeout: int = 300,
) -> Dict[str, Any]:
    """
    Run complete MMPBSA preparation + execution pipeline.

    Steps 1-5 only. Analysis (parse + comparison) is in 01h mmpbsa_analysis.py.

    Args:
        scored_mol2:      Scored mol2 from 01c (multi-pose)
        receptor_mol2:    Receptor mol2 from 00b (for reference)
        receptor_pdb:     Receptor PDB from 00b (for tleap)
        output_dir:       Output directory
        mode:             "single_point" or "md"
        pose_selection:   "best_score" or "pose_index"
        pose_index:       1-based index (if pose_selection="pose_index")
        receptor_ff:      Protein force field
        ligand_ff:        Ligand force field
        charge_method:    "bcc" or "gas"
        mmpbsa_idecomp:   1=per-residue, 3=pairwise
        mmpbsa_igb:       GB model (2=OBC-II)
        mmpbsa_saltcon:   Salt concentration in M
        md_params:        MD parameters dict (if mode="md")
        antechamber_timeout: Timeout for antechamber

    Returns:
        Dict with success, all output paths, timings
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t_pipeline_start = time.time()
    pipeline_log = {"steps": []}

    logger.info("=" * 60)
    logger.info("  MMPBSA Per-Residue Decomposition (01g)")
    logger.info("=" * 60)
    logger.info(f"  Mode:        {mode}")
    logger.info(f"  Pose:        {pose_selection}")
    logger.info(f"  Receptor FF: {receptor_ff}")
    logger.info(f"  Ligand FF:   {ligand_ff}")
    logger.info(f"  MMPBSA:      idecomp={mmpbsa_idecomp}, igb={mmpbsa_igb}")
    logger.info(f"  Salt:        {mmpbsa_saltcon} M")

    # ─────────────────────────────────────────────────────────────
    # Step 1: Extract best pose
    # ─────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("─── Step 1: Pose Extraction ───")

    pose_mol2 = output_dir / "pose_selected.mol2"
    result_pose = extract_pose_from_mol2(
        scored_mol2=scored_mol2,
        output_mol2=pose_mol2,
        selection=pose_selection,
        pose_index=pose_index,
    )
    if not result_pose["success"]:
        return {"success": False, "error": f"Pose extraction: {result_pose['error']}"}

    pipeline_log["steps"].append({"step": "pose_extraction", **result_pose})

    # ─────────────────────────────────────────────────────────────
    # Step 2: Parametrize ligand
    # ─────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("─── Step 2: Ligand Parametrization ───")

    param_dir = output_dir / "ligand_params"
    result_param = parametrize_ligand(
        ligand_mol2=pose_mol2,
        output_dir=param_dir,
        charge_method=charge_method,
        ligand_ff=ligand_ff,
        timeout=antechamber_timeout,
    )
    if not result_param["success"]:
        return {"success": False, "error": f"Parametrization: {result_param['error']}"}

    pipeline_log["steps"].append({"step": "ligand_param",
                                  **{k: v for k, v in result_param.items() if k != "success"}})

    # ─────────────────────────────────────────────────────────────
    # Step 3: Build topologies
    # ─────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("─── Step 3: Topology Construction ───")

    topo_dir = output_dir / "topologies"
    do_solvate = (mode == "md")
    md_p = md_params or {}

    result_topo = build_topologies(
        receptor_pdb=receptor_pdb,
        ligand_mol2=result_param["mol2_path"],
        ligand_frcmod=result_param["frcmod_path"],
        output_dir=topo_dir,
        receptor_ff=receptor_ff,
        ligand_ff=ligand_ff,
        solvate=do_solvate,
        solvent_padding=md_p.get("solvent_padding", 10.0),
        ionic_strength=md_p.get("ionic_strength", 0.15),
    )
    if not result_topo["success"]:
        return {"success": False, "error": f"Topology: {result_topo['error']}"}

    pipeline_log["steps"].append({"step": "topologies",
                                  "n_atoms": result_topo.get("n_atoms_complex")})

    # ─────────────────────────────────────────────────────────────
    # Step 4: Generate trajectory
    # ─────────────────────────────────────────────────────────────
    logger.info("")

    solvated_prmtop = None

    if mode == "single_point":
        logger.info("─── Step 4: Single-Frame Trajectory ───")

        traj_mdcrd = output_dir / "trajectory" / "complex.mdcrd"
        result_traj = generate_single_frame_trajectory(
            prmtop=result_topo["complex_prmtop"],
            inpcrd=result_topo["complex_inpcrd"],
            output_mdcrd=traj_mdcrd,
        )
        if not result_traj["success"]:
            return {"success": False, "error": f"Trajectory: {result_traj['error']}"}

        n_frames = 1
        trajectory_path = result_traj["mdcrd"]

    elif mode == "md":
        logger.info("─── Step 4: Molecular Dynamics (OpenMM) ───")

        md_dir = output_dir / "md"
        min_params = md_p.get("minimization", {})

        result_md = run_md_openmm(
            solvated_prmtop=result_topo["solvated_prmtop"],
            solvated_inpcrd=result_topo["solvated_inpcrd"],
            output_dir=md_dir,
            production_ns=md_p.get("production_ns", 5.0),
            dt=md_p.get("dt", 0.002),
            temperature=md_p.get("temperature", 300.0),
            pressure=md_p.get("pressure", 1.0),
            save_interval_ps=md_p.get("save_interval_ps", 10.0),
            restraint_k=min_params.get("restraint_k", 10.0),
            restrained_steps=min_params.get("restrained_steps", 5000),
            unrestrained_steps=min_params.get("unrestrained_steps", 5000),
            heating_ps=md_p.get("heating_ps", 100.0),
            equilibration_ps=md_p.get("equilibration_ps", 500.0),
        )
        if not result_md["success"]:
            return {"success": False, "error": f"MD: {result_md['error']}"}

        pipeline_log["steps"].append({"step": "md",
                                      **{k: v for k, v in result_md.items() if k != "success"}})

        # Strip solvent
        logger.info("")
        logger.info("─── Step 4b: Strip Solvent ───")

        strip_dir = output_dir / "trajectory"
        result_strip = strip_solvent_from_trajectory(
            solvated_prmtop=result_topo["solvated_prmtop"],
            trajectory=result_md["trajectory_dcd"],
            output_dir=strip_dir,
            gas_prmtop=result_topo["complex_prmtop"],
            mmpbsa_interval_ps=md_p.get("mmpbsa_interval_ps", 100.0),
            save_interval_ps=md_p.get("save_interval_ps", 10.0),
        )
        if not result_strip["success"]:
            return {"success": False, "error": f"Strip: {result_strip['error']}"}

        n_frames = result_strip["n_frames"]
        trajectory_path = result_strip["dry_mdcrd"]
        solvated_prmtop = result_topo.get("solvated_prmtop")

    else:
        return {"success": False, "error": f"Unknown mode: {mode}"}

    pipeline_log["steps"].append({"step": "trajectory", "n_frames": n_frames, "mode": mode})

    # ─────────────────────────────────────────────────────────────
    # Step 5: MMPBSA.py
    # ─────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("─── Step 5: MMPBSA Per-Residue Decomposition ───")

    mmpbsa_dir = output_dir / "mmpbsa"
    result_mmpbsa = run_mmpbsa(
        complex_prmtop=result_topo["complex_prmtop"],
        receptor_prmtop=result_topo["receptor_prmtop"],
        ligand_prmtop=result_topo["ligand_prmtop"],
        trajectory=trajectory_path,
        output_dir=mmpbsa_dir,
        n_frames=n_frames,
        idecomp=mmpbsa_idecomp,
        igb=mmpbsa_igb,
        saltcon=mmpbsa_saltcon,
        solvated_prmtop=solvated_prmtop,
    )
    if not result_mmpbsa["success"]:
        return {"success": False, "error": f"MMPBSA: {result_mmpbsa['error']}"}

    pipeline_log["steps"].append({
        "step": "mmpbsa",
        "runtime_sec": result_mmpbsa["runtime_sec"],
    })

    # ─────────────────────────────────────────────────────────────
    # Summary
    # ─────────────────────────────────────────────────────────────
    total_time = time.time() - t_pipeline_start

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  MMPBSA Execution Complete")
    logger.info(f"  Mode:     {mode}")
    logger.info(f"  Frames:   {n_frames}")
    logger.info(f"  Time:     {total_time:.0f}s ({total_time / 60:.1f} min)")
    logger.info(f"  Output:   {output_dir}")
    logger.info(f"  Decomp:   {result_mmpbsa['decomp_file']}")
    logger.info("")
    logger.info(f"  Next: python 02_scripts/01h_mmpbsa_analysis.py \\")
    logger.info(f"          --config 03_configs/01h_mmpbsa_analysis.yaml \\")
    logger.info(f"          --campaign <campaign_config.yaml>")
    logger.info("=" * 60)

    # Save pipeline log
    log_path = output_dir / "01g_pipeline_log.json"
    pipeline_log["total_time_sec"] = round(total_time, 1)
    pipeline_log["mode"] = mode
    pipeline_log["n_frames"] = n_frames
    pipeline_log["timestamp"] = datetime.now().isoformat()
    with open(log_path, "w") as f:
        json.dump(pipeline_log, f, indent=2, default=str)

    return {
        "success": True,
        "mode": mode,
        "n_frames": n_frames,
        "decomp_file": result_mmpbsa["decomp_file"],
        "results_file": result_mmpbsa["results_file"],
        "complex_prmtop": result_topo["complex_prmtop"],
        "receptor_prmtop": result_topo["receptor_prmtop"],
        "ligand_prmtop": result_topo["ligand_prmtop"],
        "trajectory": trajectory_path,
        "pose_grid_score": result_pose.get("grid_score"),
        "total_time_sec": round(total_time, 1),
        "pipeline_log": str(log_path),
    }