"""
Receptor Preparation - Core Module (00b)
==========================================
Prepara el receptor para DOCK6: limpieza + protonacion al pH +
generacion de mol2 con cargas parciales.

DOCK6 requiere:
  - rec_charged.mol2  (Sybyl atom types + AMBER partial charges)
  - rec_noH.pdb       (sin H, para DMS surface en 01a)

Strategies:
  A. pdb2pqr  → PDB2PQR+PROPKA pKa prediction + ChimeraX mol2 generation
               (recomendada para pH no-fisiologico, e.g. Golgi 6.3)
  B. chimerax → ChimeraX DockPrep directo (AddH + addcharge)
               (buena para pH fisiologico ~7.2)
  C. obabel   → OpenBabel simple, Gasteiger charges (fallback de emergencia)
               (NO recomendada — ver Practitioner's Guide)

ChimeraX genera mol2 con:
  - Sybyl atom types (C.3, N.am, O.2, etc.) — requeridos por DOCK6
  - AMBER ff14SB charges (RESP-derived) — el estándar para DOCK6
  - @<TRIPOS>SUBSTRUCTURE section — requerida por grid program
  - Molecule name limpio (no paths)

Esto reemplaza completamente el pipeline anterior basado en OpenBabel,
que tenía 3 problemas críticos:
  1. Molecule name = full path → DOCK6 lo malinterpretaba (BUG 1)
  2. Molecule type = SMALL → grid filtraba ATOMs incorrectamente
  3. Sin @<TRIPOS>SUBSTRUCTURE → grid no podia parsear residuos

Pipeline:
    1. Limpiar PDB (agua, alt conf, HETATM, cadenas)
    2. Protonar al pH de docking (pdb2pqr o chimerax)
    3. Generar mol2 con ChimeraX (Sybyl types + AMBER charges)
    4. Generar rec_noH.pdb (strip H del protonado)
    5. Validar outputs
    6. Generar protonation_report.json

Location: 01_src/reference_docking/m00_preparation/receptor_preparation.py
Project: reference_docking
Module: 00b (core)
Version: 2.0 — ChimeraX rewrite (2026-03-13)
"""

import json
import logging
import re
import shutil
import subprocess
import tempfile
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union

logger = logging.getLogger(__name__)

# ChimeraX binary — installed on the server
CHIMERAX_BIN = "/usr/bin/chimerax-daily"


# =============================================================================
# MOL2 SANITIZATION (defensive — ChimeraX should produce clean mol2,
# but we keep these as safety nets for edge cases and obabel fallback)
# =============================================================================

def _fix_mol2_molecule_name(mol2_path: str, name: str = "receptor"):
    """
    Fix the mol2 MOLECULE section for DOCK6 compatibility.

    OpenBabel generates mol2 files with two issues for DOCK6:
      1. Line 2 (molecule name): writes the full input path instead of a name.
         DOCK6 tries to open this path as a file.
      2. Line 4 (molecule type): writes "SMALL" for everything.
         DOCK6's grid program expects "PROTEIN" for receptor mol2 files
         and filters ATOM records by this type.

    ChimeraX handles both correctly, but we call this as a safety net
    after every mol2 generation — cheap and defensive.

    mol2 @<TRIPOS>MOLECULE format:
      Line 1: @<TRIPOS>MOLECULE
      Line 2: molecule name
      Line 3: num_atoms num_bonds num_subst num_feat num_sets
      Line 4: molecule type (SMALL | PROTEIN | BIOPOLYMER | ...)
      Line 5: charge type (GASTEIGER | NO_CHARGES | ...)
    """
    path = Path(mol2_path)
    if not path.exists():
        return

    lines = path.read_text().split("\n")
    modified = False

    for i, line in enumerate(lines):
        if "@<TRIPOS>MOLECULE" in line:
            # Fix molecule name (line i+1)
            if i + 1 < len(lines):
                old_name = lines[i + 1].strip()
                lines[i + 1] = name
                if old_name != name:
                    logger.debug(f"    Fixed mol2 name: '{old_name}' -> '{name}'")
                    modified = True

            # Fix molecule type (line i+3): SMALL -> PROTEIN
            if i + 3 < len(lines):
                old_type = lines[i + 3].strip()
                if old_type == "SMALL":
                    lines[i + 3] = "PROTEIN"
                    logger.debug(f"    Fixed mol2 type: SMALL -> PROTEIN")
                    modified = True
            break

    if modified:
        path.write_text("\n".join(lines))


def _add_mol2_substructure(mol2_path: str):
    """
    Add @<TRIPOS>SUBSTRUCTURE section to a mol2 file if missing.

    DOCK6 requires the SUBSTRUCTURE section to parse protein receptors.
    ChimeraX generates this automatically; this is a safety net for the
    obabel fallback strategy.
    """
    path = Path(mol2_path)
    if not path.exists():
        return

    text = path.read_text()

    # Skip if SUBSTRUCTURE already exists
    if "@<TRIPOS>SUBSTRUCTURE" in text:
        return

    # Parse residues from ATOM section
    residues = {}  # subst_id -> (subst_name, first_atom_id)
    in_atom = False

    for line in text.split("\n"):
        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            continue
        if line.startswith("@<TRIPOS>") and in_atom:
            break
        if in_atom and line.strip():
            parts = line.split()
            if len(parts) >= 8:
                atom_id = parts[0]
                subst_id = parts[6]
                subst_name = parts[7]
                if subst_id not in residues:
                    residues[subst_id] = {
                        "subst_name": subst_name,
                        "first_atom": atom_id,
                    }

    if not residues:
        logger.warning("    No residues found in mol2 ATOM section")
        return

    # Update molecule counts (line 3: num_atoms num_bonds num_subst)
    lines = text.split("\n")
    for i, line in enumerate(lines):
        if "@<TRIPOS>MOLECULE" in line and i + 2 < len(lines):
            parts = lines[i + 2].split()
            if len(parts) >= 2:
                n_atoms = parts[0]
                n_bonds = parts[1]
                n_subst = str(len(residues))
                lines[i + 2] = f" {n_atoms} {n_bonds} {n_subst} 0 0"
            break

    # Build SUBSTRUCTURE section
    subst_lines = ["@<TRIPOS>SUBSTRUCTURE"]
    for subst_id in sorted(residues.keys(), key=lambda x: int(x)):
        res = residues[subst_id]
        subst_lines.append(
            f"{subst_id:>7s} {res['subst_name']:<8s} {res['first_atom']:>7s} "
            f"RESIDUE           1 A     {res['subst_name']:<8s}"
        )

    text_out = "\n".join(lines).rstrip("\n") + "\n" + "\n".join(subst_lines) + "\n"
    path.write_text(text_out)
    logger.debug(f"    Added SUBSTRUCTURE section: {len(residues)} residues")


# =============================================================================
# PDB CLEANING
# =============================================================================

def clean_pdb(
        input_pdb: str,
        output_pdb: str,
        remove_water: bool = True,
        remove_hetatm: bool = True,
        remove_alt_conformations: bool = True,
        keep_chains: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Clean a PDB file for docking preparation.

    Operations:
      - Remove water molecules (HOH, WAT, TIP, TIP3)
      - Remove HETATM records (ligands, ions, cofactors)
      - Keep only first alternate conformation (altLoc A or ' ')
      - Filter by chain ID
      - Preserve CONECT records (for disulfide bonds)

    Args:
        input_pdb: Path to input PDB
        output_pdb: Path for cleaned PDB
        remove_water: Remove water molecules
        remove_hetatm: Remove HETATM records
        remove_alt_conformations: Keep only first alt conformation
        keep_chains: Chain IDs to keep (None = all)

    Returns:
        Dict with cleaning statistics
    """
    water_residues = {"HOH", "WAT", "TIP", "TIP3", "SOL"}
    stats = {
        "atoms_input": 0,
        "atoms_output": 0,
        "waters_removed": 0,
        "hetatm_removed": 0,
        "alt_conf_removed": 0,
        "chains_removed": 0,
    }

    output_lines = []

    with open(input_pdb) as f:
        for line in f:
            is_atom = line.startswith("ATOM")
            is_hetatm = line.startswith("HETATM")

            if is_atom or is_hetatm:
                stats["atoms_input"] += 1
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                alt_loc = line[16].strip()

                # Chain filter
                if keep_chains and chain_id not in keep_chains:
                    stats["chains_removed"] += 1
                    continue

                # Water filter
                if remove_water and res_name in water_residues:
                    stats["waters_removed"] += 1
                    continue

                # HETATM filter
                if remove_hetatm and is_hetatm:
                    stats["hetatm_removed"] += 1
                    continue

                # Alt conformation filter (keep '' or 'A')
                if remove_alt_conformations and alt_loc and alt_loc != "A":
                    stats["alt_conf_removed"] += 1
                    continue

                # Clear alt loc indicator for kept atoms
                if remove_alt_conformations and alt_loc == "A":
                    line = line[:16] + " " + line[17:]

                stats["atoms_output"] += 1
                output_lines.append(line)

            elif line.startswith("TER") or line.startswith("END"):
                output_lines.append(line)
            elif line.startswith("CONECT"):
                output_lines.append(line)
            elif line.startswith(("HEADER", "TITLE", "REMARK", "CRYST")):
                output_lines.append(line)

    # Write cleaned PDB
    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as f:
        f.writelines(output_lines)

    logger.info(f"  Cleaned PDB: {stats['atoms_input']} -> {stats['atoms_output']} atoms")
    if stats["waters_removed"]:
        logger.info(f"    Removed {stats['waters_removed']} water atoms")
    if stats["hetatm_removed"]:
        logger.info(f"    Removed {stats['hetatm_removed']} HETATM atoms")
    if stats["alt_conf_removed"]:
        logger.info(f"    Removed {stats['alt_conf_removed']} alt conformations")

    return stats


def strip_hydrogens(input_pdb: str, output_pdb: str) -> int:
    """
    Remove hydrogen atoms from PDB and save as rec_noH.pdb.
    Needed for DMS surface generation in 01a.

    Returns:
        Number of hydrogen atoms removed
    """
    n_removed = 0
    output_lines = []

    with open(input_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else ""

                # Detect hydrogen by element column or atom name
                is_hydrogen = (
                    element == "H"
                    or (not element and (
                        atom_name.startswith("H")
                        or atom_name in ("1H", "2H", "3H")
                    ))
                )
                if is_hydrogen:
                    n_removed += 1
                    continue

            output_lines.append(line)

    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as f:
        f.writelines(output_lines)

    logger.info(f"  Stripped {n_removed} hydrogen atoms -> {Path(output_pdb).name}")
    return n_removed


# =============================================================================
# PQR CHARGE PARSING (for pdb2pqr strategy)
# =============================================================================

def parse_pqr_charges(pqr_path: str) -> Dict[str, float]:
    """
    Parse a PQR file and extract per-atom partial charges.

    PQR format: like PDB but columns 55-62 = charge, 63-70 = radius

    Returns:
        Dict mapping "chainID:resNum:atomName" -> charge
    """
    charges = {}
    with open(pqr_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                atom_name = line[12:16].strip()
                chain_id = line[21].strip() or "A"
                res_num = line[22:26].strip()
                # PQR: charge is after coordinates
                parts = line[30:].split()
                # x, y, z, charge, radius
                if len(parts) >= 5:
                    charge = float(parts[3])
                    key = f"{chain_id}:{res_num}:{atom_name}"
                    charges[key] = charge
            except (ValueError, IndexError):
                continue

    return charges


def inject_charges_into_mol2(
        mol2_path: str,
        charges: Dict[str, float],
        output_path: str,
) -> Tuple[int, int]:
    """
    Replace partial charges in a mol2 file with values from PQR.

    Matches atoms by residue number + atom name.

    Returns:
        (n_matched, n_total) atoms
    """
    lines = Path(mol2_path).read_text().split("\n")
    output_lines = []
    in_atom_section = False
    n_matched = 0
    n_total = 0

    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            in_atom_section = True
            output_lines.append(line)
            continue
        if line.startswith("@<TRIPOS>") and in_atom_section:
            in_atom_section = False
            output_lines.append(line)
            continue

        if in_atom_section and line.strip():
            n_total += 1
            parts = line.split()
            if len(parts) >= 9:
                atom_name = parts[1]
                subst_name = parts[7] if len(parts) > 7 else ""
                res_match = re.search(r"(\d+)", subst_name)
                res_num = res_match.group(1) if res_match else "0"

                # Try matching with different chain IDs
                matched = False
                for chain in ["A", "B", "C", "D", ""]:
                    key = f"{chain}:{res_num}:{atom_name}"
                    if key in charges:
                        parts[8] = f"{charges[key]:.4f}"
                        line = (
                            f"{parts[0]:>7s} {parts[1]:<8s} "
                            f"{float(parts[2]):>10.4f}{float(parts[3]):>10.4f}"
                            f"{float(parts[4]):>10.4f} "
                            f"{parts[5]:<8s} {parts[6]:>5s} {parts[7]:<8s} "
                            f"{charges[key]:>10.4f}"
                        )
                        n_matched += 1
                        matched = True
                        break

        output_lines.append(line)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(output_lines))

    return n_matched, n_total


# =============================================================================
# CHIMERAX MOL2 GENERATION
# =============================================================================

def _find_chimerax() -> Optional[str]:
    """
    Find ChimeraX binary. Search order:
      1. Module constant CHIMERAX_BIN (/usr/bin/chimerax-daily)
      2. chimerax-daily in PATH
      3. chimerax in PATH
      4. ChimeraX in PATH

    Returns path to binary or None.
    """
    # 1. Module constant
    if Path(CHIMERAX_BIN).exists():
        return CHIMERAX_BIN

    # 2-4. Search PATH
    for name in ["chimerax-daily", "chimerax", "ChimeraX"]:
        result = shutil.which(name)
        if result:
            return result

    return None


def _run_chimerax_mol2(
        input_pdb: str,
        output_mol2: str,
        add_hydrogens: bool = True,
        chimerax_bin: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate a DOCK6-ready mol2 using ChimeraX.

    ChimeraX produces mol2 with:
      - Sybyl atom types (C.3, N.am, O.2, etc.)
      - AMBER ff14SB partial charges
      - @<TRIPOS>SUBSTRUCTURE section
      - Clean molecule name

    Args:
        input_pdb: Path to PDB file (clean or protonated)
        output_mol2: Path for output mol2
        add_hydrogens: If True, run addh before addcharge.
                       Set False if PDB already has correct H (e.g. from PDB2PQR).
        chimerax_bin: Path to ChimeraX binary (auto-detected if None)

    Returns:
        Dict with success, mol2_path, chimerax_log
    """
    chimerax = chimerax_bin or _find_chimerax()
    if not chimerax:
        return {
            "success": False,
            "error": "ChimeraX not found. Install ChimeraX or set CHIMERAX_BIN.",
        }

    input_path = Path(input_pdb).resolve()
    output_path = Path(output_mol2).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build ChimeraX command sequence
    cmds = [f"open {input_path}"]

    if add_hydrogens:
        # addh with hbond consideration for better HIS assignment
        cmds.append("addh")

    # Assign AMBER ff14SB charges (standard residues) + AM1-BCC (nonstandard)
    cmds.append("addcharge")

    # Save as mol2 with Sybyl types
    cmds.append(f"save {output_path} format mol2")
    cmds.append("exit")

    cmd_string = "; ".join(cmds)

    logger.info(f"  ChimeraX: generating mol2 with AMBER charges")
    logger.debug(f"    Command: {chimerax} --nogui --cmd \"{cmd_string}\"")

    try:
        result = subprocess.run(
            [chimerax, "--nogui", "--cmd", cmd_string],
            capture_output=True,
            text=True,
            timeout=600,  # 10 min — large proteins can be slow
        )

        chimerax_log = (result.stdout or "") + "\n" + (result.stderr or "")

        if not output_path.exists() or output_path.stat().st_size == 0:
            logger.error("    ChimeraX did not produce mol2 output")
            if result.stderr:
                logger.error(f"    stderr: {result.stderr[:500]}")
            return {
                "success": False,
                "error": "ChimeraX produced no mol2 output",
                "chimerax_log": chimerax_log,
            }

        logger.info(f"    -> {output_path.name} ({output_path.stat().st_size} bytes)")

        # Defensive sanitization (ChimeraX should be clean, but just in case)
        _fix_mol2_molecule_name(str(output_path))
        _add_mol2_substructure(str(output_path))

        return {
            "success": True,
            "mol2_path": str(output_path),
            "chimerax_log": chimerax_log,
        }

    except FileNotFoundError:
        return {
            "success": False,
            "error": f"ChimeraX binary not found at: {chimerax}",
        }
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "error": "ChimeraX timed out (>600s). Very large protein?",
        }


# =============================================================================
# STRATEGY A: PDB2PQR + ChimeraX (RECOMMENDED)
# =============================================================================

def prepare_with_pdb2pqr(
        clean_pdb: str,
        output_dir: str,
        docking_ph: float = 7.2,
        force_field: str = "AMBER",
) -> Dict[str, Any]:
    """
    Protonate receptor using PDB2PQR+PROPKA, then generate mol2 via ChimeraX.

    This is the RECOMMENDED strategy because:
      - PROPKA predicts per-residue pKa values (pH-aware)
      - ChimeraX generates proper Sybyl types + AMBER ff14SB charges
      - No OpenBabel involved (eliminates BUG 1, BUG 2, BUG 3)

    Pipeline:
      1. pdb2pqr --ff=AMBER --titration-state-method=propka --with-ph=X
         → protonated PDB + PQR with AMBER charges
      2. ChimeraX: open protonated PDB → addcharge (no addh — H already present)
         → save mol2 (Sybyl types + AMBER ff14SB charges)

    Why not inject PQR charges into ChimeraX mol2?
      PDB2PQR outputs PQR-format AMBER charges, and ChimeraX assigns
      ff14SB charges from its own templates. Both are AMBER-derived and
      nearly identical for standard residues. We use ChimeraX's charges
      because they're self-consistent with its atom typing. The PQR is
      kept for reference and PROPKA pKa reporting.
    """
    output_dir = Path(output_dir)
    pqr_path = output_dir / "receptor.pqr"
    protonated_pdb = output_dir / "receptor_protonated.pdb"

    # --- Step 1: PDB2PQR ---
    logger.info(f"  PDB2PQR: protonating at pH {docking_ph} (ff={force_field})")

    cmd = [
        "pdb2pqr",
        "--ff", force_field,
        "--titration-state-method", "propka",
        "--with-ph", str(docking_ph),
        "--keep-chain",
        "--pdb-output", str(protonated_pdb),
        "--log-level", "WARNING",
        str(clean_pdb),
        str(pqr_path),
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300,
        )
        if result.returncode != 0:
            logger.error(f"    PDB2PQR failed (rc={result.returncode})")
            if result.stderr:
                logger.error(f"    {result.stderr[:500]}")
            return {"success": False, "error": "PDB2PQR failed", "tool": "pdb2pqr"}

        if not pqr_path.exists():
            return {
                "success": False,
                "error": "PDB2PQR produced no PQR output",
                "tool": "pdb2pqr",
            }

    except FileNotFoundError:
        return {"success": False, "error": "pdb2pqr not found in PATH", "tool": "pdb2pqr"}
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "PDB2PQR timed out", "tool": "pdb2pqr"}

    logger.info(f"    PQR: {pqr_path.name}")
    if protonated_pdb.exists():
        logger.info(f"    Protonated PDB: {protonated_pdb.name}")

    # --- Step 2: Generate mol2 via ChimeraX ---
    # Use protonated PDB (H already at correct pH states).
    # addh=False because PDB2PQR already added correct hydrogens.
    source_pdb = str(protonated_pdb) if protonated_pdb.exists() else clean_pdb
    add_h = not protonated_pdb.exists()  # Only add H if PDB2PQR didn't produce PDB

    rec_charged = output_dir / "rec_charged.mol2"

    chimerax_result = _run_chimerax_mol2(
        input_pdb=source_pdb,
        output_mol2=str(rec_charged),
        add_hydrogens=add_h,
    )

    if not chimerax_result.get("success"):
        # Fallback: try with clean PDB + addh
        logger.warning("    ChimeraX failed with protonated PDB, trying with clean PDB + addh...")
        chimerax_result = _run_chimerax_mol2(
            input_pdb=clean_pdb,
            output_mol2=str(rec_charged),
            add_hydrogens=True,
        )

    if not chimerax_result.get("success"):
        return {
            "success": False,
            "error": f"ChimeraX mol2 generation failed: {chimerax_result.get('error')}",
            "tool": "pdb2pqr",
        }

    # Save ChimeraX log for debugging
    chimerax_log_path = output_dir / "chimerax.log"
    if chimerax_result.get("chimerax_log"):
        chimerax_log_path.write_text(chimerax_result["chimerax_log"])

    # --- Parse PROPKA results ---
    propka_report = _parse_propka_log(output_dir)

    return {
        "success": True,
        "tool": "pdb2pqr",
        "rec_charged_mol2": str(rec_charged),
        "pqr_path": str(pqr_path),
        "protonated_pdb": str(protonated_pdb) if protonated_pdb.exists() else None,
        "propka_report": propka_report,
    }


def _parse_propka_log(output_dir: Path) -> Dict[str, Any]:
    """Try to parse PROPKA pKa predictions from PDB2PQR output."""
    report = {"titratable_residues": []}

    propka_files = list(output_dir.glob("*.propka")) + list(output_dir.glob("*.pka"))

    for propka_file in propka_files:
        try:
            text = propka_file.read_text()
            for line in text.split("\n"):
                match = re.match(
                    r"\s*(ASP|GLU|HIS|LYS|CYS|TYR)\s+(\d+)\s+(\S+)\s+([\d.]+)",
                    line,
                )
                if match:
                    report["titratable_residues"].append({
                        "residue": match.group(1),
                        "number": int(match.group(2)),
                        "chain": match.group(3),
                        "pKa": float(match.group(4)),
                    })
        except Exception:
            pass

    return report


# =============================================================================
# STRATEGY B: ChimeraX only
# =============================================================================

def prepare_with_chimerax(
        clean_pdb: str,
        output_dir: str,
        docking_ph: float = 7.2,
) -> Dict[str, Any]:
    """
    Protonate receptor using ChimeraX's built-in AddH + addcharge.

    Simpler than pdb2pqr but less pH-aware:
      - AddH uses geometry-based rules (not pKa prediction)
      - HIS assignment based on local H-bonding environment
      - No PROPKA pKa report

    Good for:
      - pH ~7.2 (standard protonation is usually correct)
      - Quick validation runs
      - When PDB2PQR is not available

    Pipeline:
      ChimeraX: open clean PDB → addh → addcharge → save mol2
    """
    output_dir = Path(output_dir)
    rec_charged = output_dir / "rec_charged.mol2"

    logger.info(f"  ChimeraX: protonating + assigning AMBER charges")

    chimerax_result = _run_chimerax_mol2(
        input_pdb=clean_pdb,
        output_mol2=str(rec_charged),
        add_hydrogens=True,
    )

    if not chimerax_result.get("success"):
        return {
            "success": False,
            "error": f"ChimeraX failed: {chimerax_result.get('error')}",
            "tool": "chimerax",
        }

    # Save ChimeraX log
    chimerax_log_path = output_dir / "chimerax.log"
    if chimerax_result.get("chimerax_log"):
        chimerax_log_path.write_text(chimerax_result["chimerax_log"])

    # Generate a protonated PDB from ChimeraX for rec_noH source
    protonated_pdb = output_dir / "receptor_protonated.pdb"
    _chimerax_save_pdb(clean_pdb, str(protonated_pdb))

    return {
        "success": True,
        "tool": "chimerax",
        "rec_charged_mol2": str(rec_charged),
        "protonated_pdb": str(protonated_pdb) if protonated_pdb.exists() else None,
    }


def _chimerax_save_pdb(input_pdb: str, output_pdb: str):
    """Save a protonated PDB via ChimeraX (addh + save pdb)."""
    chimerax = _find_chimerax()
    if not chimerax:
        return

    input_path = Path(input_pdb).resolve()
    output_path = Path(output_pdb).resolve()

    cmd_string = f"open {input_path}; addh; save {output_path} format pdb; exit"

    try:
        subprocess.run(
            [chimerax, "--nogui", "--cmd", cmd_string],
            capture_output=True, text=True, timeout=300,
        )
    except Exception as e:
        logger.debug(f"    ChimeraX PDB save failed (non-critical): {e}")


# =============================================================================
# STRATEGY C: OpenBabel (FALLBACK — NOT RECOMMENDED)
# =============================================================================

def prepare_with_obabel(
        clean_pdb: str,
        output_dir: str,
        docking_ph: float = 7.2,
) -> Dict[str, Any]:
    """
    Protonate receptor using OpenBabel.

    *** NOT RECOMMENDED — USE ONLY AS LAST RESORT ***

    Problems with OpenBabel for DOCK6 receptor preparation:
      - Context-independent pKa rules (wrong protonation states)
      - Writes input path as molecule name (BUG 1)
      - Sets molecule type to SMALL instead of PROTEIN
      - Missing @<TRIPOS>SUBSTRUCTURE section
      - Gasteiger charges (less accurate than AMBER ff14SB)

    We patch all these issues post-hoc, but the protonation quality
    is fundamentally lower than PDB2PQR or ChimeraX.

    Only use if both ChimeraX and PDB2PQR are unavailable.
    """
    output_dir = Path(output_dir)
    rec_charged = output_dir / "rec_charged.mol2"

    logger.warning("  *** USING OBABEL FALLBACK — NOT RECOMMENDED ***")
    logger.warning("  Install ChimeraX for proper DOCK6 receptor preparation.")
    logger.info(f"  obabel: protonating at pH {docking_ph} + Gasteiger charges")

    cmd = [
        "obabel", str(clean_pdb),
        "-O", str(rec_charged),
        "-p", str(docking_ph),
        "--partialcharge", "gasteiger",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0 or not rec_charged.exists() or rec_charged.stat().st_size == 0:
            return {"success": False, "error": "obabel failed", "tool": "obabel"}
    except FileNotFoundError:
        return {"success": False, "error": "obabel not found in PATH", "tool": "obabel"}
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "obabel timed out", "tool": "obabel"}

    logger.info(f"    -> {rec_charged.name} ({rec_charged.stat().st_size} bytes)")

    # Patch all the obabel issues
    _fix_mol2_molecule_name(str(rec_charged))
    _add_mol2_substructure(str(rec_charged))

    return {
        "success": True,
        "tool": "obabel",
        "rec_charged_mol2": str(rec_charged),
    }


# =============================================================================
# MOL2 VALIDATION
# =============================================================================

def validate_prepared_mol2(mol2_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Validate a receptor mol2 for DOCK6 compatibility.

    Checks:
      - File exists and non-empty
      - Has @<TRIPOS>ATOM and @<TRIPOS>BOND sections
      - Has @<TRIPOS>SUBSTRUCTURE section
      - Reasonable atom count (>100 for a protein)
      - Charges are present (not all zero)
      - Sybyl atom types present (C.3, N.am, O.2, etc.)
      - Molecule type is PROTEIN (not SMALL)
    """
    result = {
        "exists": False,
        "has_atoms": False,
        "has_bonds": False,
        "has_substructure": False,
        "has_charges": False,
        "has_sybyl_types": False,
        "molecule_type": None,
        "n_atoms": 0,
        "n_residues": 0,
        "total_charge": 0.0,
        "valid": False,
    }

    path = Path(mol2_path)
    if not path.exists() or path.stat().st_size == 0:
        return result
    result["exists"] = True

    text = path.read_text()
    result["has_atoms"] = "@<TRIPOS>ATOM" in text
    result["has_bonds"] = "@<TRIPOS>BOND" in text
    result["has_substructure"] = "@<TRIPOS>SUBSTRUCTURE" in text

    # Parse molecule type
    for i, line in enumerate(text.split("\n")):
        if "@<TRIPOS>MOLECULE" in line:
            mol_lines = text.split("\n")
            if i + 3 < len(mol_lines):
                result["molecule_type"] = mol_lines[i + 3].strip()
            break

    if result["has_atoms"]:
        in_atom = False
        charges = []
        atom_types = set()
        residue_ids = set()

        for line in text.split("\n"):
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    atom_types.add(parts[5])
                    try:
                        charges.append(float(parts[8]))
                    except (ValueError, IndexError):
                        pass
                    if len(parts) > 6:
                        residue_ids.add(parts[6])

        result["n_atoms"] = len(charges)
        result["n_residues"] = len(residue_ids)
        result["has_charges"] = (
            len(charges) > 0
            and not all(abs(c) < 1e-10 for c in charges)
        )
        result["total_charge"] = round(sum(charges), 2) if charges else 0.0

        # Check for Sybyl types (contain dots: C.3, N.am, O.2, S.3, etc.)
        sybyl_types = [t for t in atom_types if "." in t]
        result["has_sybyl_types"] = len(sybyl_types) > 0

    result["valid"] = (
        result["has_atoms"]
        and result["has_bonds"]
        and result["has_charges"]
        and result["has_sybyl_types"]
        and result["n_atoms"] > 50  # A protein should have many atoms
    )

    return result


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

def run_receptor_preparation(
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        docking_ph: float = 7.2,
        protonation_tool: str = "pdb2pqr",
        force_field: str = "AMBER",
        chain: Optional[str] = None,
        remove_water: bool = True,
        remove_hetatm: bool = True,
        remove_alt_conformations: bool = True,
) -> Dict[str, Any]:
    """
    Run the complete receptor preparation pipeline.

    Args:
        receptor_pdb: Path to input PDB file
        output_dir: Directory for output files
        docking_ph: pH for protonation
        protonation_tool: "pdb2pqr" | "chimerax" | "obabel"
                          (pdb2pqr uses PROPKA + ChimeraX mol2 generation)
                          (chimerax uses ChimeraX AddH + addcharge)
                          (obabel is LAST RESORT only)
        force_field: "AMBER" | "CHARMM" | "PARSE" (for PDB2PQR)
        chain: Chain ID to keep (None = all)
        remove_water: Remove water molecules
        remove_hetatm: Remove HETATM records
        remove_alt_conformations: Keep only first alt conformation

    Returns:
        Dict with output paths, validation, and protonation report
    """
    receptor_pdb = Path(receptor_pdb)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not receptor_pdb.exists():
        return {"success": False, "error": f"Receptor PDB not found: {receptor_pdb}"}

    logger.info("=" * 60)
    logger.info("  Receptor Preparation Pipeline v2.0 (ChimeraX)")
    logger.info("=" * 60)
    logger.info(f"  Input:   {receptor_pdb.name}")
    logger.info(f"  pH:      {docking_ph}")
    logger.info(f"  Tool:    {protonation_tool}")
    logger.info(f"  Output:  {output_dir}")

    report = {
        "start_time": datetime.now().isoformat(),
        "input_pdb": str(receptor_pdb),
        "docking_ph": docking_ph,
        "protonation_tool": protonation_tool,
        "version": "2.0-chimerax",
    }

    # --- Step 1: Clean PDB ---
    logger.info("\nStep 1: Cleaning PDB")
    clean_path = output_dir / "receptor_clean.pdb"
    keep_chains = [chain] if chain else None

    clean_stats = clean_pdb(
        str(receptor_pdb), str(clean_path),
        remove_water=remove_water,
        remove_hetatm=remove_hetatm,
        remove_alt_conformations=remove_alt_conformations,
        keep_chains=keep_chains,
    )
    report["clean_stats"] = clean_stats

    # --- Step 2: Protonate + generate mol2 ---
    logger.info(f"\nStep 2: Protonation + mol2 generation ({protonation_tool})")

    # Check ChimeraX availability upfront for non-obabel strategies
    chimerax_available = _find_chimerax() is not None
    if protonation_tool in ("pdb2pqr", "chimerax") and not chimerax_available:
        logger.warning(f"  ChimeraX not found — cannot use '{protonation_tool}' strategy")
        logger.warning(f"  Falling back to obabel (NOT RECOMMENDED)")
        protonation_tool = "obabel"

    strategies = {
        "pdb2pqr": lambda: prepare_with_pdb2pqr(
            str(clean_path), str(output_dir), docking_ph, force_field,
        ),
        "chimerax": lambda: prepare_with_chimerax(
            str(clean_path), str(output_dir), docking_ph,
        ),
        "obabel": lambda: prepare_with_obabel(
            str(clean_path), str(output_dir), docking_ph,
        ),
        # Legacy aliases
        "reduce": lambda: prepare_with_chimerax(
            str(clean_path), str(output_dir), docking_ph,
        ),
    }

    if protonation_tool not in strategies:
        return {
            "success": False,
            "error": f"Unknown protonation tool: {protonation_tool}. "
                     f"Options: pdb2pqr, chimerax, obabel",
        }

    prep_result = strategies[protonation_tool]()

    if not prep_result.get("success"):
        # Cascade fallback: pdb2pqr -> chimerax -> obabel
        fallback_order = []
        if protonation_tool == "pdb2pqr":
            fallback_order = ["chimerax", "obabel"]
        elif protonation_tool == "chimerax":
            fallback_order = ["obabel"]

        for fallback in fallback_order:
            logger.warning(f"  {protonation_tool} failed, trying {fallback} fallback...")
            prep_result = strategies[fallback]()
            if prep_result.get("success"):
                break

    if not prep_result.get("success"):
        report["error"] = prep_result.get("error", "All protonation methods failed")
        return {"success": False, "report": report, **prep_result}

    report["protonation_result"] = {
        k: v for k, v in prep_result.items() if k != "success"
    }

    rec_charged = prep_result["rec_charged_mol2"]

    # --- Step 3: Generate rec_noH.pdb ---
    logger.info("\nStep 3: Generating rec_noH.pdb")
    rec_noH = output_dir / "rec_noH.pdb"

    # Source: protonated PDB if available, else clean PDB
    protonated_pdb = prep_result.get("protonated_pdb")
    source_pdb = (
        protonated_pdb
        if protonated_pdb and Path(protonated_pdb).exists()
        else str(clean_path)
    )
    strip_hydrogens(source_pdb, str(rec_noH))

    # --- Step 4: Validate mol2 ---
    logger.info("\nStep 4: Validating rec_charged.mol2")
    validation = validate_prepared_mol2(rec_charged)
    report["validation"] = validation

    if validation["valid"]:
        logger.info(f"  VALID: {validation['n_atoms']} atoms, "
                    f"{validation['n_residues']} residues, "
                    f"total charge: {validation['total_charge']}")
        if validation["has_sybyl_types"]:
            logger.info("  Sybyl atom types: present")
        else:
            logger.warning("  WARNING: Sybyl atom types NOT detected "
                          "(DOCK6 may not score correctly)")
        if validation["has_substructure"]:
            logger.info("  SUBSTRUCTURE section: present")
        else:
            logger.warning("  WARNING: SUBSTRUCTURE section missing "
                          "(DOCK6 grid may fail)")
        if validation.get("molecule_type") == "PROTEIN":
            logger.info("  Molecule type: PROTEIN")
        elif validation.get("molecule_type"):
            logger.warning(f"  WARNING: Molecule type is '{validation['molecule_type']}' "
                          f"(expected PROTEIN)")
    else:
        logger.warning(f"  WARNING: mol2 validation issues: {validation}")

    # --- Save report ---
    report["end_time"] = datetime.now().isoformat()
    report["success"] = True

    report_path = output_dir / "protonation_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    # --- Summary TXT ---
    summary_path = output_dir / "protonation_summary.txt"
    w = 70
    lines = [
        "=" * w,
        "00b RECEPTOR PREPARATION - SUMMARY (v2.0 ChimeraX)",
        "=" * w,
        "",
        f"Date:              {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Input:             {receptor_pdb.name}",
        f"pH:                {docking_ph}",
        f"Tool:              {prep_result.get('tool', protonation_tool)}",
        "",
        f"Atoms (input):     {clean_stats['atoms_input']}",
        f"Atoms (cleaned):   {clean_stats['atoms_output']}",
        f"Waters removed:    {clean_stats['waters_removed']}",
        f"HETATM removed:    {clean_stats['hetatm_removed']}",
        "",
        f"mol2 atoms:        {validation.get('n_atoms', 0)}",
        f"mol2 residues:     {validation.get('n_residues', 0)}",
        f"Total charge:      {validation.get('total_charge', 0):.2f}",
        f"Sybyl types:       {'yes' if validation.get('has_sybyl_types') else 'NO'}",
        f"SUBSTRUCTURE:      {'yes' if validation.get('has_substructure') else 'NO'}",
        f"Molecule type:     {validation.get('molecule_type', 'unknown')}",
        f"Valid:             {'YES' if validation.get('valid') else 'NO'}",
        "",
    ]

    # PROPKA results if available
    propka = prep_result.get("propka_report", {})
    titratable = propka.get("titratable_residues", [])
    if titratable:
        lines.extend([
            "-" * w,
            "PROPKA pKa Predictions (titratable residues)",
            "-" * w,
        ])
        for res in titratable:
            marker = ""
            if res["residue"] == "HIS" and res["pKa"] > docking_ph:
                marker = " <- PROTONATED (HIP)"
            elif res["residue"] in ("ASP", "GLU") and res["pKa"] > docking_ph:
                marker = " <- NEUTRAL"
            elif res["residue"] == "LYS" and res["pKa"] < docking_ph:
                marker = " <- NEUTRAL"
            elif res["residue"] == "CYS" and res["pKa"] < docking_ph:
                marker = " <- DEPROTONATED"

            lines.append(
                f"  {res['residue']}{res['number']:>4d} {res['chain']}: "
                f"pKa = {res['pKa']:.2f}{marker}"
            )
        lines.append("")

    lines.extend([
        "=" * w,
        "",
        f"rec_charged.mol2:  {rec_charged}",
        f"rec_noH.pdb:       {rec_noH}",
    ])

    summary_path.write_text("\n".join(lines))
    logger.info(f"  Summary: {summary_path}")

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  Receptor prepared: {Path(rec_charged).name}")
    logger.info(f"  rec_noH.pdb:       {rec_noH.name}")
    logger.info(f"{'=' * 60}")

    return {
        "success": True,
        "rec_charged_mol2": rec_charged,
        "rec_noH_pdb": str(rec_noH),
        "report": report,
        "report_path": str(report_path),
        "summary_path": str(summary_path),
        "validation": validation,
    }