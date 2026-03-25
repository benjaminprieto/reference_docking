"""
PLIP Interaction Analysis (03a) v1.0
=======================================
Analyze protein-ligand interactions from crystallographic coordinates
using PLIP (Protein-Ligand Interaction Profiler).

Takes the protonated receptor (from 00b) and protonated ligand (from 00c)
at the correct docking pH, combines them into a complex, and runs PLIP
to detect H-bonds, salt bridges, pi-stacking, hydrophobic contacts, etc.

Pipeline:
    1. Load receptor_protonated.pdb (from 00b, correct pH)
    2. Load ligand protonated SDF/mol2 (from 00c, correct pH)
    3. Combine → complex PDB (receptor ATOM + ligand HETATM)
    4. Run PLIP analysis
    5. Parse → standardized JSON + CSV

Input:
    Required: receptor protonated PDB + ligand protonated SDF/mol2
Output:
    interactions.json      — Standardized interaction data
    interaction_summary.csv
    complex_for_plip.pdb   — Combined PDB (intermediate)
    03a_plip_analysis.log

Location: 01_src/reference_docking/m03_crystal_analysis/plip_interaction_analysis.py
Project: reference_docking
Module: 03a (core)
Version: 1.0

Dependencies: plip (conda install -c conda-forge plip)
"""

import csv
import json
import logging
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# DATA STRUCTURES
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class Interaction:
    """One protein-ligand interaction detected by PLIP."""
    interaction_type: str       # hbond, salt_bridge, pi_stack, hydrophobic, etc.
    residue: str                # receptor residue (e.g. "SER575")
    residue_number: int = 0
    chain: str = "A"
    receptor_atom: str = ""     # PDB atom name (e.g. "OG", "NE1")
    receptor_coords: List[float] = field(default_factory=list)

    ligand_atom: str = ""       # ligand atom name or type
    ligand_coords: List[float] = field(default_factory=list)

    distance: float = 0.0       # interaction distance (Angstrom)
    angle: Optional[float] = None

    # H-bond specific
    ligand_is_donor: Optional[bool] = None

    # Salt bridge specific
    ligand_charge: Optional[str] = None  # "negative" or "positive"

    # Pi-stacking specific
    stack_type: Optional[str] = None     # "P" (parallel) or "T" (T-shaped)

    # Extra metadata
    extra: Dict[str, Any] = field(default_factory=dict)


# ═══════════════════════════════════════════════════════════════════════
# COMPLEX PDB GENERATION
# ═══════════════════════════════════════════════════════════════════════

def _ligand_to_pdb(input_path: str, output_pdb: str,
                   ligand_name: str = "LIG") -> bool:
    """Convert ligand SDF/mol2 to PDB via OpenBabel."""
    cmd = ["obabel", str(input_path), "-O", output_pdb]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if result.returncode == 0 and Path(output_pdb).exists():
            return True
        logger.warning(f"  obabel failed: {result.stderr[:200]}")
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        logger.warning(f"  obabel not available: {e}")

    # Fallback: try RDKit
    try:
        from rdkit import Chem
        path = Path(input_path)
        mol = None
        if path.suffix.lower() in (".sdf", ".mol"):
            supplier = Chem.SDMolSupplier(str(path), removeHs=False)
            for m in supplier:
                if m is not None:
                    mol = m
                    break
        elif path.suffix.lower() == ".mol2":
            mol = Chem.MolFromMol2File(str(path), removeHs=False)

        if mol is not None:
            Chem.MolToPDBFile(mol, output_pdb)
            return Path(output_pdb).exists()
    except Exception as e:
        logger.warning(f"  RDKit PDB conversion failed: {e}")

    return False


def create_complex_pdb(receptor_pdb: str, ligand_path: str,
                       output_pdb: str, ligand_name: str = "LIG") -> bool:
    """
    Combine protonated receptor PDB + protonated ligand into one PDB.

    Receptor atoms → ATOM records (as-is)
    Ligand atoms → HETATM records
    """
    rec_path = Path(receptor_pdb)
    if not rec_path.exists():
        logger.error(f"  Receptor not found: {receptor_pdb}")
        return False

    # Convert ligand to PDB
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        lig_pdb = tmp.name

    if not _ligand_to_pdb(ligand_path, lig_pdb, ligand_name):
        logger.error(f"  Could not convert ligand to PDB: {ligand_path}")
        return False

    # Read receptor lines
    rec_text = rec_path.read_text()
    lig_text = Path(lig_pdb).read_text()

    # Combine
    with open(output_pdb, "w") as out:
        # Write receptor (skip END/ENDMDL)
        for line in rec_text.split("\n"):
            if line.startswith(("END", "ENDMDL")):
                continue
            out.write(line + "\n")

        # Write ligand as HETATM
        for line in lig_text.split("\n"):
            if line.startswith(("ATOM", "HETATM")):
                # Force HETATM
                if line.startswith("ATOM"):
                    line = "HETATM" + line[6:]
                out.write(line + "\n")

        out.write("END\n")

    logger.info(f"  Complex PDB: {Path(output_pdb).name}")
    return True


# ═══════════════════════════════════════════════════════════════════════
# PLIP ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def _safe_coords(atom) -> List[float]:
    """Safely extract coordinates from a PLIP atom."""
    try:
        if hasattr(atom, 'coords'):
            c = atom.coords
            return [round(float(c[0]), 3), round(float(c[1]), 3),
                    round(float(c[2]), 3)]
    except (TypeError, IndexError):
        pass
    return []


def _safe_atom_name(atom) -> str:
    """Safely extract atom name/type from a PLIP atom."""
    try:
        if hasattr(atom, 'type'):
            return str(atom.type).strip()
        if hasattr(atom, 'OBAtom'):
            res = atom.OBAtom.GetResidue()
            if res:
                return res.GetAtomID(atom.OBAtom).strip()
    except Exception:
        pass
    return ""


def run_plip(complex_pdb: str) -> List[Interaction]:
    """
    Run PLIP on a protein-ligand complex PDB.

    Returns list of Interaction objects.
    """
    try:
        from plip.structure.preparation import PDBComplex
    except ImportError:
        logger.error("  PLIP not installed!")
        logger.error("  Install: conda install -c conda-forge plip")
        return []

    logger.info(f"  Running PLIP on {Path(complex_pdb).name}")

    mol = PDBComplex()
    mol.load_pdb(complex_pdb)
    mol.analyze()

    interactions = []

    for bsid, site in sorted(mol.interaction_sets.items()):
        lig_id = bsid[0] if isinstance(bsid, tuple) else str(bsid)
        chain = bsid[1] if isinstance(bsid, tuple) and len(bsid) > 1 else "?"
        logger.info(f"  Binding site: {lig_id} (chain {chain})")

        # ── H-bonds: protein is donor ──
        for hb in site.hbonds_pdon:
            try:
                interactions.append(Interaction(
                    interaction_type="hbond",
                    residue=f"{hb.restype}{hb.resnr}",
                    residue_number=hb.resnr,
                    chain=hb.reschain,
                    receptor_atom=_safe_atom_name(hb.d) if hasattr(hb, 'd') else "",
                    receptor_coords=_safe_coords(hb.d) if hasattr(hb, 'd') else [],
                    ligand_atom=_safe_atom_name(hb.a) if hasattr(hb, 'a') else "",
                    ligand_coords=_safe_coords(hb.a) if hasattr(hb, 'a') else [],
                    distance=round(hb.distance_ad, 2) if hasattr(hb, 'distance_ad') else 0.0,
                    angle=round(hb.angle, 1) if hasattr(hb, 'angle') else None,
                    ligand_is_donor=False,
                ))
            except Exception as e:
                logger.debug(f"    Skipped hbond_pdon: {e}")

        # ── H-bonds: ligand is donor ──
        for hb in site.hbonds_ldon:
            try:
                interactions.append(Interaction(
                    interaction_type="hbond",
                    residue=f"{hb.restype}{hb.resnr}",
                    residue_number=hb.resnr,
                    chain=hb.reschain,
                    receptor_atom=_safe_atom_name(hb.a) if hasattr(hb, 'a') else "",
                    receptor_coords=_safe_coords(hb.a) if hasattr(hb, 'a') else [],
                    ligand_atom=_safe_atom_name(hb.d) if hasattr(hb, 'd') else "",
                    ligand_coords=_safe_coords(hb.d) if hasattr(hb, 'd') else [],
                    distance=round(hb.distance_ad, 2) if hasattr(hb, 'distance_ad') else 0.0,
                    angle=round(hb.angle, 1) if hasattr(hb, 'angle') else None,
                    ligand_is_donor=True,
                ))
            except Exception as e:
                logger.debug(f"    Skipped hbond_ldon: {e}")

        # ── Salt bridges: ligand negative ──
        for sb in site.saltbridge_lneg:
            try:
                # Ligand negative atoms centroid
                neg_atoms = sb.negative.atoms if hasattr(sb, 'negative') else []
                neg_names = ",".join(_safe_atom_name(a) for a in neg_atoms)
                neg_coords = []
                if neg_atoms:
                    coords_list = [_safe_coords(a) for a in neg_atoms if _safe_coords(a)]
                    if coords_list:
                        neg_coords = np.mean(coords_list, axis=0).round(3).tolist()

                # Receptor positive atoms centroid
                pos_atoms = sb.positive.atoms if hasattr(sb, 'positive') else []
                pos_coords = []
                if pos_atoms:
                    coords_list = [_safe_coords(a) for a in pos_atoms if _safe_coords(a)]
                    if coords_list:
                        pos_coords = np.mean(coords_list, axis=0).round(3).tolist()

                interactions.append(Interaction(
                    interaction_type="salt_bridge",
                    residue=f"{sb.restype}{sb.resnr}",
                    residue_number=sb.resnr,
                    chain=sb.reschain,
                    receptor_atom="Positive",
                    receptor_coords=pos_coords,
                    ligand_atom=neg_names or "Negative",
                    ligand_coords=neg_coords,
                    distance=round(sb.distance, 2),
                    ligand_charge="negative",
                ))
            except Exception as e:
                logger.debug(f"    Skipped saltbridge_lneg: {e}")

        # ── Salt bridges: ligand positive ──
        for sb in site.saltbridge_pneg:
            try:
                pos_atoms = sb.positive.atoms if hasattr(sb, 'positive') else []
                pos_names = ",".join(_safe_atom_name(a) for a in pos_atoms)
                pos_coords = []
                if pos_atoms:
                    coords_list = [_safe_coords(a) for a in pos_atoms if _safe_coords(a)]
                    if coords_list:
                        pos_coords = np.mean(coords_list, axis=0).round(3).tolist()

                neg_atoms = sb.negative.atoms if hasattr(sb, 'negative') else []
                neg_coords = []
                if neg_atoms:
                    coords_list = [_safe_coords(a) for a in neg_atoms if _safe_coords(a)]
                    if coords_list:
                        neg_coords = np.mean(coords_list, axis=0).round(3).tolist()

                interactions.append(Interaction(
                    interaction_type="salt_bridge",
                    residue=f"{sb.restype}{sb.resnr}",
                    residue_number=sb.resnr,
                    chain=sb.reschain,
                    receptor_atom="Negative",
                    receptor_coords=neg_coords,
                    ligand_atom=pos_names or "Positive",
                    ligand_coords=pos_coords,
                    distance=round(sb.distance, 2),
                    ligand_charge="positive",
                ))
            except Exception as e:
                logger.debug(f"    Skipped saltbridge_pneg: {e}")

        # ── Pi-stacking ──
        for ps in site.pistacking:
            try:
                lig_center = list(ps.ligandring.center) if hasattr(ps, 'ligandring') else []
                rec_center = list(ps.proteinring.center) if hasattr(ps, 'proteinring') else []

                interactions.append(Interaction(
                    interaction_type="pi_stack",
                    residue=f"{ps.restype}{ps.resnr}",
                    residue_number=ps.resnr,
                    chain=ps.reschain,
                    receptor_atom="Ring",
                    receptor_coords=[round(c, 3) for c in rec_center] if rec_center else [],
                    ligand_atom="Ring",
                    ligand_coords=[round(c, 3) for c in lig_center] if lig_center else [],
                    distance=round(ps.distance, 2),
                    angle=round(ps.angle, 1) if hasattr(ps, 'angle') else None,
                    stack_type=ps.type if hasattr(ps, 'type') else None,
                ))
            except Exception as e:
                logger.debug(f"    Skipped pistacking: {e}")

        # ── Hydrophobic contacts ──
        for hc in site.hydrophobic_contacts:
            try:
                interactions.append(Interaction(
                    interaction_type="hydrophobic",
                    residue=f"{hc.restype}{hc.resnr}",
                    residue_number=hc.resnr,
                    chain=hc.reschain,
                    receptor_atom=_safe_atom_name(hc.bsatom) if hasattr(hc, 'bsatom') else "",
                    receptor_coords=_safe_coords(hc.bsatom) if hasattr(hc, 'bsatom') else [],
                    ligand_atom=_safe_atom_name(hc.ligatom) if hasattr(hc, 'ligatom') else "",
                    ligand_coords=_safe_coords(hc.ligatom) if hasattr(hc, 'ligatom') else [],
                    distance=round(hc.distance, 2),
                ))
            except Exception as e:
                logger.debug(f"    Skipped hydrophobic: {e}")

        # ── Water bridges ──
        for wb in site.water_bridges:
            try:
                interactions.append(Interaction(
                    interaction_type="water_bridge",
                    residue=f"{wb.restype}{wb.resnr}",
                    residue_number=wb.resnr,
                    chain=wb.reschain,
                    receptor_atom=_safe_atom_name(wb.a) if hasattr(wb, 'a') else "",
                    receptor_coords=_safe_coords(wb.a) if hasattr(wb, 'a') else [],
                    ligand_atom=_safe_atom_name(wb.d) if hasattr(wb, 'd') else "",
                    ligand_coords=_safe_coords(wb.d) if hasattr(wb, 'd') else [],
                    distance=round(wb.distance_aw, 2) if hasattr(wb, 'distance_aw') else 0.0,
                    ligand_is_donor=not wb.protisdon if hasattr(wb, 'protisdon') else None,
                    extra={"water_coords": _safe_coords(wb.water) if hasattr(wb, 'water') else []},
                ))
            except Exception as e:
                logger.debug(f"    Skipped water_bridge: {e}")

        # ── Pi-cation ──
        for pc in getattr(site, 'pication_laro', []):
            try:
                interactions.append(Interaction(
                    interaction_type="pi_cation",
                    residue=f"{pc.restype}{pc.resnr}",
                    residue_number=pc.resnr,
                    chain=pc.reschain,
                    receptor_atom="Positive",
                    ligand_atom="Ring",
                    distance=round(pc.distance, 2),
                ))
            except Exception as e:
                logger.debug(f"    Skipped pication: {e}")

        for pc in getattr(site, 'pication_paro', []):
            try:
                interactions.append(Interaction(
                    interaction_type="pi_cation",
                    residue=f"{pc.restype}{pc.resnr}",
                    residue_number=pc.resnr,
                    chain=pc.reschain,
                    receptor_atom="Ring",
                    ligand_atom="Positive",
                    distance=round(pc.distance, 2),
                ))
            except Exception as e:
                logger.debug(f"    Skipped pication: {e}")

    # Summary
    type_counts = Counter(i.interaction_type for i in interactions)
    logger.info(f"  PLIP detected {len(interactions)} interactions:")
    for itype, n in sorted(type_counts.items()):
        logger.info(f"    {itype}: {n}")

    return interactions


# ═══════════════════════════════════════════════════════════════════════
# OUTPUT: JSON
# ═══════════════════════════════════════════════════════════════════════

def _interaction_to_dict(inter: Interaction) -> Dict[str, Any]:
    """Convert Interaction to JSON-serializable dict."""
    d = {
        "interaction_type": inter.interaction_type,
        "residue": inter.residue,
        "residue_number": inter.residue_number,
        "chain": inter.chain,
        "receptor_atom": inter.receptor_atom,
        "receptor_coords": inter.receptor_coords,
        "ligand_atom": inter.ligand_atom,
        "ligand_coords": inter.ligand_coords,
        "distance": inter.distance,
    }

    if inter.angle is not None:
        d["angle"] = inter.angle
    if inter.ligand_is_donor is not None:
        d["ligand_is_donor"] = inter.ligand_is_donor
    if inter.ligand_charge is not None:
        d["ligand_charge"] = inter.ligand_charge
    if inter.stack_type is not None:
        d["stack_type"] = inter.stack_type
    if inter.extra:
        d["extra"] = inter.extra

    return d


def write_interactions_json(interactions: List[Interaction],
                            ligand_name: str,
                            receptor_name: str,
                            output_path: Path) -> None:
    """Write standardized interactions JSON."""
    type_counts = Counter(i.interaction_type for i in interactions)

    data = {
        "version": "3.0",
        "date": datetime.now().isoformat(),
        "ligand_name": ligand_name,
        "receptor_name": receptor_name,
        "n_interactions": len(interactions),
        "counts": dict(type_counts),
        "unique_residues": sorted(set(i.residue for i in interactions)),
        "n_unique_residues": len(set(i.residue for i in interactions)),
        "coordinate_convention": {
            "ligand_coords": "Ligand atom coordinates (use for pharmacophore features)",
            "receptor_coords": "Receptor atom coordinates (for reference)",
        },
        "interactions": [_interaction_to_dict(i) for i in interactions],
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    logger.info(f"  Saved: {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# OUTPUT: CSV SUMMARY
# ═══════════════════════════════════════════════════════════════════════

def write_interaction_summary_csv(interactions: List[Interaction],
                                  output_path: Path) -> None:
    """Write interaction summary CSV."""
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Type", "Residue", "Chain", "Receptor_Atom",
            "Ligand_Atom", "Distance", "Angle",
            "Ligand_Is_Donor", "Ligand_Charge", "Stack_Type",
            "Ligand_x", "Ligand_y", "Ligand_z",
            "Receptor_x", "Receptor_y", "Receptor_z",
        ])

        for inter in interactions:
            lc = inter.ligand_coords
            rc = inter.receptor_coords
            writer.writerow([
                inter.interaction_type,
                inter.residue,
                inter.chain,
                inter.receptor_atom,
                inter.ligand_atom,
                inter.distance,
                inter.angle if inter.angle is not None else "",
                inter.ligand_is_donor if inter.ligand_is_donor is not None else "",
                inter.ligand_charge if inter.ligand_charge is not None else "",
                inter.stack_type if inter.stack_type is not None else "",
                lc[0] if len(lc) >= 3 else "",
                lc[1] if len(lc) >= 3 else "",
                lc[2] if len(lc) >= 3 else "",
                rc[0] if len(rc) >= 3 else "",
                rc[1] if len(rc) >= 3 else "",
                rc[2] if len(rc) >= 3 else "",
            ])

    logger.info(f"  Saved: {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# OUTPUT: TXT SUMMARY
# ═══════════════════════════════════════════════════════════════════════

def write_summary_txt(interactions: List[Interaction],
                      ligand_name: str, receptor_name: str,
                      output_path: Path) -> None:
    w = 78
    type_counts = Counter(i.interaction_type for i in interactions)
    unique_res = sorted(set(i.residue for i in interactions))

    lines = [
        "=" * w,
        "03a PLIP INTERACTION ANALYSIS — SUMMARY (v1.0)",
        "=" * w,
        "",
        f"Date:       {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Ligand:     {ligand_name}",
        f"Receptor:   {receptor_name}",
        f"Interactions: {len(interactions)}",
        f"Residues:   {len(unique_res)}",
        "",
    ]

    lines.append("-" * w)
    lines.append("Interaction counts:")
    lines.append("-" * w)
    for itype, n in sorted(type_counts.items()):
        lines.append(f"  {itype:<20s}: {n}")

    lines.extend(["", "-" * w, "Residues involved:", "-" * w])
    for res in unique_res:
        res_inters = [i for i in interactions if i.residue == res]
        types_str = ", ".join(sorted(set(i.interaction_type for i in res_inters)))
        lines.append(f"  {res:<12s}: {types_str}")

    lines.extend(["", "-" * w, "Interaction details:", "-" * w])
    for i, inter in enumerate(interactions):
        angle_str = f"  angle={inter.angle:.0f}°" if inter.angle else ""
        donor_str = ""
        if inter.ligand_is_donor is not None:
            donor_str = " (lig=donor)" if inter.ligand_is_donor else " (lig=acceptor)"
        if inter.ligand_charge:
            donor_str = f" (lig={inter.ligand_charge})"

        lines.append(
            f"  [{i+1:2d}] {inter.interaction_type:<15s} "
            f"{inter.ligand_atom:<10s} → "
            f"{inter.residue}:{inter.receptor_atom:<6s} "
            f"@ {inter.distance:.1f}Å{angle_str}{donor_str}"
        )

    lines.extend(["", "=" * w])
    output_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info(f"  Saved: {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════

def run_plip_analysis(
        receptor_pdb: str,
        ligand_path: str,
        output_dir: str,
        ligand_name: Optional[str] = None,
        output_name: str = "interactions",
) -> Dict[str, Any]:
    """
    Full pipeline: receptor + ligand → PLIP analysis → JSON + CSV.

    Args:
        receptor_pdb:  Protonated receptor PDB (from 00b)
        ligand_path:   Protonated ligand SDF/mol2 (from 00c)
        output_dir:    Output directory
        ligand_name:   Name for the ligand (default: filename stem)
        output_name:   Base filename for outputs

    Returns:
        dict with success, paths, interactions summary
    """
    logger.info("=" * 60)
    logger.info("PLIP INTERACTION ANALYSIS (03a) v1.0")
    logger.info("=" * 60)

    rec_file = Path(receptor_pdb)
    lig_file = Path(ligand_path)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if ligand_name is None:
        ligand_name = lig_file.stem

    logger.info(f"  Receptor: {rec_file.name}")
    logger.info(f"  Ligand:   {lig_file.name} ({ligand_name})")

    # ── Step 1: Create complex PDB ──
    logger.info("")
    logger.info("Step 1: Creating complex PDB")
    complex_pdb = str(out_dir / "complex_for_plip.pdb")

    if not create_complex_pdb(str(rec_file), str(lig_file), complex_pdb):
        return {"success": False, "error": "Failed to create complex PDB"}

    # ── Step 2: Run PLIP ──
    logger.info("")
    logger.info("Step 2: PLIP analysis")
    interactions = run_plip(complex_pdb)

    if not interactions:
        logger.warning("  PLIP returned no interactions!")
        logger.warning("  Check: is the ligand correctly placed in the complex?")

    # ── Step 3: Write outputs ──
    logger.info("")
    logger.info("Step 3: Writing outputs")

    json_path = out_dir / f"{output_name}.json"
    write_interactions_json(interactions, ligand_name, rec_file.name, json_path)

    csv_path = out_dir / f"{output_name}_summary.csv"
    write_interaction_summary_csv(interactions, csv_path)

    txt_path = out_dir / f"{output_name}_summary.txt"
    write_summary_txt(interactions, ligand_name, rec_file.name, txt_path)

    # Summary
    type_counts = Counter(i.interaction_type for i in interactions)
    unique_res = sorted(set(i.residue for i in interactions))

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  DONE — {len(interactions)} interactions, "
                f"{len(unique_res)} residues")
    logger.info(f"  JSON: {json_path}")
    logger.info("=" * 60)

    for inter in interactions:
        donor = ""
        if inter.ligand_is_donor is not None:
            donor = " (donor)" if inter.ligand_is_donor else " (acceptor)"
        if inter.ligand_charge:
            donor = f" ({inter.ligand_charge})"
        logger.info(f"  {inter.interaction_type:<15s} "
                     f"{inter.ligand_atom} → {inter.residue}:{inter.receptor_atom} "
                     f"@ {inter.distance:.1f}Å{donor}")

    return {
        "success": True,
        "json_path": str(json_path),
        "csv_path": str(csv_path),
        "txt_path": str(txt_path),
        "complex_pdb": complex_pdb,
        "n_interactions": len(interactions),
        "counts": dict(type_counts),
        "unique_residues": unique_res,
    }