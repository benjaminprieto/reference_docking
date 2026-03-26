"""
Ligand Preparation - Core Module (00a)
========================================
Prepares crystallographic ligand for DOCK6 reference scoring.

The key constraint: crystal coordinates must be preserved EXACTLY.
Antechamber recenters to origin (destroys crystal geometry).
This module provides safe alternatives:

Strategy A (default): Use reference mol2 directly
  - Verify coords match PDB HETATM records
  - Already has H + Gasteiger charges from ChimeraX/obabel
  - Fastest, safest

Strategy B: Extract from PDB + protonate
  - ChimeraX extracts HETATM as mol2 (preserves coords)
  - obabel -p adds H at target pH (preserves heavy atom coords)
  - Gasteiger charges from obabel

Strategy C: Inject AM1-BCC charges into crystal coords
  - Run antechamber on a COPY (gets AM1-BCC charges)
  - Extract charge column from antechamber output
  - Inject charges back into crystal mol2 (coords untouched)
  - Most accurate charges, crystal geometry preserved

Input:
  - Reference mol2 (campaign_dir/reference/*.mol2)
  - Receptor PDB (for HETATM extraction in Strategy B)

Output:
  - 05_results/{campaign}/00a_ligand_preparation/
    - {name}.mol2           (DOCK6-ready, crystal coords)
    - coord_validation.txt  (PDB vs mol2 comparison)

Location: 01_src/reference_docking/m00_preparation/ligand_preparation.py
Project: reference_docking
Module: 00a (core)
Version: 1.0 (2026-03-25)
"""

import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# COORDINATE EXTRACTION
# =============================================================================

def extract_mol2_heavy_coords(mol2_path: str) -> List[Tuple[str, float, float, float]]:
    """Extract heavy atom names + coords from mol2 (skip H)."""
    atoms = []
    in_atom = False
    with open(mol2_path) as f:
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    name = parts[1]
                    atom_type = parts[5] if len(parts) > 5 else ""
                    # Skip hydrogens
                    if name.startswith("H") or atom_type.startswith("H"):
                        continue
                    try:
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        atoms.append((name, x, y, z))
                    except ValueError:
                        continue
    return atoms


def extract_pdb_hetatm_coords(pdb_path: str, ligand_name: str) -> List[Tuple[str, float, float, float]]:
    """Extract heavy atom names + coords from PDB HETATM records."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM") and ligand_name in line[17:20]:
                name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else ""
                # Skip hydrogens
                if element == "H" or (not element and name.startswith("H")):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms.append((name, x, y, z))
                except (ValueError, IndexError):
                    continue
    return atoms


def validate_coordinates(
        mol2_path: str,
        pdb_path: str,
        ligand_name: str,
        tolerance: float = 0.01,
) -> Dict[str, Any]:
    """
    Validate that mol2 heavy atom coordinates match PDB HETATM records.

    Returns dict with: match (bool), n_atoms, max_deviation, details
    """
    mol2_atoms = extract_mol2_heavy_coords(mol2_path)
    pdb_atoms = extract_pdb_hetatm_coords(pdb_path, ligand_name)

    if not mol2_atoms:
        return {"match": False, "error": "No heavy atoms in mol2"}
    if not pdb_atoms:
        return {"match": False, "error": f"No HETATM records for {ligand_name} in PDB"}

    # Match by atom name
    pdb_dict = {name: (x, y, z) for name, x, y, z in pdb_atoms}
    deviations = []
    matched = 0
    unmatched = []

    for name, mx, my, mz in mol2_atoms:
        if name in pdb_dict:
            px, py, pz = pdb_dict[name]
            dev = np.sqrt((mx - px) ** 2 + (my - py) ** 2 + (mz - pz) ** 2)
            deviations.append((name, dev))
            matched += 1
        else:
            unmatched.append(name)

    max_dev = max(d for _, d in deviations) if deviations else 999.0
    mean_dev = np.mean([d for _, d in deviations]) if deviations else 999.0

    return {
        "match": max_dev <= tolerance,
        "n_mol2_heavy": len(mol2_atoms),
        "n_pdb_heavy": len(pdb_atoms),
        "n_matched": matched,
        "n_unmatched": len(unmatched),
        "max_deviation_A": round(float(max_dev), 4),
        "mean_deviation_A": round(float(mean_dev), 4),
        "unmatched_atoms": unmatched,
        "tolerance": tolerance,
    }


# =============================================================================
# STRATEGY A: USE MOL2 DIRECTLY
# =============================================================================

def prepare_from_mol2(
        reference_mol2: str,
        output_mol2: str,
        pdb_path: Optional[str] = None,
        ligand_name: str = "UDX",
) -> Dict[str, Any]:
    """
    Use reference mol2 directly (already has H + charges).
    Optionally validate coords against PDB.
    """
    ref = Path(reference_mol2)
    out = Path(output_mol2)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not ref.exists():
        return {"success": False, "error": f"Reference mol2 not found: {ref}"}

    # Copy to output
    shutil.copy2(ref, out)

    result = {
        "success": True,
        "strategy": "direct_mol2",
        "source": str(ref),
        "output": str(out),
    }

    # Validate against PDB if available
    if pdb_path and Path(pdb_path).exists():
        validation = validate_coordinates(str(ref), pdb_path, ligand_name)
        result["coord_validation"] = validation
        if validation["match"]:
            logger.info(f"  ✓ Coordinates match PDB (max dev: {validation['max_deviation_A']}Å)")
        else:
            logger.warning(f"  ⚠ Coordinate mismatch! Max deviation: {validation['max_deviation_A']}Å")

    return result


# =============================================================================
# STRATEGY B: EXTRACT FROM PDB + PROTONATE
# =============================================================================

def prepare_from_pdb(
        pdb_path: str,
        ligand_name: str,
        output_mol2: str,
        docking_ph: float = 6.3,
        chimerax_bin: str = "chimerax-daily",
) -> Dict[str, Any]:
    """
    Extract ligand from PDB with ChimeraX, protonate with obabel.

    1. ChimeraX: extract HETATM → mol2 (crystal coords exact)
    2. obabel: add H at pH (preserves heavy atom coords)
    """
    out = Path(output_mol2)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Step 1: Extract with ChimeraX
    raw_mol2 = out.parent / f"{ligand_name}_raw.mol2"
    cmd = (
        f"open {pdb_path}; "
        f"sel :{ligand_name}; "
        f"save {raw_mol2} format mol2 selectedOnly true; "
        f"exit"
    )

    try:
        proc = subprocess.run(
            [chimerax_bin, "--nogui", "--cmd", cmd],
            capture_output=True, text=True, timeout=120,
        )
        if not raw_mol2.exists():
            return {"success": False, "error": f"ChimeraX failed to extract {ligand_name}",
                    "stderr": proc.stderr[:500]}
    except FileNotFoundError:
        return {"success": False, "error": f"ChimeraX not found: {chimerax_bin}"}

    # Step 2: Protonate with obabel
    try:
        proc = subprocess.run(
            ["obabel", str(raw_mol2), "-O", str(out), "-p", str(docking_ph)],
            capture_output=True, text=True, timeout=60,
        )
        if not out.exists():
            # Fallback: use raw mol2 without protonation
            shutil.copy2(raw_mol2, out)
            logger.warning("  obabel protonation failed, using ChimeraX output as-is")
    except FileNotFoundError:
        shutil.copy2(raw_mol2, out)
        logger.warning("  obabel not found, using ChimeraX output as-is")

    # Validate
    validation = validate_coordinates(str(out), pdb_path, ligand_name)

    return {
        "success": True,
        "strategy": "pdb_extraction",
        "source": pdb_path,
        "ligand_name": ligand_name,
        "output": str(out),
        "coord_validation": validation,
    }


# =============================================================================
# STRATEGY C: INJECT AM1-BCC CHARGES
# =============================================================================

def inject_charges(
        crystal_mol2: str,
        charge_source_mol2: str,
        output_mol2: str,
) -> Dict[str, Any]:
    """
    Copy charges from antechamber mol2 into crystal mol2.

    The charge_source_mol2 has AM1-BCC charges but wrong coordinates.
    The crystal_mol2 has correct coordinates but Gasteiger charges.
    This function takes coords from crystal + charges from antechamber.
    """
    out = Path(output_mol2)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Read charges from antechamber output
    charges = {}
    in_atom = False
    with open(charge_source_mol2) as f:
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    # Map by atom index (1-based)
                    idx = int(parts[0])
                    charge = float(parts[8])
                    charges[idx] = charge

    if not charges:
        return {"success": False, "error": "No charges found in source mol2"}

    # Read crystal mol2 and replace charges
    lines = Path(crystal_mol2).read_text().split("\n")
    new_lines = []
    in_atom = False
    n_replaced = 0

    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            in_atom = True
            new_lines.append(line)
            continue
        if line.startswith("@<TRIPOS>") and in_atom:
            in_atom = False
            new_lines.append(line)
            continue

        if in_atom and line.strip():
            parts = line.split()
            if len(parts) >= 9:
                idx = int(parts[0])
                if idx in charges:
                    # Replace charge (last field in ATOM line)
                    old_charge = parts[8]
                    new_charge = f"{charges[idx]:.4f}"
                    line = line.rsplit(old_charge, 1)[0] + new_charge
                    n_replaced += 1

        new_lines.append(line)

    out.write_text("\n".join(new_lines))

    return {
        "success": True,
        "strategy": "charge_injection",
        "n_charges_replaced": n_replaced,
        "n_charges_available": len(charges),
        "output": str(out),
    }


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_ligand_preparation(
        reference_mol2: Optional[str] = None,
        pdb_path: Optional[str] = None,
        ligand_name: str = "UDX",
        output_dir: Union[str, Path] = ".",
        strategy: str = "direct",
        docking_ph: float = 6.3,
        chimerax_bin: str = "chimerax-daily",
) -> Dict[str, Any]:
    """
    Prepare crystallographic ligand for DOCK6 reference scoring.

    Args:
        reference_mol2: Path to reference mol2 (Strategy A)
        pdb_path:       Path to receptor PDB with ligand HETATM
        ligand_name:    Ligand residue name in PDB (e.g., "UDX")
        output_dir:     Output directory
        strategy:       "direct" (A), "extract" (B), or "inject" (C)
        docking_ph:     pH for protonation
        chimerax_bin:   Path to ChimeraX binary

    Returns:
        Dict with success, strategy, output path, validation
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_mol2 = str(output_dir / f"{ligand_name}.mol2")

    logger.info("=" * 60)
    logger.info("  Ligand Preparation (Reference)")
    logger.info("=" * 60)
    logger.info(f"  Strategy:    {strategy}")
    logger.info(f"  Ligand:      {ligand_name}")
    logger.info(f"  pH:          {docking_ph}")

    if strategy == "direct":
        if not reference_mol2 or not Path(reference_mol2).exists():
            return {"success": False, "error": f"Reference mol2 not found: {reference_mol2}"}
        logger.info(f"  Source:      {Path(reference_mol2).name}")
        result = prepare_from_mol2(reference_mol2, output_mol2, pdb_path, ligand_name)

    elif strategy == "extract":
        if not pdb_path or not Path(pdb_path).exists():
            return {"success": False, "error": f"PDB not found: {pdb_path}"}
        logger.info(f"  Source:      {Path(pdb_path).name}")
        result = prepare_from_pdb(pdb_path, ligand_name, output_mol2, docking_ph, chimerax_bin)

    elif strategy == "inject":
        if not reference_mol2 or not Path(reference_mol2).exists():
            return {"success": False, "error": f"Reference mol2 not found: {reference_mol2}"}
        logger.info(f"  Source:      {Path(reference_mol2).name} + AM1-BCC injection")
        # TODO: Run antechamber on copy, then inject charges
        result = {"success": False, "error": "inject strategy not yet implemented"}

    else:
        return {"success": False, "error": f"Unknown strategy: {strategy}"}

    # Save validation report
    if result.get("success") and result.get("coord_validation"):
        val = result["coord_validation"]
        report_path = output_dir / "coord_validation.txt"
        with open(report_path, "w") as f:
            f.write(f"Coordinate Validation: {ligand_name}\n")
            f.write(f"{'=' * 40}\n")
            f.write(f"Match: {val['match']}\n")
            f.write(f"Heavy atoms (mol2): {val['n_mol2_heavy']}\n")
            f.write(f"Heavy atoms (PDB):  {val['n_pdb_heavy']}\n")
            f.write(f"Matched:            {val['n_matched']}\n")
            f.write(f"Max deviation:      {val['max_deviation_A']} Å\n")
            f.write(f"Mean deviation:     {val['mean_deviation_A']} Å\n")
            f.write(f"Tolerance:          {val['tolerance']} Å\n")
        logger.info(f"  Saved: {report_path}")

    if result.get("success"):
        logger.info(f"  Output: {result['output']}")
    else:
        logger.error(f"  Failed: {result.get('error')}")

    return result
