"""
Binding Site Definition - Core Module (00e)
=============================================
Identifica el binding site y recorta el receptor para grid generation.

Proteinas grandes (>3000 atomos) hacen que sphgen (Fortran) falle
por overflow de enteros en la superficie DMS. Este modulo recorta
el receptor a una esfera alrededor del binding site, produciendo
un PDB reducido que DMS/sphgen pueden manejar.

Metodos de definicion del binding site:
  A. reference_ligand  -> residuos dentro de contact_cutoff del ligando
  B. residues          -> lista explicita de residuos
  C. coordinates       -> centro explicito (x, y, z)

Pipeline:
    1. Definir binding site (centroide + residuos contacto)
    2. Recortar rec_noH.pdb -> rec_noH_site.pdb (atomos dentro de trim_radius)
    3. Generar reporte (residuos, centroide, estadisticas)

Input:
  - rec_noH.pdb (de 00b, receptor sin H)
  - Referencia: ligando mol2/pdb, residuos, o coordenadas

Output:
  - rec_noH_site.pdb       (PDB recortado para DMS/sphgen en 01a)
  - binding_site_report.json
  - binding_site_summary.txt

Location: 01_src/reference_docking/m00_preparation/binding_site_definition.py
Project: reference_docking
Module: 00e (core)
Version: 1.0
"""

import json
import logging
import math
import re
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union

logger = logging.getLogger(__name__)


# =============================================================================
# MOL2 COORDINATE EXTRACTION
# =============================================================================

def read_mol2_coords(mol2_path: str) -> List[Tuple[float, float, float]]:
    """
    Extract atomic coordinates from a mol2 file.

    Returns:
        List of (x, y, z) tuples
    """
    coords = []
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
                if len(parts) >= 5:
                    try:
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        coords.append((x, y, z))
                    except (ValueError, IndexError):
                        continue

    return coords


def compute_centroid(coords: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    """Compute geometric centroid of a list of coordinates."""
    if not coords:
        raise ValueError("No coordinates provided")
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    return (cx, cy, cz)


# =============================================================================
# PDB UTILITIES
# =============================================================================

def _parse_pdb_atom(line: str) -> Optional[Dict[str, Any]]:
    """Parse an ATOM/HETATM line from a PDB file."""
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        return {
            "record": line[:6].strip(),
            "atom_name": line[12:16].strip(),
            "res_name": line[17:20].strip(),
            "chain": line[21].strip(),
            "res_num": int(line[22:26].strip()),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "line": line,
        }
    except (ValueError, IndexError):
        return None


def read_pdb_atoms(pdb_path: str, protein_only: bool = True) -> List[Dict[str, Any]]:
    """
    Read all ATOM records from a PDB file.

    Args:
        pdb_path: Path to PDB file
        protein_only: Only read ATOM records (skip HETATM)

    Returns:
        List of atom dicts
    """
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if protein_only and not line.startswith("ATOM"):
                continue
            if not protein_only and not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom = _parse_pdb_atom(line)
            if atom:
                atoms.append(atom)
    return atoms


# =============================================================================
# METHOD A: BINDING SITE FROM REFERENCE LIGAND
# =============================================================================

def find_contact_residues(
        pdb_path: str,
        ligand_coords: List[Tuple[float, float, float]],
        contact_cutoff: float = 5.0,
        protein_only: bool = True,
) -> List[Dict[str, Any]]:
    """
    Find protein residues within contact_cutoff of ligand atoms.

    Args:
        pdb_path: Path to receptor PDB
        ligand_coords: List of (x, y, z) from the reference ligand
        contact_cutoff: Distance cutoff in Angstroms
        protein_only: Only consider ATOM records (not HETATM)

    Returns:
        List of dicts with: res_name, res_num, chain, min_distance, n_contacts
    """
    atoms = read_pdb_atoms(pdb_path, protein_only=protein_only)
    cutoff_sq = contact_cutoff ** 2

    # Track contacts per residue
    residue_contacts: Dict[str, Dict[str, Any]] = {}

    for atom in atoms:
        ax, ay, az = atom["x"], atom["y"], atom["z"]
        res_key = f"{atom['chain']}:{atom['res_name']}{atom['res_num']}"

        for lx, ly, lz in ligand_coords:
            dsq = (ax - lx) ** 2 + (ay - ly) ** 2 + (az - lz) ** 2
            if dsq <= cutoff_sq:
                d = math.sqrt(dsq)
                if res_key not in residue_contacts:
                    residue_contacts[res_key] = {
                        "res_name": atom["res_name"],
                        "res_num": atom["res_num"],
                        "chain": atom["chain"],
                        "min_distance": d,
                        "n_contacts": 0,
                    }
                residue_contacts[res_key]["n_contacts"] += 1
                if d < residue_contacts[res_key]["min_distance"]:
                    residue_contacts[res_key]["min_distance"] = d
                break  # one contact per atom is enough

    result = sorted(
        residue_contacts.values(),
        key=lambda r: r["res_num"],
    )

    logger.info(f"  Contact residues (cutoff={contact_cutoff} A): {len(result)}")
    return result


def binding_site_from_ligand(
        pdb_path: str,
        ligand_mol2: str,
        contact_cutoff: float = 5.0,
) -> Dict[str, Any]:
    """
    Method A: Define binding site from reference ligand contacts.

    Returns:
        Dict with: centroid, residues, ligand_centroid, n_ligand_atoms
    """
    # Read ligand coordinates
    lig_coords = read_mol2_coords(ligand_mol2)
    if not lig_coords:
        raise ValueError(f"No coordinates found in ligand file: {ligand_mol2}")

    lig_centroid = compute_centroid(lig_coords)
    logger.info(f"  Ligand: {Path(ligand_mol2).name} ({len(lig_coords)} atoms)")
    logger.info(f"  Ligand centroid: ({lig_centroid[0]:.3f}, {lig_centroid[1]:.3f}, {lig_centroid[2]:.3f})")

    # Find contact residues
    residues = find_contact_residues(pdb_path, lig_coords, contact_cutoff)

    if not residues:
        raise ValueError(
            f"No protein residues found within {contact_cutoff} A of ligand. "
            f"Check that ligand coordinates overlap with receptor."
        )

    # Compute binding site centroid (from CA atoms of contact residues)
    res_specs = [f"{r['res_name']}{r['res_num']}" for r in residues]
    site_centroid = _centroid_from_residues(pdb_path, residues)

    return {
        "method": "reference_ligand",
        "centroid": site_centroid,
        "ligand_centroid": lig_centroid,
        "n_ligand_atoms": len(lig_coords),
        "contact_cutoff": contact_cutoff,
        "residues": residues,
        "residue_ids": res_specs,
    }


def _centroid_from_residues(
        pdb_path: str,
        residues: List[Dict[str, Any]],
) -> Tuple[float, float, float]:
    """Compute centroid from CA atoms of specified residues."""
    target_set = {(r["res_name"], r["res_num"], r["chain"]) for r in residues}
    ca_coords = []
    fallback_coords = []

    atoms = read_pdb_atoms(pdb_path, protein_only=True)
    for atom in atoms:
        key = (atom["res_name"], atom["res_num"], atom["chain"])
        if key in target_set:
            if atom["atom_name"] == "CA":
                ca_coords.append((atom["x"], atom["y"], atom["z"]))
            elif not any(f[0] == key[0] and f[1] == key[1] for f in fallback_coords):
                fallback_coords.append((atom["x"], atom["y"], atom["z"]))

    coords = ca_coords if ca_coords else fallback_coords
    if not coords:
        raise ValueError("No CA atoms found for contact residues")

    return compute_centroid(coords)


# =============================================================================
# METHOD B: BINDING SITE FROM RESIDUE LIST
# =============================================================================

def binding_site_from_residues(
        pdb_path: str,
        residue_ids: List[str],
        chain: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Method B: Define binding site from explicit residue list.

    Args:
        pdb_path: Path to receptor PDB
        residue_ids: e.g. ["GLU529", "TRP555", "CYS574"]
        chain: Chain ID filter

    Returns:
        Dict with: centroid, residues, residue_ids
    """
    # Parse residue specs
    specs = []
    for r in residue_ids:
        match = re.match(r"([A-Z]{3})(\d+)", r)
        if match:
            specs.append({"res_name": match.group(1), "res_num": int(match.group(2))})
        else:
            logger.warning(f"  Cannot parse residue: {r} (expected GLU529)")

    if not specs:
        raise ValueError(f"No valid residues parsed from: {residue_ids}")

    # Find atoms and compute centroid
    residues_found = []
    atoms = read_pdb_atoms(pdb_path, protein_only=True)

    for spec in specs:
        found = False
        for atom in atoms:
            if atom["res_name"] == spec["res_name"] and atom["res_num"] == spec["res_num"]:
                if chain and atom["chain"] != chain:
                    continue
                if not found:
                    residues_found.append({
                        "res_name": spec["res_name"],
                        "res_num": spec["res_num"],
                        "chain": atom["chain"],
                        "min_distance": 0.0,
                        "n_contacts": 0,
                    })
                    found = True

    if not residues_found:
        raise ValueError(f"No matching residues found in PDB for: {residue_ids}")

    centroid = _centroid_from_residues(pdb_path, residues_found)

    return {
        "method": "residues",
        "centroid": centroid,
        "residues": residues_found,
        "residue_ids": residue_ids,
    }


# =============================================================================
# METHOD C: BINDING SITE FROM COORDINATES
# =============================================================================

def binding_site_from_coordinates(
        center: List[float],
) -> Dict[str, Any]:
    """
    Method C: Define binding site from explicit coordinates.

    Returns:
        Dict with: centroid, method
    """
    if len(center) != 3:
        raise ValueError(f"Center must be [x, y, z], got: {center}")

    return {
        "method": "coordinates",
        "centroid": tuple(center),
        "residues": [],
        "residue_ids": [],
    }


# =============================================================================
# PDB TRIMMING
# =============================================================================

def trim_pdb_by_radius(
        input_pdb: str,
        output_pdb: str,
        center: Tuple[float, float, float],
        radius: float = 25.0,
        keep_whole_residues: bool = True,
) -> Dict[str, Any]:
    """
    Trim a PDB file to atoms within a radius of a center point.

    If keep_whole_residues=True, includes ALL atoms of any residue
    that has at least one atom within the radius. This avoids broken
    residues in the trimmed PDB.

    Args:
        input_pdb: Path to input PDB (rec_noH.pdb)
        output_pdb: Path for trimmed output PDB
        center: (x, y, z) center point
        radius: Radius in Angstroms
        keep_whole_residues: Keep all atoms of residues that touch the sphere

    Returns:
        Dict with trimming statistics
    """
    cx, cy, cz = center
    radius_sq = radius ** 2

    # First pass: find residues within radius
    residues_in_sphere = set()
    atoms_in_sphere = 0
    total_atoms = 0

    with open(input_pdb) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            total_atoms += 1
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                dsq = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
                if dsq <= radius_sq:
                    atoms_in_sphere += 1
                    if keep_whole_residues:
                        chain = line[21].strip()
                        res_num = line[22:26].strip()
                        res_name = line[17:20].strip()
                        residues_in_sphere.add((chain, res_num, res_name))
            except (ValueError, IndexError):
                continue

    # Second pass: write trimmed PDB
    output_lines = []
    atoms_written = 0

    with open(input_pdb) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if keep_whole_residues:
                    chain = line[21].strip()
                    res_num = line[22:26].strip()
                    res_name = line[17:20].strip()
                    if (chain, res_num, res_name) in residues_in_sphere:
                        output_lines.append(line)
                        atoms_written += 1
                else:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        dsq = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
                        if dsq <= radius_sq:
                            output_lines.append(line)
                            atoms_written += 1
                    except (ValueError, IndexError):
                        continue
            elif line.startswith(("TER", "END", "REMARK", "CRYST", "HEADER")):
                output_lines.append(line)

    # Add END if not present
    if output_lines and not output_lines[-1].startswith("END"):
        output_lines.append("END\n")

    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as f:
        f.writelines(output_lines)

    stats = {
        "total_atoms": total_atoms,
        "atoms_in_sphere": atoms_in_sphere,
        "residues_in_sphere": len(residues_in_sphere),
        "atoms_written": atoms_written,
        "center": list(center),
        "radius": radius,
        "keep_whole_residues": keep_whole_residues,
    }

    logger.info(f"  Trimmed PDB: {total_atoms} -> {atoms_written} atoms "
                f"({len(residues_in_sphere)} residues within {radius} A)")

    return stats


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

def run_binding_site_definition(
        receptor_noH_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        method: str = "reference_ligand",
        reference_mol2: Optional[str] = None,
        residue_ids: Optional[List[str]] = None,
        center: Optional[List[float]] = None,
        chain: Optional[str] = None,
        contact_cutoff: float = 5.0,
        trim_radius: float = 25.0,
        keep_whole_residues: bool = True,
) -> Dict[str, Any]:
    """
    Run the binding site definition pipeline.

    Identifies the binding site, finds contact residues, computes
    centroid, and trims the receptor PDB for grid generation.

    Args:
        receptor_noH_pdb: Path to receptor PDB without H (from 00b)
        output_dir: Directory for outputs
        method: "reference_ligand" | "residues" | "coordinates"
        reference_mol2: Path to reference ligand mol2 (method A)
        residue_ids: List of residue IDs (method B)
        center: [x, y, z] coordinates (method C)
        chain: Chain ID filter
        contact_cutoff: Angstroms for residue contacts (method A)
        trim_radius: Angstroms for PDB trimming
        keep_whole_residues: Keep complete residues in trimmed PDB

    Returns:
        Dict with: success, rec_noH_site_pdb, centroid, residues, report
    """
    receptor_noH_pdb = str(receptor_noH_pdb)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not Path(receptor_noH_pdb).exists():
        return {"success": False, "error": f"Receptor PDB not found: {receptor_noH_pdb}"}

    logger.info("=" * 60)
    logger.info("  Binding Site Definition")
    logger.info("=" * 60)
    logger.info(f"  Receptor:    {Path(receptor_noH_pdb).name}")
    logger.info(f"  Method:      {method}")
    logger.info(f"  Trim radius: {trim_radius} A")

    report = {
        "start_time": datetime.now().isoformat(),
        "receptor_pdb": receptor_noH_pdb,
        "method": method,
    }

    # --- Step 1: Define binding site ---
    logger.info("\nStep 1: Identifying binding site")

    try:
        if method == "reference_ligand":
            if not reference_mol2 or not Path(reference_mol2).exists():
                return {"success": False,
                        "error": f"Reference ligand not found: {reference_mol2}"}
            site = binding_site_from_ligand(
                receptor_noH_pdb, reference_mol2, contact_cutoff,
            )

        elif method == "residues":
            if not residue_ids:
                return {"success": False, "error": "No residues specified"}
            site = binding_site_from_residues(receptor_noH_pdb, residue_ids, chain)

        elif method == "coordinates":
            if not center or len(center) != 3:
                return {"success": False,
                        "error": f"Invalid center: {center}"}
            site = binding_site_from_coordinates(center)

        else:
            return {"success": False,
                    "error": f"Unknown method: {method}"}

    except ValueError as e:
        return {"success": False, "error": str(e)}

    centroid = site["centroid"]
    report["binding_site"] = {
        "centroid": list(centroid),
        "method": site["method"],
        "n_contact_residues": len(site.get("residues", [])),
        "residue_ids": site.get("residue_ids", []),
    }

    if site.get("ligand_centroid"):
        report["binding_site"]["ligand_centroid"] = list(site["ligand_centroid"])
        report["binding_site"]["n_ligand_atoms"] = site.get("n_ligand_atoms", 0)
        report["binding_site"]["contact_cutoff"] = site.get("contact_cutoff", 0)

    logger.info(f"  Centroid: ({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})")
    logger.info(f"  Contact residues: {len(site.get('residues', []))}")

    # --- Step 2: Trim receptor PDB ---
    logger.info(f"\nStep 2: Trimming receptor PDB (radius={trim_radius} A)")

    rec_site_pdb = output_dir / "rec_noH_site.pdb"
    trim_stats = trim_pdb_by_radius(
        input_pdb=receptor_noH_pdb,
        output_pdb=str(rec_site_pdb),
        center=centroid,
        radius=trim_radius,
        keep_whole_residues=keep_whole_residues,
    )
    report["trim_stats"] = trim_stats

    # --- Step 3: Save report ---
    report["end_time"] = datetime.now().isoformat()
    report["success"] = True

    report_path = output_dir / "binding_site_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    # --- Summary TXT ---
    summary_path = output_dir / "binding_site_summary.txt"
    w = 70
    lines = [
        "=" * w,
        "00e BINDING SITE DEFINITION - SUMMARY",
        "=" * w,
        "",
        f"Date:              {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Receptor:          {Path(receptor_noH_pdb).name}",
        f"Method:            {method}",
        "",
        f"Centroid:          ({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})",
        f"Contact residues:  {len(site.get('residues', []))}",
        f"Trim radius:       {trim_radius} A",
        "",
        f"Atoms (original):  {trim_stats['total_atoms']}",
        f"Atoms (trimmed):   {trim_stats['atoms_written']}",
        f"Residues (kept):   {trim_stats['residues_in_sphere']}",
        "",
    ]

    if site.get("ligand_centroid"):
        lc = site["ligand_centroid"]
        lines.extend([
            f"Reference ligand:  {Path(reference_mol2).name}",
            f"Ligand centroid:   ({lc[0]:.3f}, {lc[1]:.3f}, {lc[2]:.3f})",
            f"Ligand atoms:      {site.get('n_ligand_atoms', 0)}",
            f"Contact cutoff:    {contact_cutoff} A",
            "",
        ])

    # Residue table
    residues = site.get("residues", [])
    if residues:
        lines.extend([
            "-" * w,
            f"{'Residue':<12} {'Num':>5} {'Chain':>6} {'Min Dist':>10} {'Contacts':>10}",
            "-" * w,
        ])
        for r in residues:
            lines.append(
                f"{r['res_name']:<12} {r['res_num']:>5} {r['chain']:>6} "
                f"{r['min_distance']:>10.2f} {r['n_contacts']:>10}"
            )
        lines.append("")

    lines.extend([
        "=" * w,
        "",
        f"rec_noH_site.pdb:  {rec_site_pdb}",
    ])

    summary_path.write_text("\n".join(lines))
    logger.info(f"  Summary: {summary_path}")

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  Binding site defined: {len(residues)} contact residues")
    logger.info(f"  Trimmed PDB: {trim_stats['atoms_written']} atoms")
    logger.info(f"  rec_noH_site.pdb: {rec_site_pdb}")
    logger.info("=" * 60)

    return {
        "success": True,
        "rec_noH_site_pdb": str(rec_site_pdb),
        "centroid": list(centroid),
        "residues": residues,
        "residue_ids": site.get("residue_ids", []),
        "trim_stats": trim_stats,
        "report_path": str(report_path),
        "summary_path": str(summary_path),
    }