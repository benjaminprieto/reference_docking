"""
Pharmit Pharmacophore Generator (06a) v5.0
============================================
Generate ranked sub-pocket pharmacophore queries from PLIP interactions
AND DOCK6 footprint per-residue energy decomposition.

v5.0 change: PLIP + Footprint IS the pharmacophore.
Previously, only PLIP interactions generated features. This missed
vdW-dominated contacts (e.g., TRP392/TRP495 CH-π stacking with xylose)
that PLIP doesn't detect geometrically but the footprint confirms
energetically.

Pipeline:
    1. Load PLIP interactions JSON → pharmacophore points (if available)
    2. Load DOCK6 footprint (residue_consensus.csv from 04b)
    3. For footprint residues NOT covered by PLIP: generate new features
       from ligand-residue contact geometry
    4. Merge, deduplicate, classify sub-pockets
    5. Rank by footprint energy, filter repulsive
    6. Generate Pharmit JSONs by strategy

Feature generation from footprint:
    For each residue with energy < cutoff and no PLIP feature:
        - Find receptor residue atoms
        - Find closest ligand atoms (within contact_distance)
        - Centroid of closest ligand atoms = feature coordinate
        - Feature type inferred from residue chemistry + energy profile:
            Aromatic residues (TRP, PHE, TYR, HIS) + vdW dominant → Aromatic
            Aliphatic (VAL, LEU, ILE, PRO) + vdW dominant → Hydrophobic
            Positive (ARG, LYS) + ES dominant → NegativeIon (ligand)
            Negative (ASP, GLU) + ES dominant → PositiveIon (ligand)
            Polar (SER, THR) → HydrogenAcceptor

Input:
    Required: DOCK6 residue_consensus.csv (from reference_docking 04b)
    Required: Reference ligand mol2 (crystal coordinates)
    Required: Receptor mol2 (for residue coordinates)
    Optional: PLIP interactions JSON (from 03a — adds geometric detail)

Output:
    pharmacophore_{strategy}.json  — Pharmit queries (5 strategies)
    pharmacophore_ranking.csv      — Decision table
    pharmacophore_ranking.txt      — Human-readable summary

Location: 01_src/reference_docking/m06_pharmit/pharmit_pharmacophore.py
Project: reference_docking
Module: 06a (core)
Version: 5.0 — Footprint + PLIP (2026-03-26)
"""

import csv
import json
import logging
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

# Sub-pocket classification by receptor residue (PDB numbering)
XYLOSE_RESIDUES = {"TRP495", "ASP494", "TRP392", "SER575", "TYR565", "CYS574", "GLU529"}
URACIL_RESIDUES = {"ASP361", "THR390", "ARG363"}
RIBOSE_RESIDUES = {"HIS335", "VAL333", "VAL334", "THR390"}
PHOSPHATE_RESIDUES = {"ARG598", "LYS599"}
CATALYTIC_RESIDUES = {"GLU529"}

# Residue 3-letter code → chemistry class
AROMATIC_RESIDUES = {"TRP", "PHE", "TYR", "HIS"}
ALIPHATIC_RESIDUES = {"VAL", "LEU", "ILE", "ALA", "PRO", "MET", "CYS", "GLY"}
POSITIVE_RESIDUES = {"ARG", "LYS"}
NEGATIVE_RESIDUES = {"ASP", "GLU"}
POLAR_RESIDUES = {"SER", "THR", "ASN", "GLN"}

# PLIP interaction type → Pharmit feature type
PLIP_TO_PHARMIT = {
    ("hbond", False): "HydrogenAcceptor",
    ("hbond", True): "HydrogenDonor",
    ("salt_bridge", "negative"): "NegativeIon",
    ("salt_bridge", "positive"): "PositiveIon",
    ("pi_stack", None): "Aromatic",
    ("hydrophobic", None): "Hydrophobic",
    ("pi_cation", None): "Aromatic",
    ("water_bridge", None): "HydrogenAcceptor",
    ("halogen_bond", None): "HydrogenAcceptor",
}

# Strategies for Pharmit query generation
STRATEGIES = {
    "xylose": {"include": ["xylose", "catalytic"], "phosphate_as": None},
    "uracil": {"include": ["uracil", "ribose"], "phosphate_as": None},
    "combined": {"include": ["xylose", "uracil", "ribose", "catalytic"], "phosphate_as": None},
    "analogues": {"include": ["xylose", "uracil", "ribose", "phosphate", "catalytic"], "phosphate_as": "NegativeIon"},
    "druglike": {"include": ["xylose", "uracil", "ribose", "phosphate", "catalytic"], "phosphate_as": "HydrogenAcceptor"},
}


# ═══════════════════════════════════════════════════════════════════════
# DATA STRUCTURES
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class PharmacophorePoint:
    """One pharmacophore feature point."""
    index: int
    pharmit_type: str  # HydrogenAcceptor, HydrogenDonor, Aromatic, etc.
    x: float
    y: float
    z: float

    # Source
    source: str = "footprint"  # "plip", "footprint", or "plip+footprint"

    # Residue evidence
    residues: List[str] = field(default_factory=list)
    interaction_types: List[str] = field(default_factory=list)
    n_interactions: int = 0
    plip_details: List[Dict] = field(default_factory=list)

    # Direction vector
    svector: Dict[str, float] = field(default_factory=lambda: {"x": 0.0, "y": 0.0, "z": 1.0})

    # Footprint energy (from 04b)
    dock6_energy: float = 0.0
    dock6_vdw: float = 0.0
    dock6_es: float = 0.0
    dock6_freq: float = 0.0

    # Confidence: vdW fraction of total energy
    vdw_fraction: float = 0.0

    # Sub-pocket
    sub_pocket: str = "unknown"

    # Ranking
    rank: int = 0
    priority: str = "DISABLED"


# ═══════════════════════════════════════════════════════════════════════
# MOL2 PARSING
# ═══════════════════════════════════════════════════════════════════════

def parse_mol2_atoms(mol2_path: str) -> List[Dict[str, Any]]:
    """
    Parse atom coordinates from a mol2 file.

    Returns:
        List of dicts with keys: atom_id, atom_name, x, y, z,
        atom_type, residue_id, residue_name
    """
    atoms = []
    in_atom = False
    with open(mol2_path, "r") as f:
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    atoms.append({
                        "atom_id": int(parts[0]),
                        "atom_name": parts[1],
                        "x": float(parts[2]),
                        "y": float(parts[3]),
                        "z": float(parts[4]),
                        "atom_type": parts[5],
                        "residue_id": int(parts[6]),
                        "residue_name": parts[7],
                    })
    return atoms


def build_pdb_to_mol2_residue_map(
        receptor_atoms: List[Dict],
        receptor_pdb_path: Optional[str] = None,
) -> Dict[str, str]:
    """
    Build a mapping from PDB residue names (TRP392) to mol2 residue names (TRP141).

    ChimeraX mol2 uses sequential numbering. This function maps by matching
    the 3-letter residue code and atom coordinates between PDB and mol2.

    If no PDB is available, returns identity mapping (PDB name = mol2 name).
    """
    if not receptor_pdb_path or not Path(receptor_pdb_path).exists():
        return {}

    # Read PDB CA/heavy atoms with their residue info
    pdb_residues = {}  # "TRP392.A" -> (x, y, z) of first heavy atom
    with open(receptor_pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name != "CA":
                    continue
                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                chain = line[21].strip() or "A"
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                key = f"{res_name}{res_num}"
                if key not in pdb_residues:
                    pdb_residues[key] = (x, y, z)

    # Build mol2 residue → CA coords
    # ChimeraX mol2: residue_name="TRP" (no number), residue_id=141 (sequential int)
    # Key must combine both: "TRP141"
    mol2_residues = {}
    for atom in receptor_atoms:
        mol2_key = f"{atom['residue_name']}{atom['residue_id']}"
        if atom["atom_name"] == "CA" and mol2_key not in mol2_residues:
            mol2_residues[mol2_key] = (atom["x"], atom["y"], atom["z"])

    # Match by coordinate proximity (CA positions should be identical)
    mapping = {}
    for pdb_key, (px, py, pz) in pdb_residues.items():
        best_mol2 = None
        best_dist = 0.5  # tolerance
        for mol2_key, (mx, my, mz) in mol2_residues.items():
            # Must have same 3-letter code
            pdb_res3 = "".join(c for c in pdb_key if c.isalpha())
            mol2_res3 = "".join(c for c in mol2_key if c.isalpha())
            if pdb_res3 != mol2_res3:
                continue
            dist = ((px - mx)**2 + (py - my)**2 + (pz - mz)**2) ** 0.5
            if dist < best_dist:
                best_dist = dist
                best_mol2 = mol2_key
        if best_mol2:
            mapping[pdb_key] = best_mol2

    return mapping


def get_residue_atoms(receptor_atoms: List[Dict], residue_name: str,
                      chain: str = "A",
                      pdb_to_mol2_map: Optional[Dict[str, str]] = None,
                      ) -> List[Dict]:
    """
    Get all atoms for a given residue from the receptor.

    residue_name: format "RES###" e.g. "TRP392" or "ARG598"
    pdb_to_mol2_map: mapping from PDB residue names to mol2 residue names
    """
    # Try mapping PDB name to mol2 name
    mol2_name = residue_name
    if pdb_to_mol2_map:
        mol2_name = pdb_to_mol2_map.get(residue_name, residue_name)

    # Parse residue 3-letter code and number from mol2_name (e.g. "TRP141")
    res3 = ""
    resnum = ""
    for i, ch in enumerate(mol2_name):
        if ch.isdigit():
            res3 = mol2_name[:i]
            resnum = mol2_name[i:]
            break

    if not res3 or not resnum:
        return []

    # Match by residue_name + residue_id in mol2
    # ChimeraX mol2: residue_name="TRP", residue_id=141 (separate fields)
    matches = []
    try:
        resnum_int = int(resnum)
    except ValueError:
        return []

    for atom in receptor_atoms:
        if atom["residue_name"] == res3 and atom["residue_id"] == resnum_int:
            matches.append(atom)

    return matches


def get_ligand_atoms_near_residue(
        ligand_atoms: List[Dict],
        residue_atoms: List[Dict],
        contact_distance: float = 4.5,
) -> List[Dict]:
    """
    Find ligand atoms within contact_distance of any residue atom.

    Returns ligand atoms sorted by minimum distance to residue.
    """
    if not residue_atoms or not ligand_atoms:
        return []

    res_coords = np.array([[a["x"], a["y"], a["z"]] for a in residue_atoms])
    contacts = []

    for lig_atom in ligand_atoms:
        lig_pos = np.array([lig_atom["x"], lig_atom["y"], lig_atom["z"]])
        min_dist = np.min(np.linalg.norm(res_coords - lig_pos, axis=1))
        if min_dist <= contact_distance:
            contacts.append((lig_atom, min_dist))

    contacts.sort(key=lambda x: x[1])
    return [c[0] for c in contacts]


# ═══════════════════════════════════════════════════════════════════════
# FEATURE TYPE INFERENCE FROM FOOTPRINT
# ═══════════════════════════════════════════════════════════════════════

def infer_feature_type(residue_name: str, vdw: float, es: float) -> str:
    """
    Infer Pharmit feature type from receptor residue chemistry
    and energy decomposition.

    The feature type describes the LIGAND atom, not the receptor.
    Example: ARG (positive receptor) interacting with ligand →
    ligand must have negative character → NegativeIon or HydrogenAcceptor.

    Args:
        residue_name: e.g. "TRP392" (just need the 3-letter code)
        vdw: van der Waals energy contribution
        es: electrostatic energy contribution

    Returns:
        Pharmit feature type string
    """
    # Extract 3-letter residue code
    res3 = ""
    for i, ch in enumerate(residue_name):
        if ch.isdigit():
            res3 = residue_name[:i].upper()
            break
    if not res3:
        res3 = residue_name[:3].upper()

    total = abs(vdw) + abs(es) if (abs(vdw) + abs(es)) > 0 else 1.0
    vdw_frac = abs(vdw) / total

    # Aromatic residues with vdW dominant → stacking
    if res3 in AROMATIC_RESIDUES and vdw_frac > 0.5:
        return "Aromatic"

    # Aromatic residues with ES dominant → H-bond to ring
    if res3 in AROMATIC_RESIDUES and vdw_frac <= 0.5:
        return "HydrogenAcceptor"

    # Aliphatic residues → hydrophobic contact
    if res3 in ALIPHATIC_RESIDUES:
        return "Hydrophobic"

    # Positive receptor residues (ARG, LYS) → ligand is negative/acceptor
    if res3 in POSITIVE_RESIDUES:
        if abs(es) > 3.0:
            return "NegativeIon"
        return "HydrogenAcceptor"

    # Negative receptor residues (ASP, GLU) → ligand is positive/donor
    if res3 in NEGATIVE_RESIDUES:
        if abs(es) > 3.0:
            return "PositiveIon"
        return "HydrogenDonor"

    # Polar residues → H-bond
    if res3 in POLAR_RESIDUES:
        return "HydrogenAcceptor"

    # Default
    return "Hydrophobic"


def classify_sub_pocket(residue_name: str) -> str:
    """Classify a residue into a sub-pocket."""
    # Strip chain suffix if present (e.g. "TRP392.A" → "TRP392")
    res_base = residue_name.split(".")[0]

    if res_base in XYLOSE_RESIDUES:
        return "xylose"
    if res_base in URACIL_RESIDUES:
        return "uracil"
    if res_base in RIBOSE_RESIDUES:
        return "ribose"
    if res_base in PHOSPHATE_RESIDUES:
        return "phosphate"
    if res_base in CATALYTIC_RESIDUES:
        return "catalytic"
    return "other"


def compute_direction_vector(lig_coords: np.ndarray,
                             rec_coords: np.ndarray,
                             pharmit_type: str) -> Dict[str, float]:
    """Compute direction vector for Pharmit."""
    default = {"x": 0.0, "y": 0.0, "z": 1.0}
    if lig_coords is None or rec_coords is None:
        return default
    if len(lig_coords) == 0 or len(rec_coords) == 0:
        return default

    lig_center = np.mean(lig_coords, axis=0) if lig_coords.ndim > 1 else lig_coords
    rec_center = np.mean(rec_coords, axis=0) if rec_coords.ndim > 1 else rec_coords

    if pharmit_type in ("HydrogenAcceptor", "NegativeIon"):
        vec = rec_center - lig_center
    elif pharmit_type == "HydrogenDonor":
        vec = lig_center - rec_center
    elif pharmit_type == "Aromatic":
        return default
    else:
        vec = rec_center - lig_center

    norm = np.linalg.norm(vec)
    if norm < 1e-6:
        return default
    vec /= norm
    return {"x": round(float(vec[0]), 6),
            "y": round(float(vec[1]), 6),
            "z": round(float(vec[2]), 6)}


# ═══════════════════════════════════════════════════════════════════════
# STEP 1: LOAD PLIP (optional)
# ═══════════════════════════════════════════════════════════════════════

def load_plip_points(json_path: str,
                     group_tolerance: float = 1.0,
                     include_hydrophobic: bool = False,
                     ) -> List[PharmacophorePoint]:
    """
    Load PLIP interactions and convert to pharmacophore points.
    Returns empty list if PLIP JSON not available.
    """
    if not json_path or not Path(json_path).exists():
        logger.info("  PLIP JSON not available — skipping")
        return []

    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    interactions = data.get("interactions", [])
    if not interactions:
        logger.info("  PLIP JSON has no interactions")
        return []

    logger.info(f"  Loaded {len(interactions)} PLIP interactions")

    # Group by ligand coordinates (within tolerance)
    groups = []
    for inter in interactions:
        coords = inter.get("ligand_coords", inter.get("coords", []))
        if not coords or len(coords) < 3:
            continue
        itype = inter.get("interaction_type", "")
        if itype == "hydrophobic" and not include_hydrophobic:
            continue

        pos = np.array(coords[:3], dtype=np.float64)

        merged = False
        for g_centroid, g_interactions in groups:
            if np.linalg.norm(pos - g_centroid) <= group_tolerance:
                g_interactions.append(inter)
                all_coords = [np.array(
                    i.get("ligand_coords", i.get("coords", [0, 0, 0]))[:3]
                ) for i in g_interactions]
                g_centroid[:] = np.mean(all_coords, axis=0)
                merged = True
                break
        if not merged:
            groups.append((pos.copy(), [inter]))

    # Convert to PharmacophorePoints
    points = []
    for idx, (centroid, inters) in enumerate(groups):
        # Determine feature type from dominant interaction
        itypes = [_classify_plip_interaction(i) for i in inters]
        type_counts = Counter(itypes)
        pharmit_type = type_counts.most_common(1)[0][0]

        residues = list(set(
            i.get("receptor_residue", i.get("residue", "UNK"))
            for i in inters
        ))

        plip_details = []
        for i in inters:
            plip_details.append({
                "type": i.get("interaction_type", ""),
                "residue": i.get("receptor_residue", i.get("residue", "")),
                "receptor_atom": i.get("receptor_atom", ""),
                "distance": i.get("distance", 0.0),
            })

        pt = PharmacophorePoint(
            index=idx,
            pharmit_type=pharmit_type,
            x=round(float(centroid[0]), 4),
            y=round(float(centroid[1]), 4),
            z=round(float(centroid[2]), 4),
            source="plip",
            residues=residues,
            interaction_types=list(set(i.get("interaction_type", "") for i in inters)),
            n_interactions=len(inters),
            plip_details=plip_details,
        )
        points.append(pt)

    logger.info(f"  PLIP → {len(points)} pharmacophore points")
    return points


def _classify_plip_interaction(inter: Dict) -> str:
    """Determine Pharmit feature type from a PLIP interaction."""
    itype = inter.get("interaction_type", "")
    if itype == "hbond":
        is_donor = inter.get("ligand_is_donor", False)
        return PLIP_TO_PHARMIT.get(("hbond", is_donor), "HydrogenAcceptor")
    if itype == "salt_bridge":
        charge = inter.get("ligand_charge", "negative")
        return PLIP_TO_PHARMIT.get(("salt_bridge", charge), "NegativeIon")
    if itype == "pi_stack":
        return "Aromatic"
    if itype == "hydrophobic":
        return "Hydrophobic"
    if itype == "pi_cation":
        return "Aromatic"
    return "HydrogenAcceptor"


# ═══════════════════════════════════════════════════════════════════════
# STEP 2: LOAD FOOTPRINT + GENERATE FEATURES
# ═══════════════════════════════════════════════════════════════════════

def load_footprint_as_pharmacophore(
        residue_csv_path: str,
        ligand_mol2_path: str,
        receptor_mol2_path: str,
        energy_cutoff: float = -0.5,
        contact_distance: float = 4.5,
        existing_residues: Optional[set] = None,
        receptor_pdb_path: Optional[str] = None,
) -> List[PharmacophorePoint]:
    """
    Generate pharmacophore points from DOCK6 footprint data.

    For each residue with energy < cutoff that is NOT already covered
    by PLIP (existing_residues), create a new pharmacophore point at
    the centroid of the closest ligand atoms.

    Args:
        residue_csv_path:  Path to residue_consensus.csv (from 04b)
        ligand_mol2_path:  Path to ligand mol2 (crystal coordinates)
        receptor_mol2_path: Path to receptor mol2
        energy_cutoff:     Only include residues below this energy
        contact_distance:  Max distance for ligand-residue contact (Å)
        existing_residues: Residues already covered by PLIP features

    Returns:
        List of PharmacophorePoints from footprint
    """
    if existing_residues is None:
        existing_residues = set()

    # Load footprint data
    df = pd.read_csv(residue_csv_path)
    logger.info(f"  Footprint: {len(df)} residues")

    # Parse mol2 files
    lig_atoms = parse_mol2_atoms(ligand_mol2_path)
    rec_atoms = parse_mol2_atoms(receptor_mol2_path)
    logger.info(f"  Ligand: {len(lig_atoms)} atoms, Receptor: {len(rec_atoms)} atoms")

    if not lig_atoms or not rec_atoms:
        logger.error("  Cannot parse mol2 files")
        return []

    # Build PDB→mol2 residue name mapping (ChimeraX uses sequential numbering)
    pdb_to_mol2 = build_pdb_to_mol2_residue_map(rec_atoms, receptor_pdb_path)
    if pdb_to_mol2:
        logger.info(f"  PDB→mol2 residue mapping: {len(pdb_to_mol2)} residues mapped")

    # Use explicit column names matching residue_consensus.csv format
    # Same approach as assign_energy_from_footprint (hardcoded, no auto-detect)
    residue_col = "residue_id" if "residue_id" in df.columns else df.columns[0]
    energy_col = "mean_total" if "mean_total" in df.columns else None
    vdw_col = "mean_vdw" if "mean_vdw" in df.columns else ("ref_vdw" if "ref_vdw" in df.columns else None)
    es_col = "mean_es" if "mean_es" in df.columns else ("ref_es" if "ref_es" in df.columns else None)
    freq_col = "frac_contributing" if "frac_contributing" in df.columns else None

    if energy_col is None:
        # Try computing from vdW + ES
        if vdw_col and es_col:
            df["_total"] = df[vdw_col].fillna(0) + df[es_col].fillna(0)
            energy_col = "_total"
        else:
            logger.error(f"  Cannot find energy column. Columns: {list(df.columns)}")
            return []

    logger.info(f"  Columns: residue={residue_col}, energy={energy_col}, "
                f"vdw={vdw_col}, es={es_col}, freq={freq_col}")

    # Filter significant residues
    df_sig = df[df[energy_col].astype(float) < energy_cutoff].copy()
    logger.info(f"  Significant residues (< {energy_cutoff} kcal/mol): {len(df_sig)}")

    # Generate points
    points = []
    n_plip_covered = 0
    n_no_contact = 0

    for _, row in df_sig.iterrows():
        res_id = str(row[residue_col])

        # Strip chain suffix for comparison
        res_base = res_id.split(".")[0]

        # Skip if already covered by PLIP
        if res_base in existing_residues or res_id in existing_residues:
            n_plip_covered += 1
            continue

        vdw = float(row[vdw_col]) if vdw_col and pd.notna(row.get(vdw_col)) else 0.0
        es = float(row[es_col]) if es_col and pd.notna(row.get(es_col)) else 0.0
        total = float(row[energy_col])
        freq = float(row[freq_col]) if freq_col and pd.notna(row.get(freq_col)) else 1.0

        # Get receptor residue atoms (use PDB→mol2 mapping for ChimeraX numbering)
        res_rec_atoms = get_residue_atoms(rec_atoms, res_base, pdb_to_mol2_map=pdb_to_mol2)
        if not res_rec_atoms:
            # Try with chain suffix
            res_rec_atoms = get_residue_atoms(rec_atoms, res_id.replace(".", ""), pdb_to_mol2_map=pdb_to_mol2)
        if not res_rec_atoms:
            logger.debug(f"    {res_id}: no receptor atoms found")
            n_no_contact += 1
            continue

        # Find closest ligand atoms
        contact_lig_atoms = get_ligand_atoms_near_residue(
            lig_atoms, res_rec_atoms, contact_distance,
        )
        if not contact_lig_atoms:
            logger.debug(f"    {res_id}: no ligand contacts within {contact_distance}Å")
            n_no_contact += 1
            continue

        # Feature coordinate = centroid of closest ligand atoms (top 3)
        top_lig = contact_lig_atoms[:3]
        lig_coords = np.array([[a["x"], a["y"], a["z"]] for a in top_lig])
        centroid = np.mean(lig_coords, axis=0)

        # Feature type from residue chemistry + energy profile
        pharmit_type = infer_feature_type(res_base, vdw, es)

        # Direction vector
        rec_coords = np.array([[a["x"], a["y"], a["z"]] for a in res_rec_atoms])
        svector = compute_direction_vector(lig_coords, rec_coords, pharmit_type)

        # vdW fraction
        total_abs = abs(vdw) + abs(es)
        vdw_frac = abs(vdw) / total_abs if total_abs > 0 else 0.5

        # Sub-pocket
        sub_pocket = classify_sub_pocket(res_id)

        pt = PharmacophorePoint(
            index=len(points),
            pharmit_type=pharmit_type,
            x=round(float(centroid[0]), 4),
            y=round(float(centroid[1]), 4),
            z=round(float(centroid[2]), 4),
            source="footprint",
            residues=[res_base],
            interaction_types=[f"vdW({vdw:.1f})+ES({es:.1f})"],
            n_interactions=1,
            dock6_energy=round(total, 4),
            dock6_vdw=round(vdw, 4),
            dock6_es=round(es, 4),
            dock6_freq=round(freq, 4),
            vdw_fraction=round(vdw_frac, 4),
            sub_pocket=sub_pocket,
            svector=svector,
        )
        points.append(pt)

    logger.info(f"  Footprint → {len(points)} new features "
                f"({n_plip_covered} already in PLIP, {n_no_contact} no contacts)")
    return points


# ═══════════════════════════════════════════════════════════════════════
# STEP 3: MERGE PLIP + FOOTPRINT
# ═══════════════════════════════════════════════════════════════════════

def merge_points(plip_points: List[PharmacophorePoint],
                 footprint_points: List[PharmacophorePoint],
                 merge_distance: float = 2.0,
                 ) -> List[PharmacophorePoint]:
    """
    Merge PLIP and footprint points. If a footprint point is within
    merge_distance of a PLIP point, enrich the PLIP point with energy.
    Otherwise, add the footprint point as new.
    """
    merged = list(plip_points)

    for fp_pt in footprint_points:
        fp_pos = np.array([fp_pt.x, fp_pt.y, fp_pt.z])

        # Check overlap with existing points
        found_overlap = False
        for existing in merged:
            ex_pos = np.array([existing.x, existing.y, existing.z])
            if np.linalg.norm(fp_pos - ex_pos) <= merge_distance:
                # Enrich existing point with footprint energy
                existing.dock6_energy = fp_pt.dock6_energy
                existing.dock6_vdw = fp_pt.dock6_vdw
                existing.dock6_es = fp_pt.dock6_es
                existing.dock6_freq = fp_pt.dock6_freq
                existing.vdw_fraction = fp_pt.vdw_fraction
                existing.source = "plip+footprint"
                # Add residue if not already present
                for r in fp_pt.residues:
                    if r not in existing.residues:
                        existing.residues.append(r)
                found_overlap = True
                break

        if not found_overlap:
            fp_pt.index = len(merged)
            merged.append(fp_pt)

    # Assign sub-pockets for all points
    for pt in merged:
        if pt.sub_pocket == "unknown":
            for res in pt.residues:
                sp = classify_sub_pocket(res)
                if sp != "other":
                    pt.sub_pocket = sp
                    break

    # Re-index
    for i, pt in enumerate(merged):
        pt.index = i

    n_plip = sum(1 for p in merged if p.source == "plip")
    n_fp = sum(1 for p in merged if p.source == "footprint")
    n_both = sum(1 for p in merged if p.source == "plip+footprint")
    logger.info(f"  Merged: {len(merged)} total "
                f"(PLIP-only: {n_plip}, footprint-only: {n_fp}, both: {n_both})")

    return merged


# ═══════════════════════════════════════════════════════════════════════
# STEP 4: ASSIGN ENERGY TO PLIP-ONLY POINTS
# ═══════════════════════════════════════════════════════════════════════

def assign_energy_from_footprint(
        points: List[PharmacophorePoint],
        residue_csv_path: str,
) -> None:
    """
    For PLIP-only points that didn't get merged with footprint,
    look up their residues in the footprint CSV and assign energy.
    """
    if not Path(residue_csv_path).exists():
        return

    df = pd.read_csv(residue_csv_path)

    # Build residue → energy lookup
    # Use explicit column names matching residue_consensus.csv format
    energy_map = {}
    rid_col = "residue_id" if "residue_id" in df.columns else df.columns[0]
    vdw_col = "mean_vdw" if "mean_vdw" in df.columns else "ref_vdw"
    es_col_name = "mean_es" if "mean_es" in df.columns else "ref_es"
    total_col = "mean_total" if "mean_total" in df.columns else None
    freq_col = "frac_contributing" if "frac_contributing" in df.columns else None

    for _, row in df.iterrows():
        rid = str(row[rid_col])
        vdw = float(row[vdw_col]) if vdw_col in df.columns else 0.0
        es = float(row[es_col_name]) if es_col_name in df.columns else 0.0
        total = float(row[total_col]) if total_col and total_col in df.columns else vdw + es
        freq = float(row[freq_col]) if freq_col and freq_col in df.columns else 1.0
        energy_map[rid] = {"vdw": vdw, "es": es, "total": total, "freq": freq}

    for pt in points:
        if pt.dock6_energy != 0:
            continue  # Already has energy
        total_vdw, total_es, freqs = 0.0, 0.0, []
        for res in pt.residues:
            for key in [res, f"{res}.A", res.split(".")[0]]:
                if key in energy_map:
                    e = energy_map[key]
                    total_vdw += e["vdw"]
                    total_es += e["es"]
                    freqs.append(e["freq"])
                    break
        if total_vdw != 0 or total_es != 0:
            pt.dock6_vdw = round(total_vdw, 4)
            pt.dock6_es = round(total_es, 4)
            pt.dock6_energy = round(total_vdw + total_es, 4)
            pt.dock6_freq = round(float(np.mean(freqs)), 4) if freqs else 0.0
            total_abs = abs(total_vdw) + abs(total_es)
            pt.vdw_fraction = round(abs(total_vdw) / total_abs, 4) if total_abs > 0 else 0.5


# ═══════════════════════════════════════════════════════════════════════
# STEP 5: RANK + FILTER
# ═══════════════════════════════════════════════════════════════════════

def rank_and_filter(points: List[PharmacophorePoint],
                    energy_cutoff: float = 0.0,
                    n_required: int = 3,
                    n_optional: int = 2) -> List[PharmacophorePoint]:
    """
    Filter repulsive features and rank by energy.
    """
    # Filter repulsive
    filtered = [p for p in points if p.dock6_energy <= energy_cutoff or p.dock6_energy == 0]
    removed = len(points) - len(filtered)
    if removed:
        logger.info(f"  Filtered {removed} repulsive features (energy > {energy_cutoff})")

    # Sort: most negative energy first
    has_energy = any(p.dock6_energy != 0 for p in filtered)
    if has_energy:
        filtered.sort(key=lambda p: p.dock6_energy if p.dock6_energy != 0 else 0)
    else:
        filtered.sort(key=lambda p: p.n_interactions, reverse=True)

    # Assign rank and priority
    for i, pt in enumerate(filtered):
        pt.rank = i + 1
        if i < n_required:
            pt.priority = "REQUIRED"
        elif i < n_required + n_optional:
            pt.priority = "OPTIONAL"
        else:
            pt.priority = "DISABLED"

    return filtered


# ═══════════════════════════════════════════════════════════════════════
# STEP 6: GENERATE PHARMIT JSON
# ═══════════════════════════════════════════════════════════════════════

def generate_pharmit_json(
        points: List[PharmacophorePoint],
        strategy_name: str,
        n_required: int = 3,
        n_optional: int = 2,
        radius: float = 1.0,
) -> Dict[str, Any]:
    """Generate a Pharmit-compatible JSON for a given strategy."""
    strategy = STRATEGIES.get(strategy_name)
    if not strategy:
        return {"points": [], "_meta": {"error": f"Unknown strategy: {strategy_name}"}}

    include_pockets = strategy["include"]
    phosphate_as = strategy.get("phosphate_as")

    selected = [p for p in points if p.sub_pocket in include_pockets]

    # Re-rank within selection
    has_energy = any(p.dock6_energy != 0 for p in selected)
    if has_energy:
        selected.sort(key=lambda p: p.dock6_energy if p.dock6_energy != 0 else 0)
    else:
        selected.sort(key=lambda p: p.n_interactions, reverse=True)

    json_points = []
    n_enabled = 0
    for i, pt in enumerate(selected):
        name = pt.pharmit_type

        # Convert phosphate type if strategy says so
        if pt.sub_pocket == "phosphate" and phosphate_as is not None:
            name = phosphate_as

        enabled = i < (n_required + n_optional)
        if enabled:
            n_enabled += 1

        has_vec = name in ("HydrogenDonor", "HydrogenAcceptor", "Aromatic")

        json_points.append({
            "name": name,
            "hasvec": has_vec,
            "x": pt.x, "y": pt.y, "z": pt.z,
            "radius": radius,
            "enabled": enabled,
            "vector_on": 0,
            "svector": pt.svector,
            "minsize": "", "maxsize": "",
            "selected": False,
        })

    return {
        "points": json_points,
        "_meta": {
            "strategy": strategy_name,
            "n_features": len(json_points),
            "n_enabled": n_enabled,
            "sub_pockets": include_pockets,
            "phosphate_as": phosphate_as,
        },
    }


# ═══════════════════════════════════════════════════════════════════════
# OUTPUT WRITERS
# ═══════════════════════════════════════════════════════════════════════

def write_ranking_csv(points: List[PharmacophorePoint], path: Path):
    """Write decision table CSV."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Rank", "Priority", "Source", "Feature", "Sub_Pocket",
            "x", "y", "z",
            "Residues", "N_Interactions", "Interaction_Types",
            "DOCK6_Energy", "DOCK6_vdW", "DOCK6_ES", "DOCK6_Freq",
            "vdW_Fraction", "Confidence", "Details",
        ])
        for pt in points:
            # Confidence based on vdW fraction and energy
            if pt.dock6_energy < -5.0 and pt.vdw_fraction > 0.6:
                confidence = "HIGH"
            elif pt.dock6_energy < -1.0:
                confidence = "MEDIUM"
            elif pt.dock6_energy != 0:
                confidence = "LOW"
            else:
                confidence = "UNKNOWN"

            details = "; ".join(
                f"{d['type']}→{d['residue']}:{d['receptor_atom']}@{d['distance']:.1f}A"
                for d in pt.plip_details
            ) if pt.plip_details else f"footprint: vdW={pt.dock6_vdw:.1f} ES={pt.dock6_es:.1f}"

            writer.writerow([
                pt.rank, pt.priority, pt.source, pt.pharmit_type, pt.sub_pocket,
                pt.x, pt.y, pt.z,
                ";".join(pt.residues), pt.n_interactions,
                ";".join(pt.interaction_types),
                pt.dock6_energy, pt.dock6_vdw, pt.dock6_es, pt.dock6_freq,
                pt.vdw_fraction, confidence, details,
            ])
    logger.info(f"  Saved: {path}")


def write_ranking_txt(points: List[PharmacophorePoint],
                      ligand_name: str, strategies: List[str],
                      path: Path):
    """Write human-readable ranking summary."""
    w = 100
    has_energy = any(p.dock6_energy != 0 for p in points)

    lines = [
        "=" * w,
        "06a PHARMIT PHARMACOPHORE — FEATURE RANKING (v5.0)",
        "=" * w,
        "",
        f"Date:         {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Ligand:       {ligand_name}",
        f"Sources:      PLIP interactions + DOCK6 footprint",
        f"Ranking by:   {'DOCK6 footprint energy' if has_energy else 'PLIP interaction count'}",
        f"Features:     {len(points)}",
        f"Strategies:   {', '.join(strategies)}",
        "",
    ]

    by_pri = Counter(p.priority for p in points)
    by_src = Counter(p.source for p in points)
    lines.append(f"  REQUIRED: {by_pri.get('REQUIRED', 0)}  |  "
                 f"PLIP-only: {by_src.get('plip', 0)}  |  "
                 f"Footprint-only: {by_src.get('footprint', 0)}  |  "
                 f"Both: {by_src.get('plip+footprint', 0)}")
    lines.append(f"  OPTIONAL: {by_pri.get('OPTIONAL', 0)}")
    lines.append(f"  DISABLED: {by_pri.get('DISABLED', 0)}")
    lines.append("")

    lines.append("-" * w)
    hdr = (f"{'Rk':>2s}  {'Pri':<8s}  {'Src':<10s}  {'Feature':<18s}  "
           f"{'Pocket':<10s}  {'Residues':<20s}  {'Energy':>8s}  "
           f"{'vdW':>7s}  {'ES':>7s}  {'vdW%':>5s}  {'Conf':<6s}")
    lines.append(hdr)
    lines.append("-" * w)

    for pt in points:
        star = {"REQUIRED": "***", "OPTIONAL": "** ", "DISABLED": "   "}[pt.priority]
        res_str = ", ".join(pt.residues)
        if len(res_str) > 18:
            res_str = res_str[:18] + ".."
        e_str = f"{pt.dock6_energy:+8.2f}" if pt.dock6_energy != 0 else "       -"
        vdw_str = f"{pt.dock6_vdw:+7.2f}" if pt.dock6_vdw != 0 else "      -"
        es_str = f"{pt.dock6_es:+7.2f}" if pt.dock6_es != 0 else "      -"
        vdw_pct = f"{pt.vdw_fraction * 100:4.0f}%" if pt.vdw_fraction > 0 else "    -"

        if pt.dock6_energy < -5.0 and pt.vdw_fraction > 0.6:
            conf = "HIGH"
        elif pt.dock6_energy < -1.0:
            conf = "MED"
        else:
            conf = "LOW"

        lines.append(
            f"{pt.rank:2d}  {star} {pt.priority:<8s}  {pt.source:<10s}  "
            f"{pt.pharmit_type:<18s}  {pt.sub_pocket:<10s}  {res_str:<20s}  "
            f"{e_str}  {vdw_str}  {es_str}  {vdw_pct}  {conf:<6s}"
        )

    lines.extend(["", "=" * w])
    path.write_text("\n".join(lines), encoding="utf-8")
    logger.info(f"  Saved: {path}")


# ═══════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════

def generate_pharmit_pharmacophore(
        output_dir: str,
        residue_csv_path: str,
        ligand_mol2_path: str,
        receptor_mol2_path: str,
        plip_json_path: Optional[str] = None,
        receptor_pdb_path: Optional[str] = None,
        output_name: str = "pharmacophore",
        include_hydrophobic: bool = False,
        energy_cutoff: float = 0.0,
        footprint_cutoff: float = -0.5,
        contact_distance: float = 4.5,
        n_required: int = 3,
        n_optional: int = 2,
        radius: float = 1.0,
        group_tolerance: float = 1.0,
        merge_distance: float = 2.0,
        strategies: Optional[List[str]] = None,
        ligand_name: str = "UDX",
) -> Dict[str, Any]:
    """
    Full pipeline: Footprint + PLIP → ranked sub-pocket Pharmit JSONs.

    PLIP + Footprint IS the pharmacophore. Footprint fills gaps where
    PLIP doesn't detect interactions (e.g., CH-π stacking, weak vdW).
    """
    logger.info("=" * 60)
    logger.info("PHARMIT PHARMACOPHORE GENERATOR (06a) v5.0")
    logger.info("  Source: DOCK6 Footprint + PLIP")
    logger.info("=" * 60)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if strategies is None:
        strategies = list(STRATEGIES.keys())

    # ── Step 1: PLIP → pharmacophore points (optional) ──
    logger.info("")
    logger.info("Step 1: PLIP interactions")
    plip_points = load_plip_points(
        plip_json_path,
        group_tolerance=group_tolerance,
        include_hydrophobic=include_hydrophobic,
    )

    # Residues already covered by PLIP
    plip_residues = set()
    for pt in plip_points:
        for r in pt.residues:
            plip_residues.add(r)
            plip_residues.add(r.split(".")[0])

    # ── Step 2: Footprint → new pharmacophore points ──
    logger.info("")
    logger.info("Step 2: Footprint → pharmacophore points")
    fp_points = load_footprint_as_pharmacophore(
        residue_csv_path=residue_csv_path,
        ligand_mol2_path=ligand_mol2_path,
        receptor_mol2_path=receptor_mol2_path,
        energy_cutoff=footprint_cutoff,
        contact_distance=contact_distance,
        existing_residues=plip_residues,
        receptor_pdb_path=receptor_pdb_path,
    )

    # ── Step 3: Merge ──
    logger.info("")
    logger.info("Step 3: Merge PLIP + Footprint")
    all_points = merge_points(plip_points, fp_points, merge_distance=merge_distance)

    # ── Step 4: Assign energy to PLIP-only points ──
    logger.info("")
    logger.info("Step 4: Assign energy from footprint")
    assign_energy_from_footprint(all_points, residue_csv_path)

    n_with_energy = sum(1 for p in all_points if p.dock6_energy != 0)
    logger.info(f"  {n_with_energy}/{len(all_points)} points have energy")

    # ── Step 5: Rank + filter ──
    logger.info("")
    logger.info("Step 5: Rank + filter")
    ranked = rank_and_filter(
        all_points,
        energy_cutoff=energy_cutoff,
        n_required=n_required,
        n_optional=n_optional,
    )
    logger.info(f"  {len(ranked)} features after filtering")

    if not ranked:
        return {"success": False, "error": "No pharmacophore features after filtering"}

    # Log top features
    logger.info("")
    logger.info("  Top features:")
    for pt in ranked[:10]:
        conf = "vdW" if pt.vdw_fraction > 0.6 else "ES" if pt.vdw_fraction < 0.4 else "mix"
        logger.info(f"    #{pt.rank} [{pt.priority}] {pt.pharmit_type:<18s} "
                     f"{pt.sub_pocket:<10s} {','.join(pt.residues):<15s} "
                     f"E={pt.dock6_energy:+.1f} ({conf}) [{pt.source}]")

    # ── Step 6: Generate Pharmit JSONs ──
    logger.info("")
    logger.info("Step 6: Generate Pharmit JSONs")
    generated = {}
    for strat in strategies:
        pharmit_json = generate_pharmit_json(
            ranked, strat,
            n_required=n_required, n_optional=n_optional, radius=radius,
        )
        n_pts = pharmit_json["_meta"]["n_features"]
        n_en = pharmit_json["_meta"]["n_enabled"]

        if n_pts > 0:
            json_path = out_dir / f"{output_name}_{strat}.json"
            with open(json_path, "w") as f:
                json.dump(pharmit_json, f, indent=2)
            generated[strat] = str(json_path)
            logger.info(f"  {strat}: {n_pts} features ({n_en} enabled) → {json_path.name}")
        else:
            logger.info(f"  {strat}: 0 features — skipped")

    # ── Output ranking files ──
    logger.info("")
    write_ranking_csv(ranked, out_dir / f"{output_name}_ranking.csv")
    write_ranking_txt(ranked, ligand_name, strategies, out_dir / f"{output_name}_ranking.txt")

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  {len(ranked)} features, {len(generated)} strategies")
    logger.info("=" * 60)

    return {
        "success": True,
        "n_features": len(ranked),
        "n_strategies": len(generated),
        "strategies": generated,
        "ranking_csv": str(out_dir / f"{output_name}_ranking.csv"),
        "ranking_txt": str(out_dir / f"{output_name}_ranking.txt"),
    }
