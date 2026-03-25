"""
Footprint Analysis - Core Module (04b)
=========================================
Per-residue vdW + ES energy decomposition via DOCK6 footprint scoring.

Pipeline:
    Step 1: Build residue mapping (mol2 sequential → PDB original numbering)
    Step 2: Parse per-residue vdW + ES from footprint TXT files
    Step 3: Cross-molecule consensus (which residues always contribute)
    Step 4: Compare each molecule's footprint vs reference (UDX)

DOCK6 6.13 footprint output format:
    - {name}_fps_scored.mol2           → summary scores per pose
    - {name}_fps_footprint_scored.txt  → per-residue tabular data
    - {name}_fps_hbond_scored.txt      → H-bond details

Residue numbering:
    ChimeraX mol2 renumbers residues sequentially (1..N).
    DOCK6 footprint inherits this sequential numbering.
    This module remaps to PDB original numbering (e.g., 141→392)
    so that footprint and contact results share the same residue IDs.

Input:
    01d_footprint_rescore/{name}/{name}_fps_footprint_scored.txt
    00b_receptor_preparation/rec_charged.mol2  (sequential numbering)
    00b_receptor_preparation/rec_noH.pdb       (PDB numbering)

Output:
    footprint_per_molecule.csv    — residue × molecule energy matrix (PDB numbering)
    residue_consensus.csv         — which residues always contribute
    vs_reference_comparison.csv   — delta vdW/ES vs reference per residue
    pharmacophore_residues.json   — residues contacted by >80% of molecules
    molecule_footprint_summary.csv — one row per molecule with totals
    residue_mapping.csv           — sequential→PDB mapping for reference

Location: 01_src/reference_docking/m04_dock6_analysis/footprint_analysis.py
Project: reference_docking
Module: 04b (DOCK6 analysis)
Version: 3.0 (2026-03-23) — adds sequential→PDB residue remapping

Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56
"""

import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# BINDING SITE ZONE CLASSIFICATION
# =============================================================================
# Sub-pocket classification from UDX PLIP analysis (06a).
# Used to group residues into functional zones for visualization.
# These are specific to the XT1/UDX system — future: make configurable.

ZONE_DEFINITIONS = {
    "phosphate": {
        "residues": {"ARG598", "LYS599"},
        "color": "#1D9E75",
        "label": "Phosphate (salt bridges)",
        "druglike": False,
        "description": "Requires charged groups (phosphate, carboxylate). Pharmit selects nucleotides.",
    },
    "xylose": {
        "residues": {"TRP392", "TRP495", "TYR565", "SER575"},
        "color": "#378ADD",
        "label": "Xylose pocket",
        "druglike": True,
        "description": "HBA from heterocycles (oxazol, furan, amide). Drug-like compatible.",
    },
    "ribose": {
        "residues": {"HIS335", "VAL333", "THR390"},
        "color": "#7F77DD",
        "label": "Ribose / base recognition",
        "druglike": True,
        "description": "Aromatic stacking (HIS335) + HBA/HBD. Drug-like compatible.",
    },
    "uracil": {
        "residues": {"ASP361", "ARG363"},
        "color": "#BA7517",
        "label": "Uracil pocket",
        "druglike": True,
        "description": "Weak energy but selective. Only extended molecules reach it.",
    },
    "catalytic": {
        "residues": {"GLU529"},
        "color": "#888780",
        "label": "Catalytic (GLU529)",
        "druglike": True,
        "description": "Functionally critical but low energy. Covered implicitly by xylose proximity.",
    },
}


def _classify_residue_zone(residue_name: str) -> str:
    """Classify a residue into a binding site zone."""
    for zone_id, zdef in ZONE_DEFINITIONS.items():
        if residue_name in zdef["residues"]:
            return zone_id
    return "other"


# =============================================================================
# BINDING SITE ZONES HTML REPORT
# =============================================================================

def generate_zones_html(
        df_consensus: pd.DataFrame,
        plip_json_path: Optional[str] = None,
        contact_csv_path: Optional[str] = None,
        campaign_id: str = "",
        n_molecules: int = 0,
) -> str:
    """
    Generate binding site zones HTML from residue consensus data.

    Cross-references:
        - DOCK6 footprint energy (from df_consensus — always available)
        - PLIP interactions (from 03a JSON — optional)
        - Crystal contact frequency (from PLIP — optional)

    Args:
        df_consensus: residue_consensus.csv as DataFrame
        plip_json_path: Path to 03a interactions.json (optional)
        contact_csv_path: Path to 05f contact_summary.csv (optional)
        campaign_id: Campaign name for title
        n_molecules: Number of molecules in campaign

    Returns:
        HTML string
    """
    from datetime import datetime

    # --- Load PLIP interactions (optional) ---
    plip_features = {}
    if plip_json_path and Path(plip_json_path).exists():
        with open(plip_json_path) as f:
            plip_data = json.load(f)
        for ix in plip_data.get("interactions", []):
            res_name = ix.get("residue_name", "")
            res_num = ix.get("residue_number", "")
            key = f"{res_name}{res_num}"
            itype = ix.get("type", "unknown")
            dist = ix.get("distance", 0)
            if key not in plip_features:
                plip_features[key] = []
            plip_features[key].append({"type": itype, "distance": dist})

    # --- Load contact frequency (optional, from 05f) ---
    contact_freq = {}
    n_contact_mols = 0
    if contact_csv_path and Path(contact_csv_path).exists():
        df_contacts = pd.read_csv(contact_csv_path)
        n_contact_mols = df_contacts["Name"].nunique()
        for res_id, grp in df_contacts.groupby("residue_id"):
            contact_freq[res_id] = {
                "n_molecules": grp["Name"].nunique(),
                "fraction": grp["Name"].nunique() / n_contact_mols if n_contact_mols > 0 else 0,
            }

    # --- Group residues by zone ---
    zone_data = {}
    for zone_id, zdef in ZONE_DEFINITIONS.items():
        # Match by residue_id prefix (e.g., "TRP392" from "TRP392.A")
        zone_rows = df_consensus[
            df_consensus["residue_id"].apply(
                lambda rid: rid.split(".")[0] in zdef["residues"]
            )
        ].copy()
        if zone_rows.empty:
            continue

        total_energy = zone_rows["mean_total"].sum()
        total_vdw = zone_rows["mean_vdw"].sum()
        total_es = zone_rows["mean_es"].sum()

        # Reference energy (from UDX crystal)
        ref_energy = 0
        if "ref_vdw" in zone_rows.columns and "ref_es" in zone_rows.columns:
            ref_energy = (zone_rows["ref_vdw"] + zone_rows["ref_es"]).sum()

        residues = []
        for _, row in zone_rows.iterrows():
            res_name_short = row["residue_name"]
            res_id = row["residue_id"]
            ref_tot = (row.get("ref_vdw", 0) or 0) + (row.get("ref_es", 0) or 0)

            # PLIP features for this residue (keyed by e.g. "TRP392")
            res_id_prefix = res_id.split(".")[0]
            plip_list = plip_features.get(res_id_prefix, [])

            # Contact freq from 05f
            cf = contact_freq.get(res_id, {})

            residues.append({
                "res_id": res_id,
                "res_name": res_name_short,
                "mean_total": row["mean_total"],
                "mean_vdw": row["mean_vdw"],
                "mean_es": row["mean_es"],
                "ref_total": ref_tot,
                "freq": row["frac_contributing"],
                "plip": plip_list,
                "contact_frac": cf.get("fraction", None),
                "contact_n": cf.get("n_molecules", None),
            })

        residues.sort(key=lambda r: r["mean_total"])

        zone_data[zone_id] = {
            "label": zdef["label"],
            "color": zdef["color"],
            "druglike": zdef["druglike"],
            "description": zdef["description"],
            "total_energy": round(total_energy, 2),
            "total_vdw": round(total_vdw, 2),
            "total_es": round(total_es, 2),
            "ref_energy": round(ref_energy, 2),
            "residues": residues,
        }

    # Also collect "other" residues with significant energy
    classified_ids = set()
    for zdef in ZONE_DEFINITIONS.values():
        classified_ids |= zdef["residues"]

    other_rows = df_consensus[
        (~df_consensus["residue_id"].apply(
            lambda rid: rid.split(".")[0] in classified_ids
        )) &
        (df_consensus["mean_total"] < -0.5)
    ].head(15)

    max_energy = min(z["total_energy"] for z in zone_data.values()) if zone_data else -1

    # --- Generate HTML ---
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Binding Site Zones: {campaign_id}</title>
<style>
body{{font-family:'Segoe UI',Arial,sans-serif;max-width:900px;margin:0 auto;padding:20px;background:#fafafa;color:#333}}
h1{{color:#1a5276;border-bottom:3px solid #1a5276;padding-bottom:10px;font-size:22px}}
h2{{color:#2c3e50;margin-top:30px;font-size:18px}}
.meta{{font-size:12px;color:#777;margin:5px 0 20px}}
.zone{{border:1px solid #ddd;border-radius:8px;padding:14px 16px;margin:0 0 12px;background:white}}
.zone-head{{display:flex;justify-content:space-between;align-items:center;margin:0 0 8px}}
.zone-name{{font-size:15px;font-weight:600}}
.zone-energy{{font-family:monospace;font-size:14px;font-weight:600}}
.zone-bar{{height:10px;border-radius:5px;margin:6px 0}}
.zone-desc{{font-size:12px;color:#666;margin:4px 0 8px}}
.badge{{display:inline-block;font-size:10px;padding:2px 8px;border-radius:6px;font-weight:600}}
.badge-ok{{background:#e8f8f5;color:#1e8449}}
.badge-no{{background:#fdedec;color:#c0392b}}
.badge-warn{{background:#fef9e7;color:#b7950b}}
.res-tbl{{width:100%;border-collapse:collapse;font-size:12px;margin:6px 0 0}}
.res-tbl th{{text-align:left;padding:4px 6px;color:#777;font-weight:500;border-bottom:1px solid #eee;font-size:11px}}
.res-tbl td{{padding:4px 6px;border-bottom:1px solid #f5f5f5}}
.res-tbl tr:hover{{background:#f0f7ff}}
.mono{{font-family:monospace;font-size:11px}}
.plip-tag{{display:inline-block;font-size:10px;padding:1px 5px;border-radius:3px;background:#eef2ff;color:#5b6abf;margin:1px}}
.bar-bg{{background:#f0f0f0;height:6px;border-radius:3px;display:inline-block;width:80px;vertical-align:middle}}
.bar-fill{{height:6px;border-radius:3px;display:inline-block}}
.summary{{display:flex;gap:15px;flex-wrap:wrap;margin:15px 0}}
.stat{{background:white;border:1px solid #ddd;border-radius:8px;padding:12px;min-width:130px;text-align:center}}
.stat-val{{font-size:20px;font-weight:bold;color:#2980b9}}
.stat-lbl{{font-size:10px;color:#777}}
.other-tbl{{width:100%;border-collapse:collapse;font-size:12px;margin:10px 0}}
.other-tbl th{{text-align:left;padding:4px 8px;color:#777;font-weight:500;border-bottom:1px solid #ddd;font-size:11px}}
.other-tbl td{{padding:4px 8px;border-bottom:1px solid #f0f0f0}}
.footer{{margin-top:30px;padding-top:10px;border-top:1px solid #ddd;font-size:11px;color:#999}}
</style></head><body>

<h1>Binding site zones: {campaign_id}</h1>
<p class="meta">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} | Molecules: {n_molecules}
{' | PLIP: ' + str(len(plip_features)) + ' residues' if plip_features else ''}
{' | Contacts (05f): ' + str(n_contact_mols) + ' molecules' if n_contact_mols else ''}</p>

<div class="summary">
  <div class="stat"><div class="stat-val">{len(zone_data)}</div><div class="stat-lbl">Zones defined</div></div>
  <div class="stat"><div class="stat-val">{sum(len(z['residues']) for z in zone_data.values())}</div><div class="stat-lbl">Key residues</div></div>
  <div class="stat"><div class="stat-val">{sum(1 for z in zone_data.values() if z['druglike'])}</div><div class="stat-lbl">Drug-like zones</div></div>
  <div class="stat"><div class="stat-val">{sum(z['total_energy'] for z in zone_data.values()):.1f}</div><div class="stat-lbl">Total energy</div></div>
</div>

<h2>Zones ranked by energy</h2>
"""

    for zone_id, zd in sorted(zone_data.items(), key=lambda x: x[1]["total_energy"]):
        bar_pct = min(100, abs(zd["total_energy"]) / abs(max_energy) * 100) if max_energy != 0 else 0
        dl_badge = '<span class="badge badge-ok">drug-like</span>' if zd["druglike"] else '<span class="badge badge-no">not drug-like</span>'

        html += f"""
<div class="zone" style="border-left:4px solid {zd['color']}">
  <div class="zone-head">
    <div class="zone-name" style="color:{zd['color']}">{zd['label']}</div>
    <div class="zone-energy" style="color:{zd['color']}">{zd['total_energy']:+.2f} kcal/mol</div>
  </div>
  <div class="zone-bar" style="background:#f0f0f0"><div style="width:{bar_pct:.0f}%;height:10px;border-radius:5px;background:{zd['color']};opacity:0.7"></div></div>
  <div class="zone-desc">{zd['description']} {dl_badge}
    {f'<br>Reference (UDX crystal): <b>{zd["ref_energy"]:+.2f}</b> kcal/mol' if abs(zd['ref_energy']) > 0.01 else ''}</div>
  <table class="res-tbl">
    <tr>
      <th>Residue</th>
      <th>Energy</th>
      <th>vdW</th>
      <th>ES</th>
      <th>Ref energy</th>
      <th>DOCK6 freq</th>
      {'<th>Contact freq (05f)</th>' if n_contact_mols else ''}
      {'<th>PLIP features</th>' if plip_features else ''}
    </tr>"""

        for res in zd["residues"]:
            ref_str = f"{res['ref_total']:+.2f}" if abs(res['ref_total']) > 0.01 else "—"
            freq_str = f"{res['freq']:.0%}"

            cf_str = ""
            if n_contact_mols and res['contact_frac'] is not None:
                cf_str = f"<td>{res['contact_frac']:.0%} ({res['contact_n']}/{n_contact_mols})</td>"
            elif n_contact_mols:
                cf_str = "<td style='color:#ccc'>—</td>"

            plip_str = ""
            if plip_features:
                tags = "".join(
                    f'<span class="plip-tag">{p["type"]} {p["distance"]:.1f}A</span>'
                    for p in res["plip"]
                )
                empty_plip = '<span style="color:#ccc">—</span>'
                plip_str = f"<td>{tags if tags else empty_plip}</td>"

            html += f"""
    <tr>
      <td><b>{res['res_id']}</b></td>
      <td class="mono">{res['mean_total']:+.3f}</td>
      <td class="mono">{res['mean_vdw']:+.3f}</td>
      <td class="mono">{res['mean_es']:+.3f}</td>
      <td class="mono">{ref_str}</td>
      <td>{freq_str}</td>
      {cf_str}
      {plip_str}
    </tr>"""

        html += """
  </table>
</div>"""

    # --- Other significant residues ---
    if len(other_rows) > 0:
        html += """
<h2>Other contributing residues (not in defined zones)</h2>
<p style="font-size:12px;color:#777">Residues with mean energy &lt; -0.5 kcal/mol not assigned to any zone.</p>
<table class="other-tbl">
  <tr><th>Residue</th><th>Energy</th><th>vdW</th><th>ES</th><th>Freq</th></tr>"""
        for _, row in other_rows.iterrows():
            html += f"""
  <tr>
    <td><b>{row['residue_id']}</b></td>
    <td class="mono">{row['mean_total']:+.3f}</td>
    <td class="mono">{row['mean_vdw']:+.3f}</td>
    <td class="mono">{row['mean_es']:+.3f}</td>
    <td>{row['frac_contributing']:.0%}</td>
  </tr>"""
        html += "\n</table>"

    html += f"""
<div class="footer">
  reference_docking | Module 04b | Footprint Analysis v3.1 | {datetime.now().strftime('%Y-%m-%d %H:%M')}
</div>
</body></html>"""

    return html


# =============================================================================
# RESIDUE MAPPING: mol2 sequential → PDB original
# =============================================================================

def build_residue_mapping(
        receptor_mol2: str,
        receptor_pdb: str,
) -> Dict[str, str]:
    """
    Build mapping from mol2 sequential numbering to PDB original numbering.

    ChimeraX mol2 SUBSTRUCTURE section lists residues sequentially (1..N).
    PDB CA atoms list residues with original numbering + chain.
    Both have the same number of residues in the same order.

    Returns:
        Dict mapping "RESseq" → "RES_pdbnum.chain"
        e.g., {"TRP141" → "TRP392.A", "HIS84" → "HIS335.A"}
    """
    # --- Parse mol2 SUBSTRUCTURE ---
    mol2_residues = []  # list of (seq_id, resname)
    with open(receptor_mol2, "r") as f:
        in_substructure = False
        for line in f:
            if "@<TRIPOS>SUBSTRUCTURE" in line:
                in_substructure = True
                continue
            if in_substructure:
                if line.startswith("@"):
                    break
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        seq_id = int(parts[0])
                        resname = parts[1]
                        mol2_residues.append((seq_id, resname))
                    except (ValueError, IndexError):
                        continue

    # --- Parse PDB CA atoms ---
    pdb_residues = []  # list of (resname, resid, chain)
    with open(receptor_pdb, "r") as f:
        for line in f:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and " CA " in line:
                resname = line[17:20].strip()
                chain = line[21].strip() or "A"
                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    continue
                pdb_residues.append((resname, resid, chain))

    # --- Build mapping ---
    if len(mol2_residues) != len(pdb_residues):
        logger.warning(f"  Residue count mismatch: mol2={len(mol2_residues)}, PDB={len(pdb_residues)}")
        logger.warning("  Falling back to sequential numbering (no remapping)")
        return {}

    mapping = {}
    mismatches = 0
    for (seq_id, mol2_name), (pdb_name, pdb_resid, chain) in zip(mol2_residues, pdb_residues):
        # Verify residue names match
        if mol2_name[:3].upper() != pdb_name[:3].upper():
            mismatches += 1
            if mismatches <= 3:
                logger.warning(f"  Residue name mismatch at seq={seq_id}: "
                               f"mol2={mol2_name}, PDB={pdb_name}{pdb_resid}.{chain}")

        seq_key = f"{mol2_name}{seq_id}"
        pdb_key = f"{pdb_name}{pdb_resid}.{chain}"
        mapping[seq_key] = pdb_key

    if mismatches > 3:
        logger.warning(f"  ... and {mismatches - 3} more mismatches")
    if mismatches > len(mol2_residues) * 0.1:
        logger.error(f"  Too many mismatches ({mismatches}/{len(mol2_residues)}), disabling remapping")
        return {}

    logger.info(f"  Residue mapping: {len(mapping)} residues (mol2 sequential → PDB original)")
    return mapping


# =============================================================================
# FOOTPRINT TXT PARSER (per-residue tabular data)
# =============================================================================

def parse_footprint_txt(
        txt_path: str,
        residue_mapping: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """
    Parse a DOCK6 footprint TXT file (*_fps_footprint_scored.txt).

    Applies residue_mapping to convert sequential→PDB numbering if provided.

    Returns:
        List of pose dicts, each with:
            - pose_id (int)
            - fps_score, fps_vdw_energy, fps_es_energy (float or None)
            - residue_footprint: list of per-residue dicts
    """
    with open(txt_path, "r") as f:
        content = f.read()

    # Split into pose blocks by the "### Molecule:" separator
    blocks = re.split(r"#{20,}\s*\n\s*###\s+Molecule:", content)

    poses = []

    for block_idx, block in enumerate(blocks):
        # Skip preamble before first molecule
        if block_idx == 0 and "resname" not in block:
            continue

        lines = block.strip().split("\n")

        # --- Parse summary fields from ## lines ---
        fps_score = None
        fps_vdw_energy = None
        fps_es_energy = None
        fps_vdw_es_energy = None
        fps_num_hbond = None

        for line in lines:
            line_s = line.strip()
            if "Footprint_Similarity_Score:" in line_s:
                try:
                    fps_score = float(line_s.split(":")[-1].strip())
                except ValueError:
                    pass
            elif "FPS_vdw_energy:" in line_s:
                try:
                    fps_vdw_energy = float(line_s.split(":")[-1].strip())
                except ValueError:
                    pass
            elif "FPS_es_energy:" in line_s:
                try:
                    fps_es_energy = float(line_s.split(":")[-1].strip())
                except ValueError:
                    pass
            elif "FPS_vdw+es_energy:" in line_s:
                try:
                    fps_vdw_es_energy = float(line_s.split(":")[-1].strip())
                except ValueError:
                    pass
            elif "FPS_num_hbond:" in line_s:
                try:
                    fps_num_hbond = int(line_s.split(":")[-1].strip())
                except ValueError:
                    pass

        # --- Parse per-residue tabular data ---
        residue_footprint = []
        in_table = False

        for line in lines:
            line_s = line.strip()

            # Detect header row
            if line_s.startswith("resname") and "resid" in line_s:
                in_table = True
                continue

            # Detect end of table (empty line or new ## block)
            if in_table and (not line_s or line_s.startswith("#")):
                in_table = False
                continue

            if in_table:
                parts = line_s.split()
                if len(parts) >= 8:
                    try:
                        resname = parts[0]
                        resid = int(parts[1])
                        vdw_ref = float(parts[2])
                        es_ref = float(parts[3])
                        hb_ref = int(float(parts[4]))
                        vdw_pose = float(parts[5])
                        es_pose = float(parts[6])
                        hb_pose = int(float(parts[7]))

                        # Sequential key (from DOCK6 output)
                        seq_key = f"{resname}{resid}"

                        # Remap to PDB numbering if mapping available
                        if residue_mapping and seq_key in residue_mapping:
                            residue_id = residue_mapping[seq_key]
                            # Parse PDB residue_id back to components
                            m_pdb = re.match(r"([A-Z]{1,4})(\d+)\.(\w+)", residue_id)
                            if m_pdb:
                                resname_out = m_pdb.group(1)
                                resid_out = int(m_pdb.group(2))
                                chain_out = m_pdb.group(3)
                            else:
                                resname_out = resname
                                resid_out = resid
                                chain_out = "A"
                                residue_id = f"{resname}{resid}.A"
                        else:
                            resname_out = resname
                            resid_out = resid
                            chain_out = "A"
                            residue_id = f"{resname}{resid}.A"

                        total_pose = vdw_pose + es_pose
                        total_ref = vdw_ref + es_ref

                        residue_footprint.append({
                            "residue_id": residue_id,
                            "residue_name": resname_out,
                            "residue_number": resid_out,
                            "chain": chain_out,
                            "vdw": round(vdw_pose, 6),
                            "es": round(es_pose, 6),
                            "total": round(total_pose, 6),
                            "ref_vdw": round(vdw_ref, 6),
                            "ref_es": round(es_ref, 6),
                            "ref_total": round(total_ref, 6),
                            "delta_vdw": round(vdw_pose - vdw_ref, 6),
                            "delta_es": round(es_pose - es_ref, 6),
                            "delta_total": round(total_pose - total_ref, 6),
                            "hb_pose": hb_pose,
                            "hb_ref": hb_ref,
                        })
                    except (ValueError, IndexError) as e:
                        logger.debug(f"  Skipping malformed line: {line_s} ({e})")
                        continue

        pose_id = block_idx if block_idx > 0 else 0
        poses.append({
            "pose_id": pose_id,
            "fps_score": fps_score,
            "fps_vdw_energy": fps_vdw_energy,
            "fps_es_energy": fps_es_energy,
            "fps_vdw_es_energy": fps_vdw_es_energy,
            "fps_num_hbond": fps_num_hbond,
            "residue_footprint": residue_footprint,
            "n_residues": len(residue_footprint),
        })

    return poses


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_footprint_analysis(
        footprint_dir: Union[str, Path],
        output_dir: Union[str, Path],
        receptor_mol2: Optional[str] = None,
        receptor_pdb: Optional[str] = None,
        pharmacophore_threshold: float = 0.8,
        energy_cutoff: float = -0.5,
        best_pose_only: bool = True,
        plip_json: Optional[str] = None,
        contact_csv: Optional[str] = None,
        campaign_id: str = "",
) -> Dict[str, Any]:
    """
    Analyze DOCK6 footprint re-scoring results.

    Args:
        footprint_dir:  Directory with {name}/{name}_fps_footprint_scored.txt
        output_dir:     Output directory for analysis files
        receptor_mol2:  Path to rec_charged.mol2 (sequential numbering, for mapping)
        receptor_pdb:   Path to rec_noH.pdb (PDB original numbering, for mapping)
        pharmacophore_threshold: Fraction of molecules that must contact a residue
        energy_cutoff:  Minimum total energy for a residue to count as "contributing"
        best_pose_only: If True, analyze only best pose per molecule

    Returns:
        Dict with: success, n_molecules, output paths
    """
    footprint_dir = Path(footprint_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  04b DOCK6 Footprint Analysis v3.0")
    logger.info("=" * 60)

    # --- Build residue mapping ---
    residue_mapping = {}
    if receptor_mol2 and receptor_pdb and Path(receptor_mol2).exists() and Path(receptor_pdb).exists():
        logger.info(f"  Building residue mapping (mol2 → PDB)...")
        residue_mapping = build_residue_mapping(receptor_mol2, receptor_pdb)
        if residue_mapping:
            # Save mapping for reference
            mapping_csv = output_dir / "residue_mapping.csv"
            mapping_rows = []
            for seq_key, pdb_key in sorted(residue_mapping.items(),
                                            key=lambda x: int(re.search(r'\d+', x[0]).group())):
                mapping_rows.append({"mol2_sequential": seq_key, "pdb_original": pdb_key})
            pd.DataFrame(mapping_rows).to_csv(mapping_csv, index=False)
            logger.info(f"  Saved: {mapping_csv}")
    else:
        logger.warning("  No receptor mol2/PDB provided — using sequential numbering")
        logger.warning("  (Residue IDs will NOT match contact mapping)")

    # --- Find molecule directories ---
    mol_dirs = sorted([
        d for d in footprint_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    ])

    all_rows = []       # flat table: Name × residue × energy
    mol_summaries = []  # one row per molecule
    n_parsed = 0
    n_failed = 0

    for d in mol_dirs:
        name = d.name

        # --- Find footprint TXT file ---
        fps_txt = d / f"{name}_fps_footprint_scored.txt"
        if not fps_txt.exists():
            txt_candidates = list(d.glob("*_fps_footprint_scored.txt"))
            if txt_candidates:
                fps_txt = txt_candidates[0]
            else:
                logger.debug(f"  No footprint TXT found for {name}, skipping")
                n_failed += 1
                continue

        poses = parse_footprint_txt(str(fps_txt), residue_mapping=residue_mapping)
        if not poses:
            logger.warning(f"  No poses parsed from TXT: {name}")
            n_failed += 1
            continue

        # Filter out empty poses
        poses = [p for p in poses if p["residue_footprint"]]
        if not poses:
            logger.warning(f"  No residue data in any pose: {name}")
            n_failed += 1
            continue

        # Select best pose (lowest FPS_Score = most similar to reference)
        if best_pose_only and len(poses) > 1:
            scored = [p for p in poses if p["fps_score"] is not None]
            if scored:
                poses = [min(scored, key=lambda p: p["fps_score"])]
            else:
                poses = [poses[0]]

        for pose in poses:
            fps_score = pose.get("fps_score")
            fps_vdw = pose.get("fps_vdw_energy")
            fps_es = pose.get("fps_es_energy")

            for res in pose["residue_footprint"]:
                all_rows.append({
                    "Name": name,
                    "pose_id": pose["pose_id"],
                    "fps_score": fps_score,
                    **res,
                })

            # Summary per molecule
            fp = pose["residue_footprint"]
            total_vdw = sum(r["vdw"] for r in fp)
            total_es = sum(r["es"] for r in fp)
            n_contributing = sum(1 for r in fp if r["total"] < energy_cutoff)

            mol_summaries.append({
                "Name": name,
                "pose_id": pose["pose_id"],
                "fps_score": fps_score,
                "fps_vdw_energy": fps_vdw,
                "fps_es_energy": fps_es,
                "fps_vdw_es_energy": pose.get("fps_vdw_es_energy"),
                "fps_num_hbond": pose.get("fps_num_hbond"),
                "total_vdw": round(total_vdw, 3),
                "total_es": round(total_es, 3),
                "total_energy": round(total_vdw + total_es, 3),
                "n_residues_total": len(fp),
                "n_residues_contributing": n_contributing,
            })

        n_parsed += 1

    logger.info(f"  Parsed: {n_parsed} molecules, {n_failed} failed")

    if not all_rows:
        logger.warning("  No footprint data found in any TXT file!")
        logger.warning("  Possible causes:")
        logger.warning("    1. Footprint re-scoring hasn't been run (run 01d first)")
        logger.warning("    2. TXT files not generated (check write_footprints=yes)")
        return {
            "success": False,
            "error": "No footprint data found. Run 01d footprint re-scoring first.",
            "n_parsed": n_parsed,
        }

    # --- Build DataFrames ---
    df_all = pd.DataFrame(all_rows)
    df_summary = pd.DataFrame(mol_summaries)

    # --- Save footprint_per_molecule.csv ---
    fps_csv = output_dir / "footprint_per_molecule.csv"
    df_all.to_csv(fps_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {fps_csv} ({len(df_all)} rows)")

    # --- Residue consensus ---
    residue_stats = []
    all_names = df_all["Name"].nunique()

    for res_id, grp in df_all.groupby("residue_id"):
        n_mol = grp["Name"].nunique()
        n_contributing = grp[grp["total"] < energy_cutoff]["Name"].nunique()

        residue_stats.append({
            "residue_id": res_id,
            "residue_name": grp.iloc[0]["residue_name"],
            "residue_number": grp.iloc[0]["residue_number"],
            "chain": grp.iloc[0]["chain"],
            "n_molecules_present": n_mol,
            "n_molecules_contributing": n_contributing,
            "frac_present": round(n_mol / all_names, 3) if all_names > 0 else 0,
            "frac_contributing": round(n_contributing / all_names, 3) if all_names > 0 else 0,
            "mean_vdw": round(grp["vdw"].mean(), 4),
            "mean_es": round(grp["es"].mean(), 4),
            "mean_total": round(grp["total"].mean(), 4),
            "std_total": round(grp["total"].std(), 4) if len(grp) > 1 else 0.0,
            "min_total": round(grp["total"].min(), 4),
            "max_total": round(grp["total"].max(), 4),
            "ref_vdw": round(grp["ref_vdw"].iloc[0], 4) if "ref_vdw" in grp.columns else 0.0,
            "ref_es": round(grp["ref_es"].iloc[0], 4) if "ref_es" in grp.columns else 0.0,
        })

    df_consensus = pd.DataFrame(residue_stats)
    df_consensus.sort_values("mean_total", ascending=True, inplace=True)
    df_consensus.reset_index(drop=True, inplace=True)

    consensus_csv = output_dir / "residue_consensus.csv"
    df_consensus.to_csv(consensus_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {consensus_csv} ({len(df_consensus)} residues)")

    # --- Pharmacophore residues ---
    pharma = df_consensus[
        df_consensus["frac_contributing"] >= pharmacophore_threshold
    ].copy()

    pharma_list = []
    for _, row in pharma.iterrows():
        pharma_list.append({
            "residue_id": row["residue_id"],
            "residue_name": row["residue_name"],
            "residue_number": int(row["residue_number"]),
            "chain": row["chain"],
            "frac_contributing": row["frac_contributing"],
            "mean_total": row["mean_total"],
        })

    pharma_json = output_dir / "pharmacophore_residues.json"
    with open(pharma_json, "w") as f:
        json.dump({
            "threshold": pharmacophore_threshold,
            "energy_cutoff": energy_cutoff,
            "n_molecules": all_names,
            "n_pharmacophore_residues": len(pharma_list),
            "numbering": "PDB_original" if residue_mapping else "mol2_sequential",
            "residues": pharma_list,
        }, f, indent=2)
    logger.info(f"  Saved: {pharma_json} ({len(pharma_list)} pharmacophore residues)")

    # --- vs reference comparison ---
    ref_csv = None
    if "delta_total" in df_all.columns:
        ref_comparison = []
        for name, grp in df_all.groupby("Name"):
            sig = grp[grp["delta_total"].abs() > 0.5]
            for _, row in sig.iterrows():
                ref_comparison.append({
                    "Name": name,
                    "residue_id": row["residue_id"],
                    "vdw": row["vdw"],
                    "es": row["es"],
                    "ref_vdw": row["ref_vdw"],
                    "ref_es": row["ref_es"],
                    "delta_vdw": row["delta_vdw"],
                    "delta_es": row["delta_es"],
                    "delta_total": row["delta_total"],
                })

        if ref_comparison:
            df_ref = pd.DataFrame(ref_comparison)
            ref_csv = output_dir / "vs_reference_comparison.csv"
            df_ref.to_csv(ref_csv, index=False, encoding="utf-8")
            logger.info(f"  Saved: {ref_csv} ({len(df_ref)} significant deltas)")
        else:
            logger.info("  No significant deltas vs reference found.")

    # --- Save molecule summary ---
    summary_csv = output_dir / "molecule_footprint_summary.csv"
    df_summary.to_csv(summary_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {summary_csv}")

    # --- Log summary ---
    logger.info("")
    logger.info(f"  Molecules analyzed:      {n_parsed}")
    logger.info(f"  Total residues tracked:  {len(df_consensus)}")
    logger.info(f"  Pharmacophore residues:  {len(pharma_list)} (>{pharmacophore_threshold * 100:.0f}%)")
    logger.info(f"  Numbering:               {'PDB original' if residue_mapping else 'mol2 sequential'}")
    if len(pharma_list) > 0:
        top5 = pharma_list[:5]
        for r in top5:
            logger.info(f"    {r['residue_id']}: {r['frac_contributing'] * 100:.0f}% "
                        f"(mean={r['mean_total']:.2f} kcal/mol)")
    logger.info("=" * 60)

    # --- Generate binding site zones HTML ---
    zones_html_path = None
    try:
        html = generate_zones_html(
            df_consensus=df_consensus,
            plip_json_path=plip_json,
            contact_csv_path=contact_csv,
            campaign_id=campaign_id,
            n_molecules=n_parsed,
        )
        zones_html_path = output_dir / "binding_site_zones.html"
        zones_html_path.write_text(html, encoding="utf-8")
        logger.info(f"  Zones report: {zones_html_path}")
    except Exception as e:
        logger.warning(f"  Zones HTML generation failed: {e}")

    return {
        "success": True,
        "n_molecules": n_parsed,
        "n_residues": len(df_consensus),
        "n_pharmacophore": len(pharma_list),
        "numbering": "PDB_original" if residue_mapping else "mol2_sequential",
        "footprint_per_molecule_csv": str(fps_csv),
        "residue_consensus_csv": str(consensus_csv),
        "pharmacophore_json": str(pharma_json),
        "vs_reference_csv": str(ref_csv) if ref_csv else None,
        "molecule_summary_csv": str(summary_csv),
        "zones_html": str(zones_html_path) if zones_html_path else None,
        "output_dir": str(output_dir),
    }