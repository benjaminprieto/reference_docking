"""
Pharmit Zone Selector - Core Module (06b)
============================================
Reads a Pharmit-exported JSON (reference/pharmit.json), maps each
auto-detected feature to a binding site zone using proximity to
LIGAND ATOMS, annotates with receptor-side evidence (PLIP + footprint),
and generates variant JSONs with different zone combinations.

v3.1: Zone by ligand atom + evidence annotation.
Zone classification uses ligand atom position (v3.0).
Evidence annotation uses PLIP interactions + DOCK6 footprint energy
to provide confidence indicators for each feature.

This allows the user to make informed decisions in Pharmit:
  - Features with HIGH confidence (strong energy, PLIP confirmed) → keep
  - Features with LOW confidence (no energy support) → try removing first
  - Features with NONE (no receptor evidence) → experimental

Ligand zones (UDX atom names):
    xylose:    C1', C2', O2', C3', O3', C4', O4', C5', O5'
    phosphate: PA, PB, O1A, O2A, O3A, O1B, O2B, O3B
    ribose:    C1D, C2D, O2D, C3D, O3D, C4D, O4D, C5D, O5D
    uracil:    N1, C2, O2, N3, C4, O4, C5, C6

Input:
    Required: reference/pharmit.json (from pharmit.csb.pitt.edu)
    Optional: PLIP interactions.json (from 03a)
    Optional: residue_consensus.csv (from reference_docking 04b)

Output:
    pharmit_{strategy}.json      — one per zone combination
    pharmit_zone_mapping.csv     — feature→zone + evidence
    pharmit_zone_mapping.html    — visual decision table with energy

Location: 01_src/reference_docking/m06_pharmit/pharmit_zone_selector.py
Project: reference_docking
Module: 06b
Version: 3.1 — Evidence-annotated (2026-03-26)
"""

import csv
import json
import logging
import math
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# ZONE DEFINITIONS — by LIGAND atom names
# =============================================================================

ATOM_ZONE_MAP = {
    # Xylose sugar
    "C1'": "xylose", "C2'": "xylose", "O2'": "xylose",
    "C3'": "xylose", "O3'": "xylose", "C4'": "xylose",
    "O4'": "xylose", "C5'": "xylose", "O5'": "xylose",
    # Phosphates
    "PA": "phosphate", "PB": "phosphate",
    "O1A": "phosphate", "O2A": "phosphate", "O3A": "phosphate",
    "O1B": "phosphate", "O2B": "phosphate", "O3B": "phosphate",
    # Ribose
    "C1D": "ribose", "C2D": "ribose", "O2D": "ribose",
    "C3D": "ribose", "O3D": "ribose", "C4D": "ribose",
    "O4D": "ribose", "C5D": "ribose", "O5D": "ribose",
    # Uracil
    "N1": "uracil", "C2": "uracil", "O2": "uracil",
    "N3": "uracil", "C4": "uracil", "O4": "uracil",
    "C5": "uracil", "C6": "uracil",
}

ZONE_PROPERTIES = {
    "xylose": {"color": "#378ADD", "druglike": True, "label": "Xylose"},
    "ribose": {"color": "#7F77DD", "druglike": True, "label": "Ribose"},
    "uracil": {"color": "#BA7517", "druglike": True, "label": "Uracil"},
    "phosphate": {"color": "#1D9E75", "druglike": False, "label": "Phosphate"},
    "catalytic": {"color": "#888780", "druglike": True, "label": "Catalytic"},
    "other": {"color": "#ccc", "druglike": False, "label": "Other"},
    "structural": {"color": "#ccc", "druglike": False, "label": "Structural"},
}

ZONE_STRATEGIES = {
    "xylose_only": {
        "zones": ["xylose"],
        "description": "Xylose pocket only (sugar OH groups)",
    },
    "xylose_ribose": {
        "zones": ["xylose", "ribose"],
        "description": "Xylose + ribose (sugar + nucleoside)",
    },
    "xylose_ribose_uracil": {
        "zones": ["xylose", "ribose", "uracil"],
        "description": "All drug-like zones (no phosphate)",
    },
    "ribose_only": {
        "zones": ["ribose"],
        "description": "Ribose/base only (nucleoside pocket)",
    },
    "all_druglike": {
        "zones": ["xylose", "ribose", "uracil"],
        "description": "All drug-like zones",
    },
}

# Receptor residue → zone (for evidence annotation, NOT for classification)
RECEPTOR_ZONE_MAP = {
    "TRP392": "xylose", "TRP495": "xylose", "TYR565": "xylose",
    "SER575": "xylose", "ASP494": "xylose", "GLU529": "catalytic",
    "HIS335": "ribose", "VAL333": "ribose", "THR390": "ribose",
    "ASP361": "uracil", "ARG363": "uracil",
    "ARG598": "phosphate", "LYS599": "phosphate",
}


# =============================================================================
# STEP 1: PARSE LIGAND ATOMS FROM PHARMIT JSON
# =============================================================================

def _parse_ligand_from_pharmit(pharmit_data: Dict) -> List[Dict]:
    """Extract ligand atom coordinates from embedded HETATM records."""
    ligand_str = pharmit_data.get("ligand", "")
    if not ligand_str:
        return []

    atoms = []
    for line in ligand_str.split("\n"):
        line = line.strip()
        if not line.startswith("HETATM"):
            continue
        try:
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            zone = ATOM_ZONE_MAP.get(atom_name, "other")
            atoms.append({"atom_name": atom_name, "x": x, "y": y, "z": z, "zone": zone})
        except (ValueError, IndexError):
            continue

    logger.info(f"  Parsed {len(atoms)} ligand atoms")
    return atoms


# =============================================================================
# STEP 1b: LOAD EVIDENCE SOURCES (for annotation, not classification)
# =============================================================================

def _load_plip_anchors(plip_json_path: Optional[str]) -> List[Dict]:
    """Load PLIP interactions as receptor-side anchors for evidence."""
    if not plip_json_path or not Path(plip_json_path).exists():
        return []

    with open(plip_json_path) as f:
        data = json.load(f)

    anchors = []
    for ix in data.get("interactions", []):
        res = ix.get("residue", "")
        chain = ix.get("chain", "A")
        lig_coords = ix.get("ligand_coords")
        if lig_coords and len(lig_coords) >= 3:
            anchors.append({
                "x": lig_coords[0], "y": lig_coords[1], "z": lig_coords[2],
                "residue": res,
                "type": ix.get("interaction_type", ""),
            })

    logger.info(f"  PLIP: {len(anchors)} interaction anchors")
    return anchors


def _load_footprint_energy(footprint_csv_path: Optional[str]) -> Dict[str, Dict]:
    """Load footprint per-residue energy as evidence."""
    if not footprint_csv_path or not Path(footprint_csv_path).exists():
        return {}

    df = pd.read_csv(footprint_csv_path)

    # Detect columns
    rid_col = "residue_id" if "residue_id" in df.columns else df.columns[0]
    vdw_col = "mean_vdw" if "mean_vdw" in df.columns else ("ref_vdw" if "ref_vdw" in df.columns else None)
    es_col = "mean_es" if "mean_es" in df.columns else ("ref_es" if "ref_es" in df.columns else None)
    total_col = "mean_total" if "mean_total" in df.columns else None

    energy = {}
    for _, row in df.iterrows():
        rid = str(row[rid_col]).split(".")[0]  # Strip chain: ARG598.A → ARG598
        vdw = float(row[vdw_col]) if vdw_col else 0.0
        es = float(row[es_col]) if es_col else 0.0
        total = float(row[total_col]) if total_col else vdw + es
        energy[rid] = {"vdw": round(vdw, 2), "es": round(es, 2), "total": round(total, 2)}

    logger.info(f"  Footprint: {len(energy)} residues with energy")
    return energy


# =============================================================================
# STEP 2: MAP FEATURES TO ZONES + ANNOTATE EVIDENCE
# =============================================================================

def _dist3d(a: Dict, b: Dict) -> float:
    return math.sqrt((a["x"] - b["x"])**2 + (a["y"] - b["y"])**2 + (a["z"] - b["z"])**2)


def map_features_to_zones(
        pharmit_points: List[Dict],
        ligand_atoms: List[Dict],
        plip_anchors: List[Dict],
        footprint_energy: Dict[str, Dict],
        max_distance: float = 2.0,
        plip_max_distance: float = 3.5,
) -> List[Dict]:
    """
    Map features to zones (by ligand atom) and annotate with evidence
    (nearest receptor residue + footprint energy).
    """
    mappings = []

    for i, pt in enumerate(pharmit_points):
        name = pt["name"]
        x, y, z = pt["x"], pt["y"], pt["z"]
        enabled = pt.get("enabled", False)

        if name == "InclusionSphere":
            mappings.append({
                "index": i, "name": name, "x": x, "y": y, "z": z,
                "enabled": enabled, "zone": "structural",
                "ligand_atom": "—", "dist": 0, "decision": "KEEP_OFF",
                "receptor_residue": "—", "plip_type": "—",
                "energy": 0, "vdw": 0, "es": 0, "confidence": "—",
            })
            continue

        # --- Zone from ligand atom ---
        best_dist = 999.0
        best_atom = None
        for atom in ligand_atoms:
            d = _dist3d({"x": x, "y": y, "z": z}, atom)
            if d < best_dist:
                best_dist = d
                best_atom = atom

        if best_atom and best_dist <= max_distance:
            zone = best_atom["zone"]
            atom_name = best_atom["atom_name"]
        else:
            zone = "other"
            atom_name = f"({best_atom['atom_name']}@{best_dist:.1f})" if best_atom else "—"

        # --- Evidence from PLIP ---
        plip_residue = "—"
        plip_type = "—"
        plip_dist = 999.0
        for anchor in plip_anchors:
            d = _dist3d({"x": x, "y": y, "z": z}, anchor)
            if d < plip_dist:
                plip_dist = d
                plip_residue = anchor["residue"]
                plip_type = anchor["type"]
        if plip_dist > plip_max_distance:
            plip_residue = "—"
            plip_type = "—"

        # --- Evidence from footprint ---
        energy = 0.0
        vdw = 0.0
        es = 0.0
        if plip_residue != "—" and plip_residue in footprint_energy:
            e = footprint_energy[plip_residue]
            energy = e["total"]
            vdw = e["vdw"]
            es = e["es"]

        # --- Confidence ---
        if energy < -5.0:
            confidence = "HIGH"
        elif energy < -1.0:
            confidence = "MEDIUM"
        elif energy < -0.1:
            confidence = "LOW"
        elif plip_residue != "—":
            confidence = "PLIP_ONLY"
        else:
            confidence = "NONE"

        # --- Decision ---
        zone_props = ZONE_PROPERTIES.get(zone, {"druglike": False})
        decision = "ZONE_ACTIVE" if zone_props["druglike"] else "DISABLE"

        mappings.append({
            "index": i, "name": name, "x": x, "y": y, "z": z,
            "enabled": enabled, "zone": zone,
            "ligand_atom": atom_name, "dist": round(best_dist, 2),
            "decision": decision,
            "receptor_residue": plip_residue, "plip_type": plip_type,
            "energy": round(energy, 2), "vdw": round(vdw, 2),
            "es": round(es, 2), "confidence": confidence,
        })

    # Log summary
    zone_counts = {}
    for m in mappings:
        z = m["zone"]
        zone_counts[z] = zone_counts.get(z, 0) + 1
    for z, n in sorted(zone_counts.items()):
        logger.info(f"    {z}: {n} features")

    return mappings


# =============================================================================
# STEP 3: GENERATE VARIANT JSONS
# =============================================================================

def generate_variants(
        pharmit_data: Dict,
        mappings: List[Dict],
        output_dir: Path,
        strategies: Optional[Dict] = None,
) -> Dict[str, str]:
    """Generate variant JSONs. Only modifies points[].enabled."""
    if strategies is None:
        strategies = ZONE_STRATEGIES

    outputs = {}
    for strat_name, strat_def in strategies.items():
        active_zones = set(strat_def["zones"])
        variant = deepcopy(pharmit_data)

        n_on = 0
        n_off = 0
        for mapping in mappings:
            idx = mapping["index"]
            pt = variant["points"][idx]

            if mapping["zone"] == "structural":
                continue

            if mapping["zone"] in active_zones and mapping["decision"] == "ZONE_ACTIVE":
                pt["enabled"] = True
                n_on += 1
            else:
                pt["enabled"] = False
                n_off += 1

        out_path = output_dir / f"pharmit_{strat_name}.json"
        with open(out_path, "w") as f:
            json.dump(variant, f, indent=4)

        outputs[strat_name] = str(out_path)
        logger.info(f"  {strat_name}: {n_on} ON, {n_off} OFF → {out_path.name}")

    return outputs


# =============================================================================
# STEP 4: MAPPING TABLE (CSV + HTML with evidence)
# =============================================================================

def write_mapping_csv(mappings: List[Dict], path: Path):
    """Save feature→zone mapping with evidence as CSV."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "index", "pharmit_type", "x", "y", "z",
            "pharmit_enabled", "zone", "ligand_atom", "dist_A",
            "receptor_residue", "plip_type",
            "energy", "vdw", "es", "confidence", "decision",
        ])
        for m in mappings:
            writer.writerow([
                m["index"], m["name"], m["x"], m["y"], m["z"],
                m["enabled"], m["zone"], m["ligand_atom"], m["dist"],
                m["receptor_residue"], m["plip_type"],
                m["energy"], m["vdw"], m["es"],
                m["confidence"], m["decision"],
            ])
    logger.info(f"  Saved: {path}")


def write_mapping_html(
        mappings: List[Dict],
        strategies: Dict,
        outputs: Dict[str, str],
        campaign_id: str,
        path: Path,
):
    """Generate HTML with zone mapping + energy evidence."""

    zone_colors = {zid: zdef["color"] for zid, zdef in ZONE_PROPERTIES.items()}

    # Confidence badge styles
    conf_badge = {
        "HIGH": '<span style="background:#e8f8f5;color:#1e8449;padding:1px 6px;border-radius:4px;font-size:10px;font-weight:600">HIGH</span>',
        "MEDIUM": '<span style="background:#fef9e7;color:#b7950b;padding:1px 6px;border-radius:4px;font-size:10px;font-weight:600">MED</span>',
        "LOW": '<span style="background:#fdedec;color:#c0392b;padding:1px 6px;border-radius:4px;font-size:10px;font-weight:600">LOW</span>',
        "PLIP_ONLY": '<span style="background:#eef2ff;color:#5b6abf;padding:1px 6px;border-radius:4px;font-size:10px;font-weight:600">PLIP</span>',
        "NONE": '<span style="background:#f0f0f0;color:#999;padding:1px 6px;border-radius:4px;font-size:10px">NONE</span>',
        "—": "",
    }

    rows_html = ""
    for m in mappings:
        color = zone_colors.get(m["zone"], "#ccc")
        dec_class = "act" if m["decision"] == "ZONE_ACTIVE" else "skip"

        # Energy cell: colored by value
        e = m["energy"]
        if e < -5.0:
            e_style = "color:#1e8449;font-weight:600"
        elif e < -1.0:
            e_style = "color:#b7950b"
        elif e < 0:
            e_style = "color:#c0392b"
        elif e == 0:
            e_style = "color:#ccc"
        else:
            e_style = "color:#c0392b;font-weight:600"

        e_str = f"{e:+.1f}" if e != 0 else "—"
        vdw_str = f"{m['vdw']:+.1f}" if m['vdw'] != 0 else ""
        es_str = f"{m['es']:+.1f}" if m['es'] != 0 else ""

        badge = conf_badge.get(m["confidence"], "")

        rows_html += f"""<tr class="{dec_class}">
  <td>{m['index']+1}</td>
  <td>{m['name']}</td>
  <td>{'ON' if m['enabled'] else 'off'}</td>
  <td style="border-left:3px solid {color};padding-left:8px"><b>{m['zone']}</b></td>
  <td class="mono">{m['ligand_atom']}</td>
  <td class="mono" style="font-size:10px;color:#555">({m['x']:.1f}, {m['y']:.1f}, {m['z']:.1f})</td>
  <td class="mono">{m['receptor_residue']}</td>
  <td class="mono">{m['plip_type']}</td>
  <td class="mono" style="{e_style}">{e_str}</td>
  <td class="mono" style="font-size:10px;color:#777">{vdw_str}</td>
  <td class="mono" style="font-size:10px;color:#777">{es_str}</td>
  <td>{badge}</td>
  <td>{m['decision']}</td>
</tr>"""

    strat_rows = ""
    for sname, sdef in strategies.items():
        zones_str = " + ".join(sdef["zones"])
        fname = Path(outputs.get(sname, "")).name if sname in outputs else "—"
        strat_rows += (f"<tr><td><b>{sname}</b></td><td>{zones_str}</td>"
                       f"<td>{sdef['description']}</td><td class='mono'>{fname}</td></tr>")

    # Count summary
    n_high = sum(1 for m in mappings if m["confidence"] == "HIGH")
    n_med = sum(1 for m in mappings if m["confidence"] == "MEDIUM")
    n_low = sum(1 for m in mappings if m["confidence"] == "LOW")
    n_plip = sum(1 for m in mappings if m["confidence"] == "PLIP_ONLY")
    n_none = sum(1 for m in mappings if m["confidence"] == "NONE")
    n_active = sum(1 for m in mappings if m["decision"] == "ZONE_ACTIVE")

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Pharmit Zone Mapping: {campaign_id}</title>
<style>
body{{font-family:'Segoe UI',Arial,sans-serif;max-width:1250px;margin:0 auto;padding:20px;background:#fafafa;color:#333}}
h1{{color:#1a5276;border-bottom:3px solid #1a5276;padding-bottom:10px;font-size:22px}}
h2{{color:#2c3e50;margin-top:30px;font-size:18px}}
.tbl{{width:100%;border-collapse:collapse;font-size:12px;margin:10px 0}}
.tbl th{{text-align:left;padding:5px 6px;background:#2c3e50;color:white;font-size:10px}}
.tbl td{{padding:4px 6px;border-bottom:1px solid #eee}}
.tbl tr.act{{background:#e8f8f5}}
.tbl tr.skip{{opacity:0.4}}
.tbl tr:hover{{opacity:1;background:#f0f7ff}}
.mono{{font-family:monospace;font-size:11px}}
.meta{{font-size:12px;color:#777;margin:5px 0 20px}}
.summary{{display:flex;gap:12px;flex-wrap:wrap;margin:15px 0}}
.stat{{background:white;border:1px solid #ddd;border-radius:8px;padding:10px;min-width:100px;text-align:center}}
.stat-val{{font-size:18px;font-weight:bold;color:#2980b9}}
.stat-lbl{{font-size:10px;color:#777}}
.guide{{font-size:12px;color:#555;background:#fff;border:1px solid #ddd;border-radius:8px;padding:12px;margin:15px 0}}
.guide b{{color:#333}}
.footer{{margin-top:30px;padding-top:10px;border-top:1px solid #ddd;font-size:11px;color:#999}}
</style></head><body>

<h1>Pharmit zone mapping: {campaign_id}</h1>
<p class="meta">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} | Module 06b v3.1 (ligand-atom zones + evidence)</p>

<div class="summary">
  <div class="stat"><div class="stat-val">{n_active}</div><div class="stat-lbl">Drug-like features</div></div>
  <div class="stat"><div class="stat-val" style="color:#1e8449">{n_high}</div><div class="stat-lbl">HIGH confidence</div></div>
  <div class="stat"><div class="stat-val" style="color:#b7950b">{n_med}</div><div class="stat-lbl">MEDIUM</div></div>
  <div class="stat"><div class="stat-val" style="color:#c0392b">{n_low}</div><div class="stat-lbl">LOW</div></div>
  <div class="stat"><div class="stat-val" style="color:#5b6abf">{n_plip}</div><div class="stat-lbl">PLIP only</div></div>
  <div class="stat"><div class="stat-val" style="color:#999">{n_none}</div><div class="stat-lbl">No evidence</div></div>
</div>

<div class="guide">
  <b>How to use:</b> Upload a variant JSON to Pharmit. If too few hits, start removing features
  with LOW/NONE confidence first — they have weakest receptor evidence. Keep HIGH features
  (strong footprint energy) as anchors.<br>
  <b>Confidence:</b> HIGH = footprint &lt; -5 kcal/mol | MEDIUM = &lt; -1 | LOW = &lt; 0 | PLIP = geometric only | NONE = no receptor contact
</div>

<h2>Feature &rarr; zone mapping with evidence</h2>

<table class="tbl">
<tr>
  <th>#</th><th>Pharmit type</th><th>Orig</th>
  <th>Zone</th><th>Ligand atom</th><th>Coordinates</th>
  <th>Receptor</th><th>PLIP type</th>
  <th>Energy</th><th>vdW</th><th>ES</th>
  <th>Confidence</th><th>Decision</th>
</tr>
{rows_html}
</table>

<h2>Generated variants</h2>
<p style="font-size:12px;color:#777">Upload each JSON to Pharmit separately. Compare hit counts and scaffolds.</p>

<table class="tbl">
<tr><th>Strategy</th><th>Active zones</th><th>Description</th><th>File</th></tr>
{strat_rows}
</table>

<div class="footer">reference_docking | Module 06b | Pharmit Zone Selector v3.1 (evidence) | {datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
</body></html>"""

    path.write_text(html, encoding="utf-8")
    logger.info(f"  Saved: {path}")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_pharmit_zone_selector(
        pharmit_json_path: str,
        output_dir: str,
        campaign_id: str = "",
        plip_json_path: Optional[str] = None,
        footprint_csv_path: Optional[str] = None,
        strategies: Optional[Dict] = None,
        max_distance: float = 2.0,
        **kwargs,
) -> Dict[str, Any]:
    """
    Map Pharmit features to zones and annotate with evidence.

    v3.1: Zones by ligand atom, evidence from PLIP + footprint.
    """
    logger.info("=" * 60)
    logger.info("  06b PHARMIT ZONE SELECTOR v3.1")
    logger.info("  Zones: ligand atom | Evidence: PLIP + footprint")
    logger.info("=" * 60)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Load Pharmit JSON ---
    logger.info(f"  Pharmit JSON: {Path(pharmit_json_path).name}")
    with open(pharmit_json_path) as f:
        pharmit_data = json.load(f)

    points = pharmit_data.get("points", [])
    n_enabled = sum(1 for p in points if p.get("enabled"))
    logger.info(f"  Features: {len(points)} ({n_enabled} enabled)")

    # --- Parse ligand atoms ---
    logger.info("")
    logger.info("  Parsing ligand:")
    ligand_atoms = _parse_ligand_from_pharmit(pharmit_data)
    if not ligand_atoms:
        return {"success": False, "error": "No ligand atoms in Pharmit JSON"}

    # --- Load evidence ---
    logger.info("")
    logger.info("  Loading evidence:")
    plip_anchors = _load_plip_anchors(plip_json_path)
    footprint_energy = _load_footprint_energy(footprint_csv_path)

    # --- Map + annotate ---
    logger.info("")
    logger.info("  Mapping features:")
    mappings = map_features_to_zones(
        points, ligand_atoms, plip_anchors, footprint_energy,
        max_distance=max_distance,
    )

    # --- Save CSV ---
    mapping_csv = out_dir / "pharmit_zone_mapping.csv"
    write_mapping_csv(mappings, mapping_csv)

    # --- Generate variants ---
    logger.info("")
    logger.info("  Generating zone variants:")
    if strategies is None:
        strategies = ZONE_STRATEGIES
    outputs = generate_variants(pharmit_data, mappings, out_dir, strategies)

    # --- Save HTML ---
    mapping_html = out_dir / "pharmit_zone_mapping.html"
    write_mapping_html(mappings, strategies, outputs, campaign_id, mapping_html)

    # --- Summary ---
    n_active = sum(1 for m in mappings if m["decision"] == "ZONE_ACTIVE")
    n_disabled = sum(1 for m in mappings if m["decision"] in ("DISABLE", "KEEP_OFF"))

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  DONE — {len(mappings)} features, {len(outputs)} variants")
    logger.info(f"  Drug-like: {n_active} | Disabled: {n_disabled}")
    logger.info("=" * 60)

    for sname, spath in outputs.items():
        logger.info(f"  Upload to Pharmit: {spath}")

    return {
        "success": True,
        "n_features": len(mappings),
        "n_active": n_active,
        "n_variants": len(outputs),
        "mapping_csv": str(mapping_csv),
        "mapping_html": str(mapping_html),
        "variants": outputs,
    }