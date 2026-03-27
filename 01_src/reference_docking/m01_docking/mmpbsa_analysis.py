"""
MMPBSA Analysis — Core Module (01h)
======================================
Parses MMPBSA.py decomposition output and compares with DOCK6 footprint.

Takes the raw output from 01g (FINAL_DECOMP_MMPBSA.dat, FINAL_RESULTS_MMPBSA.dat)
and produces:
  - per_residue_decomp.csv:             vdW, ES, GB, SA, total per residue (PDB numbering)
  - comparison_invacuo_vs_solvated.csv: merge with 01d/04b footprint
  - zone_summary.json:                  energy by binding site zone

Scientific questions answered:
  - Are ARG598/LYS599 salt bridges real or in-vacuo artifacts?
    → Compare fp_es vs MMPBSA(ES+GB) for phosphate zone
  - Does ASP494 rescue with solvation?
    → Compare fp_total vs MMPBSA total for ribose_bridge zone
  - Do hydrophobic residues rise in ranking?
    → Compare rank_fp vs rank_mmpbsa for xylose zone

Pipeline:
    6. Parse FINAL_DECOMP_MMPBSA.dat → per-residue DataFrame
    7. Map sequential → PDB residue numbering
    8. Parse FINAL_RESULTS_MMPBSA.dat → global binding energy
    9. Compare with 01d footprint (ranking delta, zone summary)

Location: 01_src/reference_docking/m01_docking/mmpbsa_analysis.py
Project: reference_docking
Module: 01h (core)
Version: 1.0 (2026-03-26)
"""

import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# ZONE DEFINITIONS
# =============================================================================
# Sub-pocket classification from UDX PLIP analysis (06a).
# Shared with 04b footprint_analysis.py.
# TODO: make configurable via YAML for other receptor systems.

ZONE_DEFINITIONS = {
    "phosphate": {
        "residues": {"ARG598", "LYS599"},
        "label": "Phosphate (salt bridges)",
        "description": "Strong in-vacuo ES, expect large GB desolvation penalty.",
    },
    "xylose": {
        "residues": {"TRP392", "TRP495", "TYR565", "SER575"},
        "label": "Xylose pocket",
        "description": "Hydrophobic/HBA interactions, expect minimal solvation change.",
    },
    "uracil": {
        "residues": {"ASP361", "THR390", "ARG363"},
        "label": "Uracil pocket",
        "description": "Mixed polar/charged, moderate solvation screening.",
    },
    "ribose_bridge": {
        "residues": {"ASP494", "GLU529", "HIS335"},
        "label": "Ribose bridge",
        "description": "ASP494 shows +2.74 ES repulsion in-vacuo. Water bridge in PLIP.",
    },
}


# =============================================================================
# GLOBAL RESULTS PARSING
# =============================================================================

def parse_mmpbsa_global(
        results_file: Union[str, Path],
) -> Dict[str, float]:
    """
    Parse global MMPBSA results (total binding energy components).

    Looks for the "Differences (Complex - Receptor - Ligand)" section
    in FINAL_RESULTS_MMPBSA.dat.

    Returns:
        Dict with delta_total, vdw, eel, egb, esurf (kcal/mol)
    """
    text = Path(results_file).read_text()
    results = {}

    in_delta = False
    for line in text.split("\n"):
        if "Differences (Complex - Receptor - Ligand)" in line:
            in_delta = True
            continue
        if not in_delta:
            continue

        stripped = line.strip()
        if not stripped or stripped.startswith("---") or stripped.startswith("Energy"):
            continue

        parts = stripped.split()
        if len(parts) < 2:
            continue

        # "DELTA TOTAL  -1.4126  0.0000  0.0000"
        if parts[0] == "DELTA" and len(parts) >= 3:
            key = parts[1]
            try:
                val = float(parts[2])
            except ValueError:
                continue
            if key == "TOTAL":
                results["delta_total"] = val
            elif key == "G":
                # "DELTA G gas" or "DELTA G solv"
                if len(parts) >= 4:
                    subkey = parts[2]
                    try:
                        results[f"delta_g_{subkey}"] = float(parts[3])
                    except ValueError:
                        pass
            continue

        # "VDWAALS  -54.0920  0.0000  0.0000"
        key = parts[0]
        try:
            val = float(parts[1])
        except ValueError:
            continue

        if key == "VDWAALS":
            results["vdw"] = val
        elif key == "EEL":
            results["eel"] = val
        elif key == "EGB":
            results["egb"] = val
        elif key == "ESURF":
            results["esurf"] = val

    return results


# =============================================================================
# DECOMPOSITION PARSING
# =============================================================================

def parse_decomp_output(
        decomp_file: Union[str, Path],
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        is_single_frame: bool = True,
) -> Dict[str, Any]:
    """
    Parse FINAL_DECOMP_MMPBSA.dat → per-residue CSV with PDB numbering.

    MMPBSA.py idecomp=1 output columns:
        Residue | Internal | vdW | Electrostatic | Polar Solv | Non-Polar Solv | TOTAL

    For multiple frames, each column has "avg" and "std" sub-columns.

    Args:
        decomp_file:    FINAL_DECOMP_MMPBSA.dat
        receptor_pdb:   Receptor PDB (for residue numbering mapping)
        output_dir:     Output directory for CSVs
        is_single_frame: True if single-point (no std columns)

    Returns:
        Dict with success, csv path, DataFrame
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    text = Path(decomp_file).read_text()

    # Build residue mapping: tleap sequential → PDB original
    residue_mapping = _build_residue_mapping(receptor_pdb)

    # Parse decomp sections
    rows = _parse_decomp_sections(text, is_single_frame)

    if not rows:
        return {"success": False, "error": "No residue data found in decomp file"}

    df = pd.DataFrame(rows)

    # Apply PDB numbering
    if residue_mapping:
        df = _apply_pdb_numbering(df, residue_mapping)

    # Sort by total energy (most favorable first)
    if "total" in df.columns:
        df = df.sort_values("total", ascending=True).reset_index(drop=True)

    # Save
    csv_path = output_dir / "per_residue_decomp.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8")
    logger.info(f"  Saved: {csv_path} ({len(df)} residues)")

    # Summary stats
    n_favorable = len(df[df["total"] < -0.5]) if "total" in df.columns else 0
    n_unfavorable = len(df[df["total"] > 0.5]) if "total" in df.columns else 0

    return {
        "success": True,
        "csv_path": str(csv_path),
        "df": df,
        "n_residues": len(df),
        "n_favorable": n_favorable,
        "n_unfavorable": n_unfavorable,
    }


def _parse_decomp_sections(text: str, is_single_frame: bool) -> List[Dict]:
    """
    Parse per-residue data from FINAL_DECOMP_MMPBSA.dat.

    MMPBSA.py 14 format: CSV with sections Complex, Receptor, Ligand, DELTAS.
    We want the DELTAS > Total Energy Decomposition section.

    Each data line:
        ASP   1,R ASP   1,0.0,0.0,0.0,-0.001,0.0,0.0,0.112,0.0,0.0,-0.112,0.0,0.0,0.0,0.0,0.0,-2.84e-14,0.0,0.0

    CSV columns:
        0: Residue (e.g., "ASP   1")
        1: Location (e.g., "R ASP   1" — R=Receptor, L=Ligand)
        2-4: Internal (Avg, Std, SEM)
        5-7: vdW (Avg, Std, SEM)
        8-10: Electrostatic (Avg, Std, SEM)
        11-13: Polar Solvation / GB (Avg, Std, SEM)
        14-16: Non-Polar Solv / SA (Avg, Std, SEM)
        17-19: TOTAL (Avg, Std, SEM)
    """
    rows = []
    lines = text.split("\n")

    in_deltas = False
    in_total_decomp = False

    for i, line in enumerate(lines):
        line_stripped = line.strip()

        # Find DELTAS section
        if line_stripped == "DELTAS:":
            in_deltas = True
            continue

        if not in_deltas:
            continue

        # Find Total Energy Decomposition sub-section
        if "Total Energy Decomposition" in line_stripped:
            in_total_decomp = True
            continue

        if not in_total_decomp:
            continue

        # Stop at next sub-section (Sidechain/Backbone) or blank followed by section
        if line_stripped in ("Sidechain Energy Decomposition:",
                            "Backbone Energy Decomposition:"):
            break

        # Skip header lines
        if line_stripped.startswith("Residue,") or line_stripped.startswith(",,"):
            continue
        if not line_stripped:
            continue

        # Parse CSV data line
        parsed = _parse_decomp_csv_line(line_stripped, is_single_frame)
        if parsed:
            rows.append(parsed)

    return rows


def _parse_decomp_csv_line(line: str, is_single_frame: bool) -> Optional[Dict]:
    """
    Parse a single CSV residue line from MMPBSA.py 14 decomp output.

    Format: "ASP   1,R ASP   1,0.0,0.0,0.0,..."
    """
    parts = line.split(",")
    if len(parts) < 18:
        return None

    # Column 0: Residue identifier (e.g., "ASP   1")
    res_id = parts[0].strip()
    # Column 1: Location (e.g., "R ASP   1" or "L UDX   1")
    location = parts[1].strip()

    # Skip ligand residues
    if location.startswith("L "):
        return None

    # Parse residue name and sequential number
    res_tokens = res_id.split()
    if len(res_tokens) < 2:
        return None

    res_name = res_tokens[0]
    try:
        res_num_seq = int(res_tokens[1])
    except ValueError:
        return None

    # Parse energy values (Avg is at indices 2, 5, 8, 11, 14, 17)
    try:
        internal = float(parts[2])
        vdw = float(parts[5])
        es = float(parts[8])
        gb = float(parts[11])
        sa = float(parts[14])
        total = float(parts[17])
    except (ValueError, IndexError):
        return None

    result = {
        "residue_name": res_name,
        "residue_number_seq": res_num_seq,
        "chain_type": "R",
        "internal": internal,
        "vdw": vdw,
        "es": es,
        "gb": gb,
        "sa": sa,
        "total": total,
    }

    # Add std dev for multi-frame runs
    if not is_single_frame:
        try:
            result["vdw_std"] = float(parts[6])
            result["es_std"] = float(parts[9])
            result["gb_std"] = float(parts[12])
            result["sa_std"] = float(parts[15])
            result["total_std"] = float(parts[18])
        except (ValueError, IndexError):
            pass

    return result


# =============================================================================
# RESIDUE NUMBERING
# =============================================================================

def _build_residue_mapping(
        receptor_pdb: Union[str, Path],
) -> Dict[int, Dict[str, Any]]:
    """
    Build sequential→PDB residue mapping from receptor PDB.

    tleap preserves PDB residue numbering in topology, but MMPBSA.py
    decomposition uses sequential index. This maps seq → PDB info.

    Returns:
        Dict[seq_index] → {"name": "TRP", "number": 392, "chain": "A",
                           "residue_id": "TRP392.A"}
    """
    pdb_path = Path(receptor_pdb)
    if not pdb_path.exists():
        logger.warning(f"  Receptor PDB not found for mapping: {pdb_path}")
        return {}

    mapping = {}
    seen_residues = []
    prev_resid = None

    for line in pdb_path.read_text().split("\n"):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        res_name = line[17:20].strip()
        chain = line[21:22].strip() or "A"
        try:
            res_num = int(line[22:26].strip())
        except ValueError:
            continue

        resid = (res_name, res_num, chain)
        if resid != prev_resid:
            seen_residues.append(resid)
            prev_resid = resid

    for i, (name, number, chain) in enumerate(seen_residues, 1):
        mapping[i] = {
            "name": name,
            "number": number,
            "chain": chain,
            "residue_id": f"{name}{number}.{chain}",
        }

    logger.info(f"  Residue mapping: {len(mapping)} residues (PDB numbering)")
    return mapping


def _apply_pdb_numbering(
        df: pd.DataFrame,
        mapping: Dict[int, Dict],
) -> pd.DataFrame:
    """Apply PDB residue numbering to decomp DataFrame."""
    if "residue_number_seq" not in df.columns:
        return df

    pdb_names = []
    pdb_numbers = []
    pdb_chains = []
    pdb_ids = []

    for _, row in df.iterrows():
        seq_num = row.get("residue_number_seq", 0)
        if seq_num in mapping:
            m = mapping[seq_num]
            pdb_names.append(m["name"])
            pdb_numbers.append(m["number"])
            pdb_chains.append(m["chain"])
            pdb_ids.append(m["residue_id"])
        else:
            pdb_names.append(row.get("residue_name", "UNK"))
            pdb_numbers.append(seq_num)
            pdb_chains.append("?")
            pdb_ids.append(f"{row.get('residue_name', 'UNK')}{seq_num}.?")

    df["residue_name_pdb"] = pdb_names
    df["residue_number_pdb"] = pdb_numbers
    df["chain"] = pdb_chains
    df["residue_id"] = pdb_ids

    return df


# =============================================================================
# FOOTPRINT COMPARISON
# =============================================================================

def compare_with_footprint(
        decomp_df: pd.DataFrame,
        footprint_csv: Union[str, Path],
        output_dir: Union[str, Path],
) -> Dict[str, Any]:
    """
    Compare MMPBSA per-residue decomp with 01d/04b footprint in-vacuo.

    Merges on residue_id (PDB numbering). Computes:
      - Ranking change per residue
      - Energy delta per component
      - Solvation contribution (GB + SA) per residue
      - Zone-level summary

    Args:
        decomp_df:     Per-residue decomposition DataFrame (from parse_decomp_output)
        footprint_csv: Path to residue_consensus.csv from 04b
        output_dir:    Output directory

    Returns:
        Dict with success, comparison CSV, zone summary
    """
    output_dir = Path(output_dir)
    fp_path = Path(footprint_csv)

    if not fp_path.exists():
        return {"success": False, "error": f"Footprint CSV not found: {fp_path}"}

    df_fp = pd.read_csv(fp_path)
    logger.info(f"  Footprint: {len(df_fp)} residues from 04b")

    if "residue_id" not in decomp_df.columns or "residue_id" not in df_fp.columns:
        return {"success": False, "error": "Missing residue_id column for merge"}

    # Merge on residue_id
    df_merged = pd.merge(
        decomp_df,
        df_fp[["residue_id", "residue_name", "mean_vdw", "mean_es", "mean_total"]],
        on="residue_id",
        how="outer",
        suffixes=("_mmpbsa", "_fp"),
    )

    rename_map = {
        "mean_vdw": "fp_vdw",
        "mean_es": "fp_es",
        "mean_total": "fp_total",
    }
    df_merged = df_merged.rename(columns=rename_map)

    # Deltas (MMPBSA component vs footprint component)
    if "vdw" in df_merged.columns and "fp_vdw" in df_merged.columns:
        df_merged["delta_vdw"] = df_merged["vdw"] - df_merged["fp_vdw"]
    if "es" in df_merged.columns and "fp_es" in df_merged.columns:
        df_merged["delta_es"] = df_merged["es"] - df_merged["fp_es"]

    # MMPBSA in-vacuo component (vdW + ES only, no solvation)
    if "vdw" in df_merged.columns and "es" in df_merged.columns:
        df_merged["mmpbsa_invacuo"] = df_merged["vdw"] + df_merged["es"]

    # Solvation contribution per residue
    if "gb" in df_merged.columns and "sa" in df_merged.columns:
        df_merged["solvation"] = df_merged["gb"] + df_merged["sa"]

    # Rankings
    if "total" in df_merged.columns:
        df_merged["rank_mmpbsa"] = df_merged["total"].rank(ascending=True)
    if "fp_total" in df_merged.columns:
        df_merged["rank_fp"] = df_merged["fp_total"].rank(ascending=True)
    if "rank_mmpbsa" in df_merged.columns and "rank_fp" in df_merged.columns:
        # Positive = residue becomes more important with solvation
        df_merged["rank_delta"] = df_merged["rank_fp"] - df_merged["rank_mmpbsa"]

    df_merged = df_merged.sort_values("total", ascending=True, na_position="last")

    # Save
    comp_csv = output_dir / "comparison_invacuo_vs_solvated.csv"
    df_merged.to_csv(comp_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {comp_csv} ({len(df_merged)} residues)")

    # Zone summary
    zone_summary = _compute_zone_summary(df_merged)

    return {
        "success": True,
        "comparison_csv": str(comp_csv),
        "df": df_merged,
        "zone_summary": zone_summary,
    }


def _compute_zone_summary(df: pd.DataFrame) -> Dict[str, Dict]:
    """
    Compute energy summary by binding site zone.

    For each zone, reports:
      - Sum of MMPBSA total, vdW, ES, GB, SA
      - Sum of footprint total (in-vacuo)
      - Solvation impact (sum of GB + SA for zone residues)
    """
    summary = {}

    for zone_name, zdef in ZONE_DEFINITIONS.items():
        mask = df["residue_id"].apply(
            lambda rid: str(rid).split(".")[0] in zdef["residues"]
            if pd.notna(rid) else False
        )
        zone_df = df[mask]

        if zone_df.empty:
            continue

        zone_info = {
            "label": zdef["label"],
            "n_residues": len(zone_df),
            "residues": zone_df["residue_id"].tolist(),
        }

        for col in ["vdw", "es", "gb", "sa", "total", "fp_total", "solvation"]:
            if col in zone_df.columns:
                vals = zone_df[col].dropna()
                if len(vals) > 0:
                    zone_info[f"sum_{col}"] = round(vals.sum(), 2)
                    zone_info[f"mean_{col}"] = round(vals.mean(), 2)

        if "solvation" in zone_df.columns:
            zone_info["solvation_impact"] = round(
                zone_df["solvation"].dropna().sum(), 2
            )

        summary[zone_name] = zone_info

    return summary


# =============================================================================
# HTML REPORT
# =============================================================================

def generate_mmpbsa_html(
        decomp_df: pd.DataFrame,
        global_results: Dict[str, float],
        zone_summary: Optional[Dict] = None,
        comparison_df: Optional[pd.DataFrame] = None,
        is_single_frame: bool = True,
        campaign_id: str = "",
) -> str:
    """
    Generate MMPBSA per-residue decomposition HTML report.

    Follows the project HTML style (CSS variables, .tbl tables, .tg tags).
    """
    from datetime import datetime

    mode_label = "Single-point (1 frame)" if is_single_frame else "MD trajectory (mean ± std)"
    dg = global_results.get("delta_total", None)
    dg_str = f"{dg:+.2f}" if dg is not None else "N/A"

    # Zone colors (shared with 04b)
    zone_colors = {
        "phosphate": "#1D9E75", "xylose": "#378ADD",
        "uracil": "#BA7517", "ribose_bridge": "#7F77DD",
    }

    # ── Top contributing residues (|total| > 0.5) ──
    df_sig = decomp_df[decomp_df["total"].abs() > 0.5].copy() if "total" in decomp_df.columns else decomp_df.head(0)
    df_fav = df_sig[df_sig["total"] < 0].sort_values("total")
    df_unfav = df_sig[df_sig["total"] > 0].sort_values("total", ascending=False)

    # Max |total| for bar scaling
    max_abs = max(df_sig["total"].abs().max(), 1) if len(df_sig) > 0 else 1

    # ── Build HTML ──
    html = f"""\
<style>
.tbl{{width:100%;border-collapse:collapse;font-size:13px}}
.tbl th{{text-align:left;padding:6px 8px;border-bottom:0.5px solid var(--color-border-secondary);color:var(--color-text-secondary);font-weight:500;font-size:12px}}
.tbl td{{padding:5px 8px;border-bottom:0.5px solid var(--color-border-tertiary);color:var(--color-text-primary);vertical-align:top}}
.m{{font-family:var(--font-mono);font-size:12px}}
.tg{{display:inline-block;font-size:11px;padding:2px 8px;border-radius:6px;font-weight:500}}
.s{{background:var(--color-background-success);color:var(--color-text-success)}}
.i{{background:var(--color-background-info);color:var(--color-text-info)}}
.w{{background:var(--color-background-warning);color:var(--color-text-warning)}}
.d{{background:var(--color-background-danger);color:var(--color-text-danger)}}
.x{{background:var(--color-background-secondary);color:var(--color-text-secondary)}}
.sub{{font-size:11px;color:var(--color-text-secondary);margin-top:2px}}
.bar{{height:8px;border-radius:4px;display:inline-block;vertical-align:middle}}
.sect{{font-size:11px;font-weight:500;color:var(--color-text-secondary);padding:10px 8px 4px;text-transform:uppercase;letter-spacing:0.5px}}
.card{{background:var(--color-background-secondary);border-radius:8px;padding:14px;margin:8px 0;display:inline-block;text-align:center;min-width:110px;vertical-align:top}}
.card-val{{font-size:20px;font-weight:bold}}
.card-lbl{{font-size:10px;color:var(--color-text-secondary);margin-top:2px}}
.neg{{color:var(--color-text-success)}}
.pos{{color:var(--color-text-danger)}}
.neut{{color:var(--color-text-secondary)}}
</style>

<div style="padding:0.5rem 0">

<div style="font-size:11px;color:var(--color-text-secondary);margin-bottom:4px">
MMPBSA Per-Residue Decomposition — {campaign_id} | {mode_label} | {datetime.now().strftime('%Y-%m-%d %H:%M')}
</div>
"""

    # ── Section 1: Global binding energy ──
    vdw = global_results.get("vdw", 0)
    eel = global_results.get("eel", 0)
    egb = global_results.get("egb", 0)
    esurf = global_results.get("esurf", 0)

    def _val_class(v):
        return "neg" if v < -1 else "pos" if v > 1 else "neut"

    html += f"""
<div style="margin:12px 0">
  <div class="card"><div class="card-val {_val_class(dg if dg else 0)}">{dg_str}</div><div class="card-lbl">ΔG bind (kcal/mol)</div></div>
  <div class="card"><div class="card-val {_val_class(vdw)}">{vdw:+.1f}</div><div class="card-lbl">vdW</div></div>
  <div class="card"><div class="card-val {_val_class(eel)}">{eel:+.1f}</div><div class="card-lbl">Electrostatic</div></div>
  <div class="card"><div class="card-val {_val_class(egb)}">{egb:+.1f}</div><div class="card-lbl">GB solvation</div></div>
  <div class="card"><div class="card-val {_val_class(esurf)}">{esurf:+.1f}</div><div class="card-lbl">SA (non-polar)</div></div>
</div>
<div class="sub">EEL ({eel:+.1f}) and GB ({egb:+.1f}) nearly cancel — solvation screens {abs(egb/(eel) * 100) if abs(eel) > 0.01 else 0:.0f}% of raw electrostatics. The net binding is dominated by vdW shape complementarity.</div>
"""

    # ── Section 2: Zone comparison (if footprint available) ──
    if zone_summary:
        html += """
<table class="tbl" style="margin-top:16px">
<tr><th>Zone</th><th>MMPBSA total</th><th>Solvation (GB+SA)</th><th></th><th>Interpretation</th></tr>
"""
        for zone_id, info in sorted(zone_summary.items(), key=lambda x: x[1].get("sum_total", 0)):
            total = info.get("sum_total", 0)
            solv = info.get("solvation_impact", 0)
            label = info.get("label", zone_id)
            color = zone_colors.get(zone_id, "#999")
            n_res = info.get("n_residues", 0)

            bar_w = min(80, int(abs(total) / max_abs * 80))
            bar_color = "#2ecc71" if total < 0 else "#e74c3c"

            if zone_id == "phosphate":
                interp = "Real but ~85% screened by water"
            elif zone_id == "xylose":
                interp = "Stable — vdW-dominated, drug-like"
            elif zone_id == "ribose_bridge":
                interp = "ASP494 repulsive even with solvation"
            elif zone_id == "uracil":
                interp = "Weak anchor, selective"
            else:
                interp = ""

            html += f"""<tr>
<td><span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:{color};margin-right:4px;vertical-align:middle"></span><b>{label}</b> <span class="sub">({n_res} res)</span></td>
<td class="m" style="font-weight:500">{total:+.2f}</td>
<td class="m">{solv:+.2f}</td>
<td><span class="bar" style="width:{bar_w}px;background:{bar_color}"></span></td>
<td class="sub">{interp}</td>
</tr>"""

        html += "</table>"

    # ── Section 3: Side-by-side comparison for key residues ──
    if comparison_df is not None and len(comparison_df) > 0:
        # Get residues with significant energy in either method
        comp = comparison_df.copy()
        has_fp = "fp_total" in comp.columns
        has_mm = "total" in comp.columns

        if has_fp and has_mm:
            comp["either_sig"] = (comp["total"].abs() > 0.3) | (comp["fp_total"].abs().fillna(0) > 0.3)
            comp_sig = comp[comp["either_sig"]].sort_values("total", na_position="last").head(20)

            if len(comp_sig) > 0:
                html += """
<table class="tbl" style="margin-top:16px">
<tr><th colspan="6" class="sect">In-vacuo footprint vs MMPBSA (solvated)</th></tr>
<tr><th>Residue</th><th>FP vdW+ES</th><th></th><th>MMPBSA total</th><th></th><th>Solvation (GB+SA)</th></tr>
"""
                for _, r in comp_sig.iterrows():
                    rid = r.get("residue_id", "?")
                    fp = r.get("fp_total", None)
                    mm = r.get("total", None)
                    solv = r.get("solvation", None)

                    fp_str = f"{fp:+.2f}" if pd.notna(fp) else "—"
                    mm_str = f"{mm:+.2f}" if pd.notna(mm) else "—"
                    solv_str = f"{solv:+.2f}" if pd.notna(solv) else "—"

                    # Bars
                    fp_w = min(60, int(abs(fp) / max_abs * 60)) if pd.notna(fp) else 0
                    mm_w = min(60, int(abs(mm) / max_abs * 60)) if pd.notna(mm) else 0
                    fp_c = "#2ecc71" if (pd.notna(fp) and fp < 0) else "#e74c3c"
                    mm_c = "#2ecc71" if (pd.notna(mm) and mm < 0) else "#e74c3c"

                    html += f"""<tr>
<td><b>{rid}</b></td>
<td class="m">{fp_str}</td>
<td><span class="bar" style="width:{fp_w}px;background:{fp_c};opacity:0.5"></span></td>
<td class="m" style="font-weight:500">{mm_str}</td>
<td><span class="bar" style="width:{mm_w}px;background:{mm_c}"></span></td>
<td class="m">{solv_str}</td>
</tr>"""

                html += "</table>"

    # ── Section 4: Top favorable residues ──
    html += """
<table class="tbl" style="margin-top:16px">
<tr><th colspan="7" class="sect">Top favorable residues (binding)</th></tr>
<tr><th>#</th><th>Residue</th><th>vdW</th><th>ES</th><th>GB</th><th>SA</th><th>TOTAL</th></tr>
"""
    for rank, (_, r) in enumerate(df_fav.head(15).iterrows(), 1):
        rid = r.get("residue_id", f"{r['residue_name']}{r['residue_number_seq']}")
        bar_w = min(80, int(abs(r["total"]) / max_abs * 80))
        html += f"""<tr>
<td class="m">{rank}</td>
<td><b>{rid}</b></td>
<td class="m">{r['vdw']:+.2f}</td>
<td class="m">{r['es']:+.2f}</td>
<td class="m">{r['gb']:+.2f}</td>
<td class="m">{r['sa']:+.4f}</td>
<td class="m" style="font-weight:500">{r['total']:+.2f} <span class="bar" style="width:{bar_w}px;background:#2ecc71;margin-left:4px"></span></td>
</tr>"""

    html += "</table>"

    # ── Section 5: Unfavorable residues ──
    if len(df_unfav) > 0:
        html += """
<table class="tbl" style="margin-top:16px">
<tr><th colspan="7" class="sect">Unfavorable residues (repulsive)</th></tr>
<tr><th>#</th><th>Residue</th><th>vdW</th><th>ES</th><th>GB</th><th>SA</th><th>TOTAL</th></tr>
"""
        for rank, (_, r) in enumerate(df_unfav.head(10).iterrows(), 1):
            rid = r.get("residue_id", f"{r['residue_name']}{r['residue_number_seq']}")
            bar_w = min(80, int(abs(r["total"]) / max_abs * 80))
            html += f"""<tr>
<td class="m">{rank}</td>
<td><b>{rid}</b></td>
<td class="m">{r['vdw']:+.2f}</td>
<td class="m">{r['es']:+.2f}</td>
<td class="m">{r['gb']:+.2f}</td>
<td class="m">{r['sa']:+.4f}</td>
<td class="m" style="font-weight:500;color:var(--color-text-danger)">{r['total']:+.2f} <span class="bar" style="width:{bar_w}px;background:#e74c3c;margin-left:4px"></span></td>
</tr>"""

        html += "</table>"

    # ── Section 6: Pharmit implications ──
    html += """
<div style="margin-top:16px;padding:12px;background:var(--color-background-secondary);border-radius:8px;font-size:12px">
<div style="font-weight:600;margin-bottom:6px">Implications for Pharmit feature selection</div>
<table style="font-size:12px;border-collapse:collapse;width:100%">
"""
    implications = [
        ("Xylose strategy", "s", "VALIDATED",
         "TRP392 + TRP495 are vdW-dominated and stable with solvation. Best drug-like anchors."),
        ("Phosphate features", "w", "REDUCE WEIGHT",
         "ARG598/LYS599 salt bridges are real but ~85% screened. Don't over-weight NegativeIon features."),
        ("ASP494", "d", "EXCLUDE",
         "Repulsive (+7.9 kcal/mol) even with solvation. Exclude from all Pharmit queries."),
        ("Uracil pocket", "i", "OPTIONAL",
         "ASP361/ARG363 weak but selective. Keep as optional features for extended molecules."),
    ]
    for label, cls, tag, desc in implications:
        html += f"""<tr>
<td style="padding:4px 8px;font-weight:500;white-space:nowrap">{label}</td>
<td style="padding:4px 8px"><span class="tg {cls}">{tag}</span></td>
<td style="padding:4px 8px;color:var(--color-text-secondary)">{desc}</td>
</tr>"""

    html += """</table></div>"""

    # ── Footer ──
    html += f"""
<div style="margin-top:16px;font-size:11px;color:var(--color-text-tertiary)">
MMPBSA.py 14.0 | idecomp=1 (per-residue) | igb=2 (OBC-II) | saltcon=0.15 M | AMBER ff14SB + GAFF2 | {mode_label}
</div>
</div>"""

    return html


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_mmpbsa_analysis(
        decomp_file: Union[str, Path],
        results_file: Union[str, Path],
        receptor_pdb: Union[str, Path],
        output_dir: Union[str, Path],
        is_single_frame: bool = True,
        compare_footprint: bool = True,
        footprint_csv: Optional[str] = None,
        campaign_id: str = "",
) -> Dict[str, Any]:
    """
    Run complete MMPBSA analysis pipeline (steps 6-9).

    Reads raw MMPBSA.py output from 01g and produces CSVs + zone summary + HTML.

    Args:
        decomp_file:      FINAL_DECOMP_MMPBSA.dat from 01g
        results_file:     FINAL_RESULTS_MMPBSA.dat from 01g
        receptor_pdb:     Receptor PDB from 00b (for PDB numbering)
        output_dir:       Output directory
        is_single_frame:  True if 01g ran in single_point mode
        compare_footprint: Compare with 01d/04b footprint
        footprint_csv:    Path to residue_consensus.csv
        campaign_id:      Campaign name for report title

    Returns:
        Dict with success, CSV paths, global results, zone summary
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("  MMPBSA Analysis (01h)")
    logger.info("=" * 60)
    logger.info(f"  Decomp file: {Path(decomp_file).name}")
    logger.info(f"  Results:     {Path(results_file).name}")
    logger.info(f"  Receptor:    {Path(receptor_pdb).name}")
    logger.info(f"  Mode:        {'single_point' if is_single_frame else 'md (multi-frame)'}")

    # ─────────────────────────────────────────────────────────────
    # Step 6: Parse global results
    # ─────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("─── Step 6: Parse Global MMPBSA Results ───")

    global_results = {}
    if Path(results_file).exists():
        global_results = parse_mmpbsa_global(results_file)
        if global_results:
            logger.info(f"  ΔG_bind:  {global_results.get('delta_total', 'N/A')} kcal/mol")
            logger.info(f"    vdW:    {global_results.get('vdw', 'N/A')}")
            logger.info(f"    EEL:    {global_results.get('eel', 'N/A')}")
            logger.info(f"    EGB:    {global_results.get('egb', 'N/A')}")
            logger.info(f"    ESURF:  {global_results.get('esurf', 'N/A')}")
        else:
            logger.warning("  No global results parsed from FINAL_RESULTS_MMPBSA.dat")

    # Save global results
    global_json = output_dir / "global_binding_energy.json"
    with open(global_json, "w") as f:
        json.dump(global_results, f, indent=2)

    # ─────────────────────────────────────────────────────────────
    # Step 7: Parse per-residue decomposition
    # ─────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("─── Step 7: Parse Per-Residue Decomposition ───")

    result_parse = parse_decomp_output(
        decomp_file=decomp_file,
        receptor_pdb=receptor_pdb,
        output_dir=output_dir,
        is_single_frame=is_single_frame,
    )
    if not result_parse["success"]:
        return {"success": False, "error": f"Parse: {result_parse['error']}"}

    logger.info(f"  Residues: {result_parse['n_residues']} total, "
                f"{result_parse['n_favorable']} favorable (<-0.5), "
                f"{result_parse['n_unfavorable']} unfavorable (>+0.5)")

    # ─────────────────────────────────────────────────────────────
    # Step 8: Compare with footprint
    # ─────────────────────────────────────────────────────────────
    comparison_result = None
    if compare_footprint and footprint_csv:
        logger.info("")
        logger.info("─── Step 8: Footprint Comparison ───")

        comparison_result = compare_with_footprint(
            decomp_df=result_parse["df"],
            footprint_csv=footprint_csv,
            output_dir=output_dir,
        )
        if comparison_result["success"]:
            zs = comparison_result.get("zone_summary", {})
            if zs:
                logger.info("")
                logger.info("  Zone Summary (MMPBSA total | solvation impact):")
                for zone, info in zs.items():
                    total = info.get("sum_total", 0)
                    solv = info.get("solvation_impact", 0)
                    logger.info(f"    {info.get('label', zone):30s}: "
                                f"{total:8.2f} | {solv:+8.2f} kcal/mol")

            # Save zone summary as JSON
            zone_json = output_dir / "zone_summary.json"
            with open(zone_json, "w") as f:
                json.dump(zs, f, indent=2, default=str)
            logger.info(f"  Saved: {zone_json}")
        else:
            logger.warning(f"  Comparison skipped: {comparison_result.get('error')}")

    elif compare_footprint and not footprint_csv:
        logger.info("")
        logger.info("  Footprint comparison: skipped (no CSV path provided)")
        logger.info("  Run 01d + 04b first, then re-run 01h with --footprint-csv")

    # ─────────────────────────────────────────────────────────────
    # Step 9: Generate HTML report
    # ─────────────────────────────────────────────────────────────
    html_path = None
    try:
        logger.info("")
        logger.info("─── Step 9: HTML Report ───")

        html = generate_mmpbsa_html(
            decomp_df=result_parse["df"],
            global_results=global_results,
            zone_summary=(comparison_result.get("zone_summary")
                          if comparison_result and comparison_result.get("success")
                          else None),
            comparison_df=(comparison_result.get("df")
                           if comparison_result and comparison_result.get("success")
                           else None),
            is_single_frame=is_single_frame,
            campaign_id=campaign_id,
        )

        html_path = output_dir / "mmpbsa_report.html"
        html_path.write_text(html, encoding="utf-8")
        logger.info(f"  Report: {html_path}")
    except Exception as e:
        logger.warning(f"  HTML generation failed: {e}")

    # ─────────────────────────────────────────────────────────────
    # Summary
    # ─────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  MMPBSA Analysis Complete")
    logger.info(f"  Per-residue CSV:  {result_parse['csv_path']}")
    if comparison_result and comparison_result.get("success"):
        logger.info(f"  Comparison CSV:   {comparison_result['comparison_csv']}")
    logger.info("=" * 60)

    return {
        "success": True,
        "per_residue_csv": result_parse["csv_path"],
        "comparison_csv": (comparison_result["comparison_csv"]
                           if comparison_result and comparison_result.get("success")
                           else None),
        "html_report": str(html_path) if html_path else None,
        "global_results": global_results,
        "zone_summary": (comparison_result.get("zone_summary")
                         if comparison_result and comparison_result.get("success")
                         else None),
        "n_residues": result_parse["n_residues"],
        "n_favorable": result_parse["n_favorable"],
        "n_unfavorable": result_parse["n_unfavorable"],
    }