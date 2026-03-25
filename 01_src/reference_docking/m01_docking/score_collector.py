"""
Score Collector - Core Module (01d)
=====================================
Parses DOCK6 scored mol2 files, extracts ALL header fields, ranks molecules,
and produces a clean CSV/Excel ready for downstream enrichment.

DOCK6 scored mol2 header fields (all captured):
    ##########                                Name:                 MOL
    ##########                    Molecular_Weight:             410.462
    ##########                DOCK_Rotatable_Bonds:                  15
    ##########                       Formal_Charge:                   0
    ##########                     HBond_Acceptors:                  10
    ##########                        HBond_Donors:                   6
    ##########                         Heavy_Atoms:                  28
    ##########                          Grid_Score:          -54.776524
    ##########                     Grid_vdw_energy:          -46.697788
    ##########                      Grid_es_energy:           -8.078735
    ##########           Internal_energy_repulsive:           17.656584

Multiple poses per file: each pose has its own ## block + MOLECULE block.
Best pose = lowest Grid_Score (most negative = best binding).

This module is DOCK6-specific. Molecular property enrichment (MW, LogP,
QED, PAINS, etc.) belongs to the external molecular_metrics package.

Location: 01_src/reference_docking/m01_docking/score_collector.py
Project: reference_docking
Module: 01d (DOCK6 engine)
Version: 2.1 — keep_all_poses for clustering (2026-03-21)
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# SCORED MOL2 PARSING
# =============================================================================

def parse_scored_mol2(mol2_path: str) -> List[Dict[str, Any]]:
    """
    Parse a DOCK6 scored mol2 file with one or more poses.

    Each pose has:
      - ## comment lines with scores (before @<TRIPOS>MOLECULE)
      - Standard mol2 MOLECULE/ATOM/BOND blocks

    Returns:
        List of dicts, one per pose, each with:
          - pose_index: 0-based
          - scores: dict of score name -> value (float, int, or str)
          - n_atoms: int
          - mol2_block: full text of this pose (## lines + mol2)
    """
    path = Path(mol2_path)
    if not path.exists() or path.stat().st_size == 0:
        return []

    text = path.read_text()

    # Split into pose blocks by looking for ## header sections
    # Each pose starts with a cluster of ## lines
    poses = []
    current_header = []
    current_mol2 = []
    in_mol2 = False
    pose_started = False

    for line in text.split("\n"):
        if line.startswith("##########"):
            # If we were in a mol2 block, save the previous pose
            if in_mol2 and current_mol2:
                poses.append({
                    "header_lines": current_header,
                    "mol2_lines": current_mol2,
                })
                current_header = []
                current_mol2 = []
                in_mol2 = False

            current_header.append(line)
            pose_started = True

        elif "@<TRIPOS>MOLECULE" in line:
            in_mol2 = True
            current_mol2.append(line)

        elif in_mol2:
            current_mol2.append(line)

        elif pose_started and not line.strip():
            # Empty line between ## block and @<TRIPOS>MOLECULE
            pass

    # Save the last pose
    if current_mol2:
        poses.append({
            "header_lines": current_header,
            "mol2_lines": current_mol2,
        })

    # Parse each pose
    results = []
    for idx, pose in enumerate(poses):
        scores = {}
        for line in pose["header_lines"]:
            # Format: ##########                    Key_Name:          value
            if ":" in line:
                parts = line.split(":", 1)
                key = parts[0].replace("#", "").strip()
                raw_val = parts[1].strip()
                # Try int first, then float, then keep as string
                try:
                    val = int(raw_val)
                except ValueError:
                    try:
                        val = float(raw_val)
                    except ValueError:
                        val = raw_val
                scores[key] = val

        # Count atoms
        n_atoms = 0
        in_atom = False
        for line in pose["mol2_lines"]:
            if "@<TRIPOS>ATOM" in line:
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                n_atoms += 1

        # Build full block text
        full_block = "\n".join(pose["header_lines"]) + "\n\n" + "\n".join(
            pose["mol2_lines"])

        results.append({
            "pose_index": idx,
            "scores": scores,
            "n_atoms": n_atoms,
            "mol2_block": full_block,
        })

    return results


def get_best_pose(poses: List[Dict], score_key: str = "Grid_Score") -> Optional[Dict]:
    """
    Get the pose with the best (lowest) score.

    Args:
        poses: List from parse_scored_mol2()
        score_key: Score field to rank by (default: Grid_Score)

    Returns:
        Best pose dict, or None if no valid scores
    """
    valid = [p for p in poses
             if isinstance(p["scores"].get(score_key), (int, float))]
    if not valid:
        return None
    return min(valid, key=lambda p: p["scores"][score_key])


def extract_single_pose_mol2(scored_mol2: str, pose_index: int,
                              output_mol2: str) -> bool:
    """
    Extract a single pose from a multi-pose scored mol2 file.

    Args:
        scored_mol2: Path to scored mol2 file
        pose_index: 0-based index of the pose to extract
        output_mol2: Path for output single-pose mol2

    Returns:
        True if successful
    """
    poses = parse_scored_mol2(scored_mol2)
    if pose_index >= len(poses):
        return False

    Path(output_mol2).parent.mkdir(parents=True, exist_ok=True)
    with open(output_mol2, "w") as f:
        f.write(poses[pose_index]["mol2_block"])
    return True


# =============================================================================
# MAIN PIPELINE
# =============================================================================

# Columns from DOCK6 header that we expose explicitly.
# Order here determines column order in output.
_DOCK6_HEADER_COLUMNS = [
    "Grid_Score",
    "Grid_vdw_energy",
    "Grid_es_energy",
    "Internal_energy_repulsive",
    "DOCK_Rotatable_Bonds",
    "Formal_Charge",
    "Heavy_Atoms",
    "Molecular_Weight",
    "HBond_Acceptors",
    "HBond_Donors",
]


def run_score_collection(
        docking_dir: Union[str, Path],
        output_dir: Union[str, Path],
        molecules_csv: Optional[Union[str, Path]] = None,
        score_key: str = "Grid_Score",
        max_molecules: int = 500,
        extract_best_pose_mol2: bool = True,
        keep_all_poses: bool = False,
        scores_filename: str = "dock6_scores.xlsx",
        mol2_dirname: str = "best_poses",
        source_label: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Collect and rank DOCK6 docking scores.

    Scans docking_dir for {name}/{name}_scored.mol2 files.
    Extracts best pose per molecule, ranks by score_key,
    merges SMILES from unique_molecules.csv, and produces
    Excel + CSV output ready for molecular_metrics enrichment.

    When keep_all_poses=True, also produces dock6_all_poses.csv
    with one row per pose per molecule (for clustering module 03a).

    Output columns:
        Rank, Name, Smile, Grid_Score, Grid_vdw_energy, Grid_es_energy,
        Internal_energy_repulsive, DOCK_Rotatable_Bonds, Formal_Charge,
        Heavy_Atoms, Molecular_Weight, HBond_Acceptors, HBond_Donors,
        n_poses, [any extra DOCK6 header fields]

    Args:
        docking_dir: Path to 01c_dock6_run output
        output_dir: Path for score collection output
        molecules_csv: Path to unique_molecules.csv (for SMILES merge)
        score_key: Score field to rank by
        max_molecules: Top N to include (0 = all)
        extract_best_pose_mol2: Extract best pose mol2 per molecule
        keep_all_poses: Export all poses to dock6_all_poses.csv (for clustering)
        scores_filename: Name for Excel output file
        mol2_dirname: Directory name for best pose mol2 files
        source_label: Optional label for source column

    Returns:
        Dict with: success, n_molecules, n_scored, output files
    """
    docking_dir = Path(docking_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not docking_dir.exists():
        return {"success": False,
                "error": f"Docking dir not found: {docking_dir}"}

    logger.info("=" * 60)
    logger.info("  Score Collection (DOCK6) v2.0")
    logger.info("=" * 60)
    logger.info(f"  Docking dir: {docking_dir}")
    logger.info(f"  Score key:   {score_key}")
    logger.info(f"  Output:      {output_dir}")

    # --- Find scored mol2 files ---
    mol_dirs = sorted([
        d for d in docking_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    ])

    all_results = []
    all_poses_results = []  # For keep_all_poses: one row per pose

    for mol_dir in mol_dirs:
        name = mol_dir.name
        scored_file = mol_dir / f"{name}_scored.mol2"

        if not scored_file.exists() or scored_file.stat().st_size == 0:
            logger.debug(f"  {name}: no scored mol2")
            continue

        # Parse all poses
        poses = parse_scored_mol2(str(scored_file))
        if not poses:
            logger.warning(f"  {name}: could not parse scored mol2")
            continue

        # Get best pose
        best = get_best_pose(poses, score_key)
        if best is None:
            logger.warning(f"  {name}: no valid {score_key}")
            continue

        scores = best["scores"]

        # Build result with ALL header fields
        result = {
            "Name": name,
            "n_poses": len(poses),
            "best_pose_index": best["pose_index"],
            "n_atoms": best["n_atoms"],
            "scored_mol2": str(scored_file),
        }

        # Add every key from the DOCK6 header (numeric and string)
        for key, val in scores.items():
            if key == "Name":
                continue  # We use directory name, not mol2 header name
            result[key] = val

        all_results.append(result)

        # --- Collect ALL poses for clustering ---
        if keep_all_poses:
            for pose in poses:
                pose_result = {
                    "Name": name,
                    "pose_index": pose["pose_index"],
                    "n_atoms": pose["n_atoms"],
                    "scored_mol2": str(scored_file),
                }
                for key, val in pose["scores"].items():
                    if key == "Name":
                        continue
                    pose_result[key] = val
                all_poses_results.append(pose_result)

        score_val = scores.get(score_key, "N/A")
        score_str = f"{score_val:.3f}" if isinstance(score_val, float) else str(score_val)
        logger.info(f"  {name}: {score_key}={score_str} ({len(poses)} poses)")

    if not all_results:
        logger.error("No scored molecules found!")
        return {"success": False, "error": "No scored molecules found",
                "n_molecules": 0, "n_scored": 0}

    # --- Build DataFrame and rank ---
    df = pd.DataFrame(all_results)
    df = df.sort_values(score_key, ascending=True).reset_index(drop=True)
    df.insert(0, "Rank", range(1, len(df) + 1))

    # --- Merge SMILES from unique_molecules.csv ---
    if molecules_csv and Path(molecules_csv).exists():
        try:
            df_mol = pd.read_csv(molecules_csv)

            # Normalize name column
            if "name" in df_mol.columns:
                df_mol = df_mol.rename(columns={"name": "Name"})

            if "Name" in df_mol.columns:
                # Identify SMILES column
                smiles_col = None
                for candidate in ["SMILES", "Smile", "smiles", "Smi",
                                  "canonical_smiles"]:
                    if candidate in df_mol.columns:
                        smiles_col = candidate
                        break

                if smiles_col:
                    # Merge only SMILES (rename to Smile for consistency)
                    merge_df = df_mol[["Name", smiles_col]].copy()
                    merge_df = merge_df.rename(columns={smiles_col: "Smile"})
                    df = df.merge(merge_df, on="Name", how="left")
                    n_with_smiles = df["Smile"].notna().sum()
                    logger.info(
                        f"  Merged SMILES from {molecules_csv} "
                        f"({n_with_smiles}/{len(df)} matched)")
                else:
                    logger.warning(
                        f"  No SMILES column found in {molecules_csv}. "
                        f"Available: {list(df_mol.columns)[:10]}")
            else:
                logger.warning(
                    f"  No Name column found in {molecules_csv}")
        except Exception as e:
            logger.warning(f"  Could not merge molecules_csv: {e}")

    # --- Apply max_molecules limit ---
    if max_molecules > 0 and len(df) > max_molecules:
        df = df.head(max_molecules)
        logger.info(f"  Trimmed to top {max_molecules}")

    # --- Add source label ---
    if source_label:
        df["Source"] = source_label

    # --- Extract best pose mol2 files ---
    mol2_output_dir = None
    if extract_best_pose_mol2:
        mol2_output_dir = output_dir / mol2_dirname
        mol2_output_dir.mkdir(parents=True, exist_ok=True)

        n_extracted = 0
        for _, row in df.iterrows():
            name = row["Name"]
            scored_file = row.get("scored_mol2", "")
            pose_idx = row.get("best_pose_index", 0)

            if scored_file and Path(scored_file).exists():
                out_mol2 = mol2_output_dir / f"{name}.mol2"
                if extract_single_pose_mol2(scored_file, int(pose_idx),
                                             str(out_mol2)):
                    n_extracted += 1

        logger.info(
            f"  Extracted {n_extracted} best pose mol2 → {mol2_output_dir}/")

    # --- Organize columns for export ---
    export_df = df.drop(
        columns=["scored_mol2", "best_pose_index"], errors="ignore")

    # Desired column order: Rank, Name, Smile, DOCK6 scores, metadata
    ordered_cols = ["Rank", "Name"]

    # Smile right after Name (if available)
    if "Smile" in export_df.columns:
        ordered_cols.append("Smile")

    # DOCK6 header columns in canonical order
    for col in _DOCK6_HEADER_COLUMNS:
        if col in export_df.columns and col not in ordered_cols:
            ordered_cols.append(col)

    # n_poses and any remaining columns
    remaining = [c for c in export_df.columns if c not in ordered_cols]
    ordered_cols.extend(remaining)

    export_df = export_df[ordered_cols]

    # --- Save Excel ---
    xlsx_path = output_dir / scores_filename
    export_df.to_excel(str(xlsx_path), index=False, sheet_name="DOCK6_Scores")
    logger.info(f"  Saved: {xlsx_path}")

    # --- Save CSV (always, for programmatic / molecular_metrics input) ---
    csv_path = output_dir / "dock6_scores.csv"
    export_df.to_csv(str(csv_path), index=False)
    logger.info(f"  Saved: {csv_path}")

    # --- Save all-poses CSV (for clustering module 03a) ---
    all_poses_csv = None
    if keep_all_poses and all_poses_results:
        df_all = pd.DataFrame(all_poses_results)
        df_all = df_all.sort_values(
            ["Name", score_key], ascending=[True, True]
        ).reset_index(drop=True)

        all_poses_csv = output_dir / "dock6_all_poses.csv"
        df_all.to_csv(str(all_poses_csv), index=False)

        n_molecules = df_all["Name"].nunique()
        n_total_poses = len(df_all)
        logger.info(
            f"  Saved: {all_poses_csv} "
            f"({n_total_poses} poses from {n_molecules} molecules)")
    elif keep_all_poses:
        logger.warning("  keep_all_poses=True but no poses collected")

    # --- Save summary TXT ---
    summary_path = output_dir / "score_collection_summary.txt"
    w = 70
    lines = [
        "=" * w,
        "01d SCORE COLLECTION - SUMMARY (DOCK6) v2.0",
        "=" * w,
        "",
        f"Date:              {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Docking dir:       {docking_dir}",
        f"Score key:         {score_key}",
        f"Molecules scored:  {len(df)}",
        "",
        f"Best score:        {df[score_key].min():.3f} ({df.iloc[0]['Name']})",
        f"Worst score:       {df[score_key].max():.3f} ({df.iloc[-1]['Name']})",
        f"Mean score:        {df[score_key].mean():.3f}",
        f"Median score:      {df[score_key].median():.3f}",
        "",
        "Columns exported:",
        f"  {', '.join(ordered_cols)}",
        "",
        "-" * w,
        f"{'Rank':<6} {'Name':<30} {score_key:>12} {'vdW':>10} "
        f"{'ES':>10} {'FC':>4} {'DRB':>4} {'HA':>4} {'Poses':>6}",
        "-" * w,
    ]

    for _, row in export_df.iterrows():
        lines.append(
            f"{int(row['Rank']):<6} {row['Name']:<30} "
            f"{row.get('Grid_Score', 0):>12.3f} "
            f"{row.get('Grid_vdw_energy', 0):>10.3f} "
            f"{row.get('Grid_es_energy', 0):>10.3f} "
            f"{row.get('Formal_Charge', ''):>4} "
            f"{row.get('DOCK_Rotatable_Bonds', ''):>4} "
            f"{row.get('Heavy_Atoms', ''):>4} "
            f"{row.get('n_poses', 0):>6}"
        )

    lines.extend(["", "=" * w])

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    logger.info(f"  Saved: {summary_path}")

    # --- Save report JSON ---
    report = {
        "timestamp": datetime.now().isoformat(),
        "version": "2.0",
        "docking_dir": str(docking_dir),
        "score_key": score_key,
        "n_scored": len(df),
        "best_score": float(df[score_key].min()),
        "best_molecule": df.iloc[0]["Name"],
        "worst_score": float(df[score_key].max()),
        "mean_score": float(df[score_key].mean()),
        "median_score": float(df[score_key].median()),
        "columns": ordered_cols,
        "output_xlsx": str(xlsx_path),
        "output_csv": str(csv_path),
        "all_poses_csv": str(all_poses_csv) if all_poses_csv else None,
        "mol2_dir": str(mol2_output_dir) if mol2_output_dir else None,
    }

    report_path = output_dir / "score_collection_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  {len(df)} molecules scored and ranked")
    logger.info(
        f"  Best: {df.iloc[0]['Name']} "
        f"({score_key}={df[score_key].min():.3f})")
    logger.info("=" * 60)

    return {
        "success": True,
        "n_molecules": len(mol_dirs),
        "n_scored": len(df),
        "best_molecule": df.iloc[0]["Name"],
        "best_score": float(df[score_key].min()),
        "output_xlsx": str(xlsx_path),
        "output_csv": str(csv_path),
        "all_poses_csv": str(all_poses_csv) if all_poses_csv else None,
        "summary_txt": str(summary_path),
        "report_json": str(report_path),
        "mol2_dir": str(mol2_output_dir) if mol2_output_dir else None,
    }