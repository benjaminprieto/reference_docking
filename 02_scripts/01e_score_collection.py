#!/usr/bin/env python3
"""
01e Score Collection - CLI (DOCK6-specific)
=============================================
Collects DOCK6 docking scores, ranks molecules, and produces Excel output.

Scans 01c output for scored mol2 files, extracts ALL header fields per pose,
selects best pose per molecule, merges SMILES, and generates ranked output.

Also merges FPS footprint data from 01d if available.

Hardcoded upstream paths:
    - 05_results/{campaign_id}/01c_dock6_run/{name}/{name}_scored.mol2
    - 05_results/{campaign_id}/01d_footprint_rescore/{name}/{name}_fps_scored.mol2
    - Output: 05_results/{campaign_id}/01e_score_collection/

Usage:
    python 02_scripts/01e_score_collection.py \\
        --config 03_configs/01e_score_collection.yaml \\
        --campaigns 04_data/campaigns/SD1_reference_pH63/campaign_config.yaml

Project: reference_docking
Module: 01e (DOCK6 engine) — renumbered from 01d (2026-03-22)
Version: 3.0
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
)

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

import pandas as pd
from reference_docking.m01_docking.score_collector import run_score_collection

logger = logging.getLogger(__name__)


def _parse_mol2_scores(mol2_path: str, prefix: str) -> dict:
    """Extract score fields from a scored mol2 header (##########  lines)."""
    scores = {}
    try:
        with open(mol2_path) as f:
            for line in f:
                if line.startswith("##########"):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        key = parts[-1]
                        try:
                            val = float(parts[-2])
                            scores[f"{prefix}_{key}"] = val
                        except ValueError:
                            pass
                if "@<TRIPOS>ATOM" in line:
                    break
    except Exception:
        pass
    return scores


def _merge_rescore_scores(result: dict, output_dir: str, fps_dir: str, gbsa_dir: str):
    """Merge FPS and GBSA scores into the main score CSV/Excel."""
    output_dir = Path(output_dir)
    csv_path = output_dir / "dock6_scores.csv"
    if not csv_path.exists():
        return

    df = pd.read_csv(csv_path)
    if "Name" not in df.columns:
        return

    # Collect rescore scores per molecule
    for _, row in df.iterrows():
        name = row["Name"]

        # FPS scores from 01d
        if fps_dir:
            fps_mol2 = Path(fps_dir) / name / f"{name}_fps_scored.mol2"
            if fps_mol2.exists():
                fps_scores = _parse_mol2_scores(str(fps_mol2), "FPS")
                for k, v in fps_scores.items():
                    df.loc[df["Name"] == name, k] = v

        # GBSA scores from 01f
        if gbsa_dir:
            gbsa_mol2 = Path(gbsa_dir) / name / f"{name}_gbsa_scored.mol2"
            if gbsa_mol2.exists():
                gbsa_scores = _parse_mol2_scores(str(gbsa_mol2), "GBSA")
                for k, v in gbsa_scores.items():
                    df.loc[df["Name"] == name, k] = v

    # Save updated CSV and Excel
    df.to_csv(csv_path, index=False)
    xlsx_path = output_dir / "dock6_scores.xlsx"
    if xlsx_path.exists():
        df.to_excel(xlsx_path, index=False)

    n_fps = len([c for c in df.columns if c.startswith("FPS_")])
    n_gbsa = len([c for c in df.columns if c.startswith("GBSA_")])
    if n_fps or n_gbsa:
        logger.info(f"  Merged rescore scores: {n_fps} FPS columns, {n_gbsa} GBSA columns")


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="01e Score Collection — Collect DOCK6 scores and produce "
                    "ranked Excel/CSV (DOCK6 engine)",
    )
    parser.add_argument("--config", "-c", type=str, required=True,
                        help="Module config YAML")
    parser.add_argument("--campaigns", type=str, required=True,
                        help="Campaign config YAML")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output directory")
    parser.add_argument("--docking-dir", type=str, default=None,
                        help="Override docking directory (default: 01c output)")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # --- Load configs ---
    cc = load_yaml(args.campaigns)
    campaign_id = cc.get("campaign_id", Path(args.campaigns).parent.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    docking_dir = args.docking_dir or str(
        Path("05_results") / campaign_id / "01c_dock6_run"
    )
    fps_dir = str(Path("05_results") / campaign_id / "01d_footprint_rescore")
    gbsa_dir = str(Path("05_results") / campaign_id / "01f_gbsa_rescore")
    molecules_csv = None  # Not used in reference docking (no screening catalog)
    output_dir = args.output or str(
        Path("05_results") / campaign_id / "01e_score_collection"
    )

    # --- Params ---
    score_key = params.get("score_key", "Grid_Score")
    max_molecules = params.get("max_molecules", 0)
    extract_best = params.get("extract_best_pose_mol2", True)
    keep_all_poses = params.get("keep_all_poses", False)
    log_level = args.log_level or params.get("log_level", "INFO")

    # --- Output filenames from campaigns config ---
    oc = cc.get("output", {})
    scores_filename = oc.get("scores_filename", "dock6_scores.xlsx")
    mol2_dirname = oc.get("mol2_dirname", "best_poses")

    # --- Setup logging ---
    logging.getLogger().setLevel(
        getattr(logging, log_level.upper(), logging.INFO))
    log_path = Path(output_dir) / "01e_score_collection.log"
    setup_log_file(log_path, log_level)

    # --- Execute ---
    logger.info("=" * 60)
    logger.info("  REFERENCE_DOCKING - Module 01e: Score Collection (DOCK6)")
    logger.info("=" * 60)
    logger.info(f"Campaign:     {campaign_id}")
    logger.info(f"Docking dir:  {docking_dir}")
    logger.info(f"FPS dir:      {fps_dir}")
    logger.info(f"GBSA dir:     {gbsa_dir}")
    logger.info(f"Score key:    {score_key}")
    logger.info(f"Output:       {output_dir}")

    result = run_score_collection(
        docking_dir=docking_dir,
        output_dir=output_dir,
        molecules_csv=molecules_csv,
        score_key=score_key,
        max_molecules=max_molecules,
        extract_best_pose_mol2=extract_best,
        keep_all_poses=keep_all_poses,
        scores_filename=scores_filename,
        mol2_dirname=mol2_dirname,
        source_label=cc.get("metadata", {}).get("source"),
    )

    # --- Merge FPS and GBSA scores from rescore steps ---
    if result.get("success"):
        _merge_rescore_scores(result, output_dir, fps_dir, gbsa_dir)

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info("")
    logger.info("DOCK6 score collection complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
