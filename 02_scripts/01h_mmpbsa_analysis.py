#!/usr/bin/env python3
"""
01h MMPBSA Analysis - CLI
============================
Parses MMPBSA.py output from 01g and compares with DOCK6 footprint.

Produces:
  - per_residue_decomp.csv:             vdW, ES, GB, SA, total (PDB numbering)
  - comparison_invacuo_vs_solvated.csv: merge with 01d/04b footprint
  - zone_summary.json:                  energy by binding site zone
  - global_binding_energy.json:         total ΔG components

Input:  05_results/{campaign}/01g_mmpbsa_decomp/mmpbsa/FINAL_DECOMP_MMPBSA.dat
Output: 05_results/{campaign}/01h_mmpbsa_analysis/

Usage:
    python 02_scripts/01h_mmpbsa_analysis.py \\
        --config 03_configs/01h_mmpbsa_analysis.yaml \\
        --campaigns 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml

    # With explicit footprint CSV:
    python 02_scripts/01h_mmpbsa_analysis.py \\
        --config 03_configs/01h_mmpbsa_analysis.yaml \\
        --campaigns 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml \\
        --footprint-csv 05_results/UDX_reference_pH63/04b_footprint_analysis/residue_consensus.csv

Project: reference_docking
Module: 01h
Version: 1.0 (2026-03-26)
"""

import argparse
import json
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from reference_docking.m01_docking.mmpbsa_analysis import run_mmpbsa_analysis

logger = logging.getLogger(__name__)


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
        description="01h MMPBSA Analysis — parse decomposition + footprint comparison",
    )
    parser.add_argument("--config", "-c", type=str, required=True, help="Module YAML")
    parser.add_argument("--campaigns", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--decomp-dir", type=str, default=None,
                        help="Override 01g output directory")
    parser.add_argument("--footprint-csv", type=str, default=None,
                        help="Explicit path to residue_consensus.csv from 04b")
    parser.add_argument("--no-compare-footprint", action="store_true",
                        help="Skip footprint comparison")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    # =========================================================================
    # LOAD CONFIGS
    # =========================================================================
    cc = load_yaml(args.campaigns)
    campaign_dir = Path(args.campaigns).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    # =========================================================================
    # RESOLVE PATHS
    # =========================================================================
    results_base = Path("05_results") / campaign_id
    output_subdir = mc.get("outputs", {}).get("subdir", "01h_mmpbsa_analysis")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    # 01g output directory
    decomp_base = Path(args.decomp_dir) if args.decomp_dir else results_base / "01g_mmpbsa_decomp"

    # Find MMPBSA output files
    decomp_file = decomp_base / "mmpbsa" / "FINAL_DECOMP_MMPBSA.dat"
    results_file = decomp_base / "mmpbsa" / "FINAL_RESULTS_MMPBSA.dat"

    if not decomp_file.exists():
        logger.error(f"Decomp file not found: {decomp_file}")
        logger.error("Run module 01g first.")
        return 1

    # Receptor PDB (from 00b)
    rec_pdb = results_base / "00b_receptor_preparation" / "receptor_protonated.pdb"
    if not rec_pdb.exists():
        rec_pdb = results_base / "00b_receptor_preparation" / "rec_noH.pdb"
    if not rec_pdb.exists():
        logger.error("Receptor PDB not found in 00b_receptor_preparation/")
        return 1

    # Detect mode from 01g pipeline log
    is_single_frame = True
    pipeline_log_path = decomp_base / "01g_pipeline_log.json"
    if pipeline_log_path.exists():
        with open(pipeline_log_path) as f:
            plog = json.load(f)
        is_single_frame = (plog.get("mode") == "single_point")
        logger.info(f"Detected 01g mode: {plog.get('mode')} "
                     f"({plog.get('n_frames', '?')} frames)")

    # Auto-detect multi-frame from DECOMP Std.Dev values
    if is_single_frame and decomp_file.exists():
        with open(decomp_file) as f:
            decomp_text = f.read()
        in_deltas = False
        for line in decomp_text.split("\n"):
            if "DELTAS:" in line:
                in_deltas = True
                continue
            if in_deltas and line.strip() and not line.startswith(",") and not line.startswith("Residue"):
                parts = line.split(",")
                if len(parts) >= 4:
                    try:
                        std_dev = float(parts[3])
                        if std_dev > 0.001:
                            is_single_frame = False
                            logger.info(f"  Auto-detected multi-frame from DECOMP (Std.Dev > 0)")
                            break
                    except (ValueError, IndexError):
                        continue

    # Footprint CSV (auto-resolve from 04b)
    footprint_csv = args.footprint_csv
    compare_fp = params.get("compare_footprint", True) and not args.no_compare_footprint

    if not footprint_csv and compare_fp:
        fp_candidates = [
            results_base / "04b_footprint_analysis" / "residue_consensus.csv",
            results_base / "04_dock6_analysis" / "04b_footprint_analysis" / "residue_consensus.csv",
        ]
        for fp in fp_candidates:
            if fp.exists():
                footprint_csv = str(fp)
                break

    # =========================================================================
    # SETUP
    # =========================================================================
    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))
    setup_log_file(output_dir / "01h_mmpbsa_analysis.log", log_level)

    logger.info("=" * 60)
    logger.info("  REFERENCE_DOCKING - Module 01h: MMPBSA Analysis")
    logger.info("=" * 60)
    logger.info(f"Campaign:   {campaign_id}")
    logger.info(f"Decomp:     {decomp_file}")
    logger.info(f"Receptor:   {rec_pdb.name}")
    logger.info(f"Footprint:  {footprint_csv or 'not found'}")
    logger.info(f"Output:     {output_dir}")

    # =========================================================================
    # EXECUTE
    # =========================================================================
    result = run_mmpbsa_analysis(
        decomp_file=str(decomp_file),
        results_file=str(results_file),
        receptor_pdb=str(rec_pdb),
        output_dir=str(output_dir),
        is_single_frame=is_single_frame,
        compare_footprint=compare_fp,
        footprint_csv=footprint_csv,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
