#!/usr/bin/env python3
"""
04b DOCK6 Footprint Analysis - CLI
======================================
Per-residue vdW + ES energy decomposition via DOCK6 footprint scoring.

Phase 2 only: Parse per-residue footprint data from 01d TXT files.
Builds mol2→PDB residue mapping so footprint and contacts share numbering.

Input:
    01d_footprint_rescore/{name}/{name}_fps_footprint_scored.txt
    00b_receptor_preparation/rec_charged.mol2  (for residue mapping)
    00b_receptor_preparation/rec_noH.pdb       (for residue mapping)

Output:
    05_results/{campaign}/04_dock6_analysis/04b_footprint_analysis/

Project: reference_docking
Module: 04b (DOCK6 analysis)
Version: 1.2 (2026-03-23) — adds residue remapping mol2→PDB
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from reference_docking.m04_dock6_analysis.footprint_analysis import run_footprint_analysis

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path, log_level="INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="04b DOCK6 Footprint Analysis — per-residue energy with PDB numbering",
    )
    parser.add_argument("--config", "-c", type=str, help="Module config YAML")
    parser.add_argument("--campaign", type=str, help="Campaign config YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None)

    args = parser.parse_args()

    # Defaults
    rescore_dir = None
    analysis_dir = None
    receptor_mol2 = None
    receptor_pdb = None
    campaign_id = "direct"
    pharmacophore_threshold = 0.8
    energy_cutoff = -0.5
    log_level = "INFO"

    # --- Campaign config ---
    if args.campaign:
        cc = load_yaml(args.campaign)
        campaign_dir = Path(args.campaign).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        base = Path("05_results") / campaign_id

        # Input: 01d footprint rescore output
        rescore_dir = str(base / "01d_footprint_rescore")

        # Output: 04_dock6_analysis
        analysis_dir = str(base / "04_dock6_analysis" / "04b_footprint_analysis")

        # Receptor files for residue mapping
        rec_mol2 = base / "00b_receptor_preparation" / "rec_charged.mol2"
        rec_pdb = base / "00b_receptor_preparation" / "rec_noH.pdb"
        if rec_mol2.exists():
            receptor_mol2 = str(rec_mol2)
        if rec_pdb.exists():
            receptor_pdb = str(rec_pdb)

    # --- Module config ---
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        pharmacophore_threshold = params.get("pharmacophore_threshold", pharmacophore_threshold)
        energy_cutoff = params.get("energy_cutoff", energy_cutoff)
        log_level = params.get("log_level", log_level)

    # --- CLI overrides ---
    if args.output:
        analysis_dir = args.output
    if args.log_level:
        log_level = args.log_level

    if not rescore_dir or not Path(rescore_dir).exists():
        logger.error(f"Footprint rescore dir not found: {rescore_dir}")
        logger.error("Run module 01d (footprint re-scoring) first.")
        return 1

    if not analysis_dir:
        analysis_dir = "05_results/04_dock6_analysis/04b_footprint_analysis"

    logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    setup_log_file(Path(analysis_dir) / "04b_footprint_analysis.log", log_level)

    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 04b: DOCK6 Footprint Analysis")
    logger.info("=" * 60)
    logger.info(f"  Campaign:      {campaign_id}")
    logger.info(f"  Footprint dir: {rescore_dir} (from 01d)")
    logger.info(f"  Receptor mol2: {receptor_mol2 or 'not found'}")
    logger.info(f"  Receptor PDB:  {receptor_pdb or 'not found'}")
    logger.info(f"  Output:        {analysis_dir}")

    result = run_footprint_analysis(
        footprint_dir=rescore_dir,
        output_dir=analysis_dir,
        receptor_mol2=receptor_mol2,
        receptor_pdb=receptor_pdb,
        pharmacophore_threshold=pharmacophore_threshold,
        energy_cutoff=energy_cutoff,
    )

    if not result.get("success"):
        logger.error(f"Analysis failed: {result.get('error')}")
        return 1

    logger.info(f"\nNext: python 02_scripts/04c_binding_modes.py "
                f"--config 03_configs/04c_binding_modes.yaml "
                f"--campaign {args.campaign or '<campaign_config.yaml>'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())