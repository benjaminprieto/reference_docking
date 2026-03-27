#!/usr/bin/env python3
"""
01f GB/SA Hawkins Re-scoring - CLI
=====================================
Re-scores DOCK6 poses with GB/SA Hawkins implicit solvation (primary).

DOCK6.13 ignores _secondary params. This module runs GB/SA as a
separate rigid rescore with gbsa_hawkins_score_primary=yes.

Input:  05_results/{campaigns}/01c_dock6_run/{name}/{name}_scored.mol2
Output: 05_results/{campaigns}/01f_gbsa_rescore/{name}/{name}_gbsa_scored.mol2

Usage:
    python 02_scripts/01f_gbsa_rescore.py \\
        --config 03_configs/01f_gbsa_rescore.yaml \\
        --campaigns 04_data/campaigns/SD1_reference_pH63/campaign_config.yaml

Project: reference_docking
Module: 01f (DOCK6)
Version: 1.0 (2026-03-25)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from reference_docking.m01_docking.gbsa_rescore import run_gbsa_rescore

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="01f GB/SA Hawkins Re-scoring")
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaigns", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--name", type=str, default=None)
    parser.add_argument("--timeout", type=int, default=None)
    parser.add_argument("--minimize", action="store_true",
                        help="Run simplex minimization before GB/SA scoring")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    cc = load_yaml(args.campaigns)
    campaign_dir = Path(args.campaigns).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    docking_dir = str(Path("05_results") / campaign_id / "01c_dock6_run")
    output_dir = args.output or str(Path("05_results") / campaign_id / "01f_gbsa_rescore")

    # Receptor mol2
    rec_mol2_path = Path("05_results") / campaign_id / "00b_receptor_preparation" / "rec_charged.mol2"
    if not rec_mol2_path.exists():
        logger.error(f"Receptor mol2 not found: {rec_mol2_path}")
        return 1

    dock6_home = params.get("dock6_home", "/opt/dock6")
    timeout = args.timeout or params.get("timeout_per_molecule", 600)
    minimize = args.minimize or params.get("minimize", False)
    solvent_dielectric = params.get("solvent_dielectric", 78.5)
    salt_concentration = params.get("salt_concentration", 0.15)
    gb_offset = params.get("gb_offset", 0.09)
    molecule_filter = [args.name] if args.name else None

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    logger.info("=" * 60)
    logger.info("  REFERENCE_DOCKING - Module 01f: GB/SA Hawkins Rescore")
    logger.info("=" * 60)
    logger.info(f"Campaign:   {campaign_id}")
    logger.info(f"Docking:    {docking_dir}")
    logger.info(f"Receptor:   {rec_mol2_path.name}")
    logger.info(f"Minimize:   {'yes' if minimize else 'no'}")
    logger.info(f"Dielectric: {solvent_dielectric}")
    logger.info(f"Salt:       {salt_concentration} M")
    logger.info(f"Output:     {output_dir}")

    result = run_gbsa_rescore(
        docking_dir=docking_dir,
        output_dir=output_dir,
        receptor_mol2=str(rec_mol2_path.resolve()),
        dock6_home=dock6_home,
        timeout_per_molecule=timeout,
        molecule_filter=molecule_filter,
        minimize=minimize,
        solvent_dielectric=solvent_dielectric,
        salt_concentration=salt_concentration,
        gb_offset=gb_offset,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0 if result["n_failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
