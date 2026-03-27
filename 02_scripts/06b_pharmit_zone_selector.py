#!/usr/bin/env python3
"""
06b Pharmit Zone Selector - CLI (v3.1)
=========================================
Maps Pharmit features to zones (by ligand atom) and annotates with
receptor-side evidence (PLIP + footprint energy).

Usage:
    python 02_scripts/06b_pharmit_zone_selector.py \\
        --config 03_configs/06b_pharmit_zone_selector.yaml \\
        --campaigns 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml \\
        --pharmit 04_data/campaigns/UDX_reference_pH63/reference/pharmit.json \\
        --plip 05_results/.../03a_plip_analysis/interactions.json \\
        --footprint ../reference_docking/05_results/.../residue_consensus.csv

Project: reference_docking
Module: 06b
Version: 3.1
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-8s | %(message)s")

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from reference_docking.m06_pharmit.pharmit_zone_selector import run_pharmit_zone_selector

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode="w", encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="06b Pharmit Zone Selector — zone mapping + evidence annotation",
    )
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaigns", type=str, default=None)
    parser.add_argument("--pharmit", type=str, default=None,
                        help="Pharmit-exported JSON")
    parser.add_argument("--plip", type=str, default=None,
                        help="PLIP interactions JSON (for evidence)")
    parser.add_argument("--footprint", type=str, default=None,
                        help="Footprint residue_consensus.csv (for evidence)")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    campaign_id = ""
    campaign_dir = None
    if args.campaigns:
        cc = load_yaml(args.campaigns)
        campaign_dir = Path(args.campaigns).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

    # --- Pharmit JSON ---
    pharmit_json = args.pharmit
    if not pharmit_json and campaign_dir:
        candidate = campaign_dir / "reference" / "pharmit.json"
        if candidate.exists():
            pharmit_json = str(candidate)
    if not pharmit_json or not Path(pharmit_json).exists():
        logger.error(f"Pharmit JSON not found: {pharmit_json}")
        return 1

    # --- PLIP (optional evidence) ---
    plip_json = args.plip
    if not plip_json and campaign_id:
        candidate = Path("05_results") / campaign_id / "03a_plip_analysis" / "interactions.json"
        if candidate.exists():
            plip_json = str(candidate)

    # --- Footprint (optional evidence) ---
    footprint_csv = args.footprint
    if not footprint_csv:
        footprint_csv = params.get("footprint_csv")
    if not footprint_csv and campaign_id:
        candidate = Path("05_results") / campaign_id / "04_dock6_analysis" / "04b_footprint_analysis" / "residue_consensus.csv"
        if candidate.exists():
            footprint_csv = str(candidate)

    # --- Output ---
    output_dir = args.output or str(
        Path("05_results") / campaign_id / "06b_pharmit_zones"
    ) if campaign_id else "06b_pharmit_zones"

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    log_path = Path(output_dir) / "06b_pharmit_zones.log"
    setup_log_file(log_path, log_level)

    # --- Execute ---
    result = run_pharmit_zone_selector(
        pharmit_json_path=pharmit_json,
        output_dir=output_dir,
        campaign_id=campaign_id,
        plip_json_path=plip_json,
        footprint_csv_path=footprint_csv,
        max_distance=params.get("max_distance", 2.0),
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info("")
    for sname, spath in result["variants"].items():
        logger.info(f"Upload to Pharmit: {spath}")

    return 0


if __name__ == "__main__":
    sys.exit(main())