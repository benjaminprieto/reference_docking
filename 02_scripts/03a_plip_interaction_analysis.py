#!/usr/bin/env python3
"""
03a PLIP Interaction Analysis - CLI
======================================
Analyze crystal protein-ligand interactions using PLIP.

Usage:
    # With campaigns (auto-resolves paths):
    python 02_scripts/03a_plip_interaction_analysis.py \
        --campaigns 04_data/campaigns/SD1_reference_pH63/campaign_config.yaml \
        --config 03_configs/03a_plip_interaction_analysis.yaml

Project: reference_docking
Module: 03a
Version: 1.0

Dependencies: plip (conda install -c conda-forge plip)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s | %(levelname)-8s | %(message)s")

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from reference_docking.m03_crystal_analysis.plip_interaction_analysis import (
    run_plip_analysis,
)

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
        description="03a PLIP Interaction Analysis — crystal protein-ligand interactions",
    )

    # Inputs
    parser.add_argument("--receptor", "-r", type=str, default=None,
                        help="Protonated receptor PDB (from 00b)")
    parser.add_argument("--ligand", "-l", type=str, default=None,
                        help="Protonated ligand SDF/mol2 (from 00c)")
    parser.add_argument("--ligand-name", type=str, default=None,
                        help="Name for the ligand (default: filename stem)")

    # Campaign mode
    parser.add_argument("--campaigns", type=str, default=None)
    parser.add_argument("--config", "-c", type=str, default=None)

    # Output
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--name", type=str, default="interactions")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # --- Config defaults ---
    output_name = args.name
    log_level = "INFO"

    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        output_name = params.get("output_name", output_name)
        log_level = params.get("log_level", log_level)

    if args.log_level:
        log_level = args.log_level
    logging.getLogger().setLevel(getattr(logging, log_level))

    # --- Resolve paths ---
    receptor_path = args.receptor
    ligand_path = args.ligand
    ligand_name = args.ligand_name
    output_dir = args.output

    if args.campaigns:
        cc = load_yaml(args.campaigns)
        campaign_dir = Path(args.campaigns).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)
        results_base = Path("05_results") / campaign_id
        docking_ph = cc.get("docking_ph", 7.2)

        # Receptor: protonated from 00b
        if receptor_path is None:
            candidates = [
                results_base / "00b_receptor_preparation" / "receptor_protonated.pdb",
                results_base / "00b_receptor_preparation" / "receptor_clean.pdb",
            ]
            for c in candidates:
                if c.exists():
                    receptor_path = str(c)
                    logger.info(f"  Receptor: {c}")
                    break

        # Ligand: reference mol2 with crystallographic coordinates
        # PLIP needs exact crystal coords — NOT protonated SDF from 00c
        # (00c may regenerate coordinates). PLIP protonates internally.
        if ligand_path is None:
            ref_mol2 = cc.get("grids", {}).get("binding_site", {}).get("reference_mol2")
            ref_name = cc.get("grids", {}).get("binding_site", {}).get(
                "reference_name", "UDX")

            if ligand_name is None:
                ligand_name = ref_name

            if ref_mol2:
                candidate = campaign_dir / ref_mol2
                if candidate.exists():
                    ligand_path = str(candidate)
                    logger.info(f"  Ligand (crystal reference): {candidate}")

        # Output
        if output_dir is None:
            output_dir = str(results_base / "03a_plip_analysis")

    # --- Validate ---
    if receptor_path is None or not Path(receptor_path).exists():
        logger.error(f"Receptor not found: {receptor_path}")
        logger.error("Run 00b first, or provide --receptor")
        sys.exit(1)
    if ligand_path is None or not Path(ligand_path).exists():
        logger.error(f"Ligand not found: {ligand_path}")
        logger.error("Run 00c first, or provide --ligand")
        sys.exit(1)

    if output_dir is None:
        output_dir = "05_results/03a_plip_analysis"

    # Log file
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    setup_log_file(out_path / "03a_plip_analysis.log", log_level)

    # --- Run ---
    result = run_plip_analysis(
        receptor_pdb=receptor_path,
        ligand_path=ligand_path,
        output_dir=output_dir,
        ligand_name=ligand_name,
        output_name=output_name,
    )

    if result["success"]:
        logger.info(f"\nInteractions JSON: {result['json_path']}")
        logger.info(f"Use in 06a: --plip {result['json_path']}")
        sys.exit(0)
    else:
        logger.error(f"Failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()