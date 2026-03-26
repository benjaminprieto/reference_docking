#!/usr/bin/env python3
"""
00a Ligand Preparation - CLI (Reference Docking)
===================================================
Prepares crystallographic ligand for DOCK6 reference scoring.
Preserves crystal coordinates exactly.

Usage:
    python 02_scripts/00a_ligand_preparation.py \\
        --config 03_configs/00a_ligand_preparation.yaml \\
        --campaign 04_data/campaigns/SD1_reference_pH63/campaign_config.yaml

Project: reference_docking
Module: 00a
Version: 1.0 (2026-03-25)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from reference_docking.m00_preparation.ligand_preparation import run_ligand_preparation

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="00a Ligand Preparation — crystal ligand for reference docking")
    parser.add_argument("--config", "-c", type=str, required=True)
    parser.add_argument("--campaign", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--strategy", type=str, choices=["direct", "extract", "inject"], default=None)
    parser.add_argument("--log-level", type=str, default=None, choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    cc = load_yaml(args.campaign)
    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    # Resolve paths
    output_dir = args.output or str(Path("05_results") / campaign_id / "00a_ligand_preparation")

    # Reference mol2
    ref_key = cc.get("reference_mol2", "reference/UDX.mol2")
    reference_mol2 = str(campaign_dir / ref_key) if (campaign_dir / ref_key).exists() else None
    if not reference_mol2:
        # Auto-detect in reference/
        ref_dir = campaign_dir / "reference"
        if ref_dir.exists():
            mol2s = list(ref_dir.glob("*.mol2"))
            if mol2s:
                reference_mol2 = str(mol2s[0])

    # Receptor PDB (for coord validation)
    pdb_path = None
    rec_cfg = cc.get("receptor", {})
    rec_pdb = rec_cfg.get("pdb", "receptor/*.pdb")
    rec_path = campaign_dir / rec_pdb
    if rec_path.exists():
        pdb_path = str(rec_path)

    ligand_name = params.get("ligand_name", "UDX")
    strategy = args.strategy or params.get("strategy", "direct")
    docking_ph = cc.get("docking_ph", 7.2)

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    logger.info("=" * 60)
    logger.info("  REFERENCE_DOCKING - Module 00a: Ligand Preparation")
    logger.info("=" * 60)
    logger.info(f"Campaign:  {campaign_id}")
    logger.info(f"Strategy:  {strategy}")
    logger.info(f"Reference: {reference_mol2 or 'None'}")
    logger.info(f"PDB:       {pdb_path or 'None'}")
    logger.info(f"Output:    {output_dir}")

    result = run_ligand_preparation(
        reference_mol2=reference_mol2,
        pdb_path=pdb_path,
        ligand_name=ligand_name,
        output_dir=output_dir,
        strategy=strategy,
        docking_ph=docking_ph,
    )

    if not result.get("success"):
        logger.error(f"Failed: {result.get('error')}")
        return 1

    logger.info("00a complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
