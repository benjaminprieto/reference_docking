#!/usr/bin/env python3
"""
01d Footprint Re-scoring - CLI (DOCK6-specific)
==================================================
Re-scores existing DOCK6 poses with per-residue energy decomposition
(VDW + ES + H-bond footprints in Cartesian space).

v2.0 (2026-03-25): Added GB/SA Hawkins implicit solvation support.
  Enable in YAML (gbsa_hawkins: true) or CLI (--gbsa-hawkins).
  Corrects in-vacuo electrostatic artifacts for charged residues.

This is a POST-DOCKING step that reads scored mol2 from 01c.
Based on Rizzo Lab protocol and DOCK6 manual §2.11.9.

Input:  05_results/{campaign}/01c_dock6_run/{name}/{name}_scored.mol2
Output: 05_results/{campaign}/01d_footprint_rescore/{name}/{name}_fps_scored.mol2

Usage:
    # Standard (in-vacuo footprint only):
    python 02_scripts/01d_footprint_rescore.py \\
        --config 03_configs/01d_footprint_rescore.yaml \\
        --campaign 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml

    # With GB/SA Hawkins implicit solvation:
    python 02_scripts/01d_footprint_rescore.py \\
        --config 03_configs/01d_footprint_rescore.yaml \\
        --campaign 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml \\
        --gbsa-hawkins

Project: reference_docking
Module: 01d (DOCK6 engine)
Version: 2.0
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

from reference_docking.m01_docking.footprint_rescore import run_footprint_rescore

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
        description="01d Footprint Re-scoring — per-residue energy decomposition",
    )
    parser.add_argument("--config", "-c", type=str, required=True, help="Module YAML")
    parser.add_argument("--campaign", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--name", type=str, default=None, help="Single molecule")
    parser.add_argument("--timeout", type=int, default=None)
    parser.add_argument("--gbsa-hawkins", action="store_true",
                        help="Enable GB/SA Hawkins implicit solvation (overrides YAML)")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    # --- Configs ---
    cc = load_yaml(args.campaign)
    campaign_dir = Path(args.campaign).parent
    campaign_id = cc.get("campaign_id", campaign_dir.name)
    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    # --- Paths ---
    docking_dir = str(Path("05_results") / campaign_id / "01c_dock6_run")
    output_dir = args.output or str(Path("05_results") / campaign_id / "01d_footprint_rescore")

    # Receptor mol2
    rec_mol2_path = Path("05_results") / campaign_id / "00b_receptor_preparation" / "rec_charged.mol2"
    if not rec_mol2_path.exists():
        logger.error(f"Receptor mol2 not found: {rec_mol2_path}")
        return 1
    receptor_mol2 = str(rec_mol2_path.resolve())

    # Reference mol2 — check multiple locations
    reference_mol2 = None
    # 1. Top-level reference_mol2 key
    ref_key = cc.get("reference_mol2")
    # 2. Nested under grids.binding_site.reference_mol2
    if not ref_key:
        ref_key = cc.get("grids", {}).get("binding_site", {}).get("reference_mol2")
    # 3. Default fallback
    if not ref_key:
        ref_key = "reference/UDX.mol2"

    ref_mol2_path = campaign_dir / ref_key
    if ref_mol2_path.exists():
        reference_mol2 = str(ref_mol2_path.resolve())
    else:
        # Try bare reference/ directory
        ref_mol2_path = campaign_dir / "reference" / "UDX.mol2"
        if ref_mol2_path.exists():
            reference_mol2 = str(ref_mol2_path.resolve())

    if not reference_mol2:
        logger.error("Reference mol2 not found. Footprint requires a reference ligand.")
        return 1

    # --- Params ---
    dock6_home = params.get("dock6_home", "/opt/dock6")
    timeout_per_molecule = args.timeout or params.get("timeout_per_molecule", 300)
    log_level = args.log_level or params.get("log_level", "INFO")
    molecule_filter = [args.name] if args.name else None

    # GB/SA Hawkins: CLI flag overrides YAML
    gbsa_hawkins = args.gbsa_hawkins or params.get("gbsa_hawkins", False)
    solvent_dielectric = params.get("solvent_dielectric", 78.5)
    salt_concentration = params.get("salt_concentration", 0.15)
    gb_offset = params.get("gb_offset", 0.09)

    # --- Logging ---
    logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    log_path = Path(output_dir) / "01d_footprint_rescore.log"
    setup_log_file(log_path, log_level)

    # --- Execute ---
    logger.info("=" * 60)
    logger.info("  MOLECULAR_DOCKING - Module 01d: Footprint Re-scoring")
    logger.info("=" * 60)
    logger.info(f"Campaign:     {campaign_id}")
    logger.info(f"Docking dir:  {docking_dir}")
    logger.info(f"Receptor:     {Path(receptor_mol2).name}")
    logger.info(f"Reference:    {Path(reference_mol2).name}")
    logger.info(f"Output:       {output_dir}")
    logger.info(f"GB/SA Hawkins: {'YES' if gbsa_hawkins else 'no'}")

    result = run_footprint_rescore(
        docking_dir=docking_dir,
        output_dir=output_dir,
        receptor_mol2=receptor_mol2,
        reference_mol2=reference_mol2,
        dock6_home=dock6_home,
        timeout_per_molecule=timeout_per_molecule,
        molecule_filter=molecule_filter,
        gbsa_hawkins=gbsa_hawkins,
        solvent_dielectric=solvent_dielectric,
        salt_concentration=salt_concentration,
        gb_offset=gb_offset,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info(f"\nNext: python 02_scripts/01e_score_collection.py "
                f"--config 03_configs/01e_score_collection.yaml "
                f"--campaign {args.campaign}")

    return 0 if result["n_failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())