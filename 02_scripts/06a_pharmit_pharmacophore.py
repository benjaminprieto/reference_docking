#!/usr/bin/env python3
"""
06a Pharmit Pharmacophore Generator - CLI (v5.0)
===================================================
Generate ranked sub-pocket Pharmit queries from DOCK6 footprint energy
and PLIP interactions.

v5.0: PLIP + Footprint IS the pharmacophore. Footprint fills gaps where
PLIP doesn't detect interactions (e.g., CH-π stacking, weak vdW).

Required inputs:
    - DOCK6 residue_consensus.csv (from reference_docking 04b)
    - Reference ligand mol2 (crystal coordinates)
    - Receptor mol2 (for residue coordinates)

Optional:
    - PLIP interactions JSON (adds geometric detail)

Usage:
    # With all inputs (recommended):
    python 02_scripts/06a_pharmit_pharmacophore.py \\
        --config 03_configs/06a_pharmit_pharmacophore.yaml \\
        --campaigns 04_data/campaigns/UDX_pharmit_pH63/campaign_config.yaml \\
        --footprint path/to/residue_consensus.csv \\
        --plip path/to/interactions.json

    # Minimal (footprint only, no PLIP):
    python 02_scripts/06a_pharmit_pharmacophore.py \\
        --config 03_configs/06a_pharmit_pharmacophore.yaml \\
        --campaigns 04_data/campaigns/UDX_pharmit_pH63/campaign_config.yaml \\
        --footprint path/to/residue_consensus.csv

    # Direct paths (no campaign config):
    python 02_scripts/06a_pharmit_pharmacophore.py \\
        --config 03_configs/06a_pharmit_pharmacophore.yaml \\
        --footprint results/residue_consensus.csv \\
        --ligand molecules/UDX.mol2 \\
        --receptor results/rec_charged.mol2 \\
        --output results/06a_pharmit/

Project: reference_docking
Module: 06a
Version: 5.0 — Footprint + PLIP (2026-03-26)
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

from reference_docking.m06_pharmit.pharmit_pharmacophore import generate_pharmit_pharmacophore

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
        description="06a Pharmit Pharmacophore — PLIP + Footprint → ranked Pharmit JSONs",
    )
    parser.add_argument("--config", "-c", type=str, required=True,
                        help="Module config YAML")
    parser.add_argument("--campaigns", type=str, default=None,
                        help="Campaign config YAML")

    # Primary inputs
    parser.add_argument("--footprint", type=str, default=None,
                        help="DOCK6 residue_consensus.csv (from reference_docking 04b)")
    parser.add_argument("--plip", type=str, default=None,
                        help="PLIP interactions JSON (optional)")

    # Explicit paths (override campaign auto-discovery)
    parser.add_argument("--ligand", type=str, default=None,
                        help="Reference ligand mol2 (crystal coordinates)")
    parser.add_argument("--receptor", type=str, default=None,
                        help="Receptor mol2 (charged, from 00b)")
    parser.add_argument("--output", "-o", type=str, default=None)

    # Strategy filter
    parser.add_argument("--strategies", type=str, default=None,
                        help="Comma-separated strategies (xylose,uracil,combined,analogues,druglike)")

    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    mc = load_yaml(args.config)
    params = mc.get("parameters", {})

    # Campaign-based path resolution
    campaign_id = None
    campaign_dir = None
    if args.campaigns:
        cc = load_yaml(args.campaigns)
        campaign_dir = Path(args.campaigns).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

    # --- Footprint CSV (required) ---
    footprint_csv = args.footprint
    if not footprint_csv and campaign_id:
        # Auto-discover from reference_docking results
        candidates = [
            # reference_docking output
            Path(params.get("reference_docking_results", "")) / "04_dock6_analysis" / "04b_footprint_analysis" / "residue_consensus.csv",
            # Same campaign in current project
            Path("05_results") / campaign_id / "04_dock6_analysis" / "04b_footprint_analysis" / "residue_consensus.csv",
            # Explicit path in config
            Path(params.get("footprint_csv", "")),
        ]
        for c in candidates:
            if c.exists():
                footprint_csv = str(c)
                break

    if not footprint_csv or not Path(footprint_csv).exists():
        logger.error("Footprint CSV not found. Provide --footprint path/to/residue_consensus.csv")
        logger.error("This comes from reference_docking module 04b.")
        return 1

    # --- Ligand mol2 ---
    ligand_mol2 = args.ligand
    if not ligand_mol2 and campaign_dir:
        cc = load_yaml(args.campaigns)
        ref_key = cc.get("reference_mol2", "reference/UDX.mol2")
        ref_path = campaign_dir / ref_key
        if ref_path.exists():
            ligand_mol2 = str(ref_path)
        else:
            # Try reference/ dir
            ref_dir = campaign_dir / "reference"
            if ref_dir.exists():
                mol2s = list(ref_dir.glob("*.mol2"))
                if mol2s:
                    ligand_mol2 = str(mol2s[0])

    if not ligand_mol2:
        # Try 00a output
        if campaign_id:
            lig_00a = Path("05_results") / campaign_id / "00a_ligand_preparation"
            if lig_00a.exists():
                mol2s = list(lig_00a.glob("*.mol2"))
                if mol2s:
                    ligand_mol2 = str(mol2s[0])

    if not ligand_mol2 or not Path(ligand_mol2).exists():
        logger.error(f"Ligand mol2 not found: {ligand_mol2}")
        return 1

    # --- Receptor mol2 ---
    receptor_mol2 = args.receptor
    if not receptor_mol2 and campaign_id:
        rec_path = Path("05_results") / campaign_id / "00b_receptor_preparation" / "rec_charged.mol2"
        if rec_path.exists():
            receptor_mol2 = str(rec_path)

    if not receptor_mol2 or not Path(receptor_mol2).exists():
        logger.error(f"Receptor mol2 not found: {receptor_mol2}")
        return 1

    # --- Receptor PDB (for PDB→mol2 residue mapping) ---
    receptor_pdb = None
    if campaign_id:
        rec_pdb_candidates = [
            Path("05_results") / campaign_id / "00b_receptor_preparation" / "receptor_protonated.pdb",
            Path("05_results") / campaign_id / "00b_receptor_preparation" / "rec_noH.pdb",
        ]
        for c in rec_pdb_candidates:
            if c.exists():
                receptor_pdb = str(c)
                break

    # --- PLIP JSON (optional) ---
    plip_json = args.plip
    if not plip_json and campaign_id:
        plip_candidates = [
            Path("05_results") / campaign_id / "03a_plip_analysis" / "interactions.json",
            campaign_dir / "interactions" / "interactions.json" if campaign_dir else Path(""),
        ]
        for c in plip_candidates:
            if c.exists():
                plip_json = str(c)
                break

    # --- Output ---
    output_dir = args.output or str(
        Path("05_results") / campaign_id / "06a_pharmit"
    ) if campaign_id else "06a_pharmit_output"

    # --- Ligand name ---
    ligand_name = Path(ligand_mol2).stem

    # --- Parameters ---
    n_required = params.get("n_required", 3)
    n_optional = params.get("n_optional", 2)
    radius = params.get("radius", 1.0)
    energy_cutoff = params.get("energy_cutoff", 0.0)
    footprint_cutoff = params.get("footprint_cutoff", -0.5)
    contact_distance = params.get("contact_distance", 4.5)
    include_hydrophobic = params.get("include_hydrophobic", False)
    log_level = args.log_level or params.get("log_level", "INFO")

    strategies = None
    if args.strategies:
        strategies = [s.strip() for s in args.strategies.split(",")]
    elif "strategies" in params:
        strategies = params["strategies"]

    # --- Logging ---
    logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))
    log_path = Path(output_dir) / "06a_pharmit.log"
    setup_log_file(log_path, log_level)

    # --- Execute ---
    logger.info("=" * 60)
    logger.info("  REFERENCE_DOCKING - Module 06a: Pharmit Pharmacophore")
    logger.info("=" * 60)
    logger.info(f"Campaign:     {campaign_id or 'N/A'}")
    logger.info(f"Footprint:    {Path(footprint_csv).name}")
    logger.info(f"Ligand:       {Path(ligand_mol2).name}")
    logger.info(f"Receptor:     {Path(receptor_mol2).name}")
    logger.info(f"PLIP:         {Path(plip_json).name if plip_json else 'not provided'}")
    logger.info(f"Output:       {output_dir}")

    result = generate_pharmit_pharmacophore(
        output_dir=output_dir,
        residue_csv_path=footprint_csv,
        ligand_mol2_path=ligand_mol2,
        receptor_mol2_path=receptor_mol2,
        plip_json_path=plip_json,
        receptor_pdb_path=receptor_pdb,
        output_name=params.get("output_name", "pharmacophore"),
        include_hydrophobic=include_hydrophobic,
        energy_cutoff=energy_cutoff,
        footprint_cutoff=footprint_cutoff,
        contact_distance=contact_distance,
        n_required=n_required,
        n_optional=n_optional,
        radius=radius,
        strategies=strategies,
        ligand_name=ligand_name,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info(f"\nDone: {result['n_features']} features, "
                f"{result['n_strategies']} strategies generated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
