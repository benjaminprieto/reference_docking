#!/usr/bin/env python3
"""
01g MMPBSA Per-Residue Decomposition - CLI
=============================================
Prepares AMBER system and runs MMPBSA.py per-residue decomposition.

Two modes:
  single_point — 1-frame MMPBSA (fast, no error bars)
  md           — OpenMM MD → trajectory MMPBSA (mean ± std)

Input:  05_results/{campaign}/01c_dock6_run/{name}/{name}_scored.mol2
Output: 05_results/{campaign}/01g_mmpbsa_decomp/

After 01g completes, run 01h for analysis:
  python 02_scripts/01h_mmpbsa_analysis.py ...

Usage:
    # Single-point (default):
    python 02_scripts/01g_mmpbsa_decomp.py \\
        --config 03_configs/01g_mmpbsa_decomp.yaml \\
        --campaigns 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml

    # With MD:
    python 02_scripts/01g_mmpbsa_decomp.py \\
        --config 03_configs/01g_mmpbsa_decomp.yaml \\
        --campaigns 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml \\
        --mode md

    # Override production length:
    python 02_scripts/01g_mmpbsa_decomp.py \\
        --config 03_configs/01g_mmpbsa_decomp.yaml \\
        --campaigns 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml \\
        --mode md --production-ns 10.0

Project: reference_docking
Module: 01g
Version: 1.0 (2026-03-26)
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-8s | %(message)s")
sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from reference_docking.m01_docking.mmpbsa_decomp import run_mmpbsa_decomp

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
        description="01g MMPBSA Per-Residue Decomposition — preparation + execution",
    )
    parser.add_argument("--config", "-c", type=str, required=True, help="Module YAML")
    parser.add_argument("--campaigns", type=str, required=True, help="Campaign YAML")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--name", type=str, default=None,
                        help="Molecule name (default: first in 01c)")
    parser.add_argument("--mode", type=str, default=None,
                        choices=["single_point", "md"],
                        help="Override mode (default from YAML)")
    parser.add_argument("--production-ns", type=float, default=None,
                        help="MD production length in ns (overrides YAML)")
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
    output_subdir = mc.get("outputs", {}).get("subdir", "01g_mmpbsa_decomp")
    output_dir = Path(args.output) if args.output else results_base / output_subdir

    # Resolve pose source: "docking" (from 01c) or "reference" (crystal mol2)
    pose_source = cc.get("pose_source", "docking")
    molecule_name = args.name
    scored_mol2 = None

    if pose_source == "reference":
        # Use reference mol2 directly (no docking required)
        ref_key = cc.get("grids", {}).get("binding_site", {}).get("reference_mol2")
        if ref_key:
            ref_path = campaign_dir / ref_key
            if ref_path.exists():
                scored_mol2 = ref_path
        if not scored_mol2 or not scored_mol2.exists():
            # Fallback: 00a output
            lig_00a = results_base / "00a_ligand_preparation"
            mol2s = sorted(lig_00a.glob("*.mol2")) if lig_00a.exists() else []
            if mol2s:
                scored_mol2 = mol2s[0]
        if not scored_mol2 or not Path(scored_mol2).exists():
            logger.error("Reference mol2 not found. Set grids.binding_site.reference_mol2 in campaign_config.")
            return 1
        molecule_name = molecule_name or Path(scored_mol2).stem
        logger.info(f"Pose source: reference ({Path(scored_mol2).name})")
    else:
        # From 01c docking results
        docking_dir = results_base / "01c_dock6_run"
        if not molecule_name:
            if docking_dir.exists():
                mol_dirs = sorted([d for d in docking_dir.iterdir() if d.is_dir()])
                if mol_dirs:
                    molecule_name = mol_dirs[0].name
                    logger.info(f"Auto-detected molecule: {molecule_name}")
        if not molecule_name:
            logger.error("No molecule found in 01c_dock6_run/. Run 01c first.")
            return 1
        scored_mol2 = docking_dir / molecule_name / f"{molecule_name}_scored.mol2"
        if not scored_mol2.exists():
            logger.error(f"Scored mol2 not found: {scored_mol2}")
            return 1

    # Receptor mol2 (from 00b)
    rec_mol2 = results_base / "00b_receptor_preparation" / "rec_charged.mol2"
    if not rec_mol2.exists():
        logger.error(f"Receptor mol2 not found: {rec_mol2}")
        return 1

    # Receptor PDB (from 00b)
    rec_pdb = results_base / "00b_receptor_preparation" / "receptor_protonated.pdb"
    if not rec_pdb.exists():
        rec_pdb = results_base / "00b_receptor_preparation" / "rec_noH.pdb"
    if not rec_pdb.exists():
        logger.error("Receptor PDB not found in 00b_receptor_preparation/")
        return 1

    # =========================================================================
    # MERGE PARAMETERS (YAML + CLI overrides)
    # =========================================================================
    mode = args.mode or cc.get("mmpbsa_mode") or params.get("mode", "single_point")
    ligand_type = cc.get("ligand_type", "small_molecule")
    peptide_sequence = cc.get("peptide_sequence")
    pose_selection = params.get("pose_selection", "best_score")
    pose_index = params.get("pose_index", 1)
    receptor_ff = params.get("receptor_ff", "ff14SB")
    ligand_ff = params.get("ligand_ff", "gaff2")
    charge_method = params.get("charge_method", "bcc")
    antechamber_timeout = params.get("antechamber_timeout", 300)

    mmpbsa_cfg = params.get("mmpbsa", {})
    mmpbsa_idecomp = mmpbsa_cfg.get("idecomp", 1)
    mmpbsa_igb = mmpbsa_cfg.get("igb", 2)
    mmpbsa_saltcon = mmpbsa_cfg.get("saltcon", 0.15)

    md_params = params.get("md", {})
    if args.production_ns is not None:
        md_params["production_ns"] = args.production_ns

    log_level = args.log_level or params.get("log_level", "INFO")
    logging.getLogger().setLevel(getattr(logging, log_level.upper()))

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================
    setup_log_file(output_dir / "01g_mmpbsa_decomp.log", log_level)

    logger.info("=" * 60)
    logger.info("  REFERENCE_DOCKING - Module 01g: MMPBSA Decomposition")
    logger.info("=" * 60)
    logger.info(f"Campaign:   {campaign_id}")
    logger.info(f"Molecule:   {molecule_name}")
    logger.info(f"Mode:       {mode}")
    logger.info(f"Ligand:     {ligand_type}")
    logger.info(f"Scored:     {scored_mol2.name}")
    logger.info(f"Receptor:   {rec_pdb.name}")
    logger.info(f"Output:     {output_dir}")

    # =========================================================================
    # EXECUTE
    # =========================================================================
    result = run_mmpbsa_decomp(
        scored_mol2=str(scored_mol2),
        receptor_mol2=str(rec_mol2),
        receptor_pdb=str(rec_pdb),
        output_dir=str(output_dir),
        mode=mode,
        pose_selection=pose_selection,
        pose_index=pose_index,
        receptor_ff=receptor_ff,
        ligand_ff=ligand_ff,
        charge_method=charge_method,
        mmpbsa_idecomp=mmpbsa_idecomp,
        mmpbsa_igb=mmpbsa_igb,
        mmpbsa_saltcon=mmpbsa_saltcon,
        md_params=md_params,
        antechamber_timeout=antechamber_timeout,
        ligand_type=ligand_type,
        peptide_sequence=peptide_sequence,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
