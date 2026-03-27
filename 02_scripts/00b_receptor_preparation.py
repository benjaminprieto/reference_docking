#!/usr/bin/env python3
"""
00b Receptor Preparation - CLI
================================
Prepara el receptor para DOCK6: limpieza + protonacion + mol2 con cargas.

OPCIONAL: Se salta si receptor.protonation.enabled=false.
Si receptor.prepared_mol2 apunta a un archivo valido, se usa directo.

Strategies (v2.0 — ChimeraX backbone):
  pdb2pqr  -> PDB2PQR+PROPKA pKa prediction + ChimeraX mol2 generation
              (RECOMENDADA para docking a pH no-fisiologico, e.g. Golgi 6.3)
  chimerax -> ChimeraX DockPrep directo (AddH + addcharge)
              (buena para pH fisiologico ~7.2)
  obabel   -> OpenBabel simple, Gasteiger charges (FALLBACK — NO RECOMENDADA)

Reads campaign_config.yaml:
    - receptor.pdb                       -> PDB del receptor
    - receptor.protonation.enabled       -> true/false
    - receptor.protonation.tool          -> pdb2pqr | chimerax | obabel
    - receptor.protonation.force_field   -> AMBER | CHARMM | PARSE
    - receptor.prepared_mol2             -> mol2 pre-preparado (skip)
    - receptor.chain                     -> chain(s) a conservar
    - receptor.remove_water              -> true/false
    - receptor.remove_hetatm             -> true/false
    - docking_ph                         -> pH de protonacion

Output: 05_results/{campaign_id}/00b_receptor_preparation/
    - rec_charged.mol2       (DOCK6-ready: Sybyl types + AMBER ff14SB charges)
    - rec_noH.pdb            (para DMS surface en 01b)
    - receptor_clean.pdb
    - chimerax.log           (log de ChimeraX, si aplica)
    - protonation_report.json
    - protonation_summary.txt
    - 00b_receptor_preparation.log

Usage:
    # Recommended (PDB2PQR + ChimeraX):
    python 02_scripts/00b_receptor_preparation.py --config 03_configs/00b_receptor_preparation.yaml --campaigns 04_data/campaigns/example_campaign/campaign_config.yaml

    # Override pH:
    python 02_scripts/00b_receptor_preparation.py \\
        --config 03_configs/00b_receptor_preparation.yaml \\
        --campaigns 04_data/campaigns/example_campaign/campaign_config.yaml \\
        --ph 6.3

    # ChimeraX only (no PDB2PQR):
    python 02_scripts/00b_receptor_preparation.py \\
        --config 03_configs/00b_receptor_preparation.yaml \\
        --campaigns 04_data/campaigns/example_campaign/campaign_config.yaml \\
        --tool chimerax

Project: reference_docking
Module: 00b
Version: 2.0 — ChimeraX rewrite (2026-03-13)
"""

import argparse
import logging
import shutil
import sys
import yaml
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
)

sys.path.insert(0, str(Path(__file__).parent.parent / "01_src"))

from reference_docking.m00_preparation.receptor_preparation import (
    run_receptor_preparation,
    validate_prepared_mol2,
)

logger = logging.getLogger(__name__)


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_log_file(log_path: Path, log_level: str = "INFO"):
    """Add file handler to root logger."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), encoding="utf-8")
    fh.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s"))
    logging.getLogger().addHandler(fh)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare receptor for DOCK6 (clean + protonate + mol2)",
    )
    # Modes
    parser.add_argument("--config", "-c", type=str, default=None,
                        help="Module config YAML")
    parser.add_argument("--campaigns", type=str, default=None,
                        help="Campaign config YAML")

    # Direct mode
    parser.add_argument("--receptor", "-r", type=str, default=None,
                        help="Path to receptor PDB (overrides campaigns)")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Output directory")

    # Parameter overrides
    parser.add_argument("--ph", type=float, default=None,
                        help="Docking pH (overrides campaigns)")
    parser.add_argument("--tool", type=str, default=None,
                        choices=["pdb2pqr", "chimerax", "obabel"],
                        help="Protonation tool (overrides campaigns)")
    parser.add_argument("--force-field", type=str, default=None,
                        choices=["AMBER", "CHARMM", "PARSE"],
                        help="Force field for PDB2PQR")
    parser.add_argument("--chain", type=str, default=None,
                        help="Chain ID to keep (e.g., 'A')")
    parser.add_argument("--keep-water", action="store_true",
                        help="Don't remove water molecules")
    parser.add_argument("--keep-hetatm", action="store_true",
                        help="Don't remove HETATM records")
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    receptor_pdb = None
    prepared_mol2 = None
    output_dir = None
    campaign_id = "direct"
    docking_ph = 7.2
    protonation_enabled = True
    protonation_tool = "pdb2pqr"
    force_field = "AMBER"
    chain = None
    remove_water = True
    remove_hetatm = True
    remove_alt_conformations = True
    log_level = "INFO"

    # --- Campaign config ---
    if args.campaigns:
        cc = load_yaml(args.campaigns)
        campaign_dir = Path(args.campaigns).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)
        docking_ph = cc.get("docking_ph", docking_ph)

        rc = cc.get("receptor", {})
        receptor_pdb = str(campaign_dir / rc.get("pdb", "receptor/receptor.pdb"))
        chain = rc.get("chain")
        remove_water = rc.get("remove_water", remove_water)
        remove_hetatm = rc.get("remove_hetatm", remove_hetatm)

        # Protonation settings
        pc = rc.get("protonation", {})
        protonation_enabled = pc.get("enabled", False)
        protonation_tool = pc.get("tool", protonation_tool)
        force_field = pc.get("force_field", force_field)

        # Pre-prepared mol2
        prepared_mol2_rel = rc.get("prepared_mol2")
        if prepared_mol2_rel:
            prepared_mol2 = str(
                Path(prepared_mol2_rel) if Path(prepared_mol2_rel).is_absolute()
                else campaign_dir / prepared_mol2_rel
            )

        output_dir = str(Path("05_results") / campaign_id / "00b_receptor_preparation")

    # --- Module config ---
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        if params.get("protonation_tool"):
            protonation_tool = params["protonation_tool"]
        force_field = params.get("force_field", force_field)
        remove_alt_conformations = params.get("remove_alt_conformations", True)
        log_level = params.get("log_level", log_level)

    # --- CLI overrides ---
    if args.receptor:
        receptor_pdb = args.receptor
    if args.output:
        output_dir = args.output
    if args.ph is not None:
        docking_ph = args.ph
    if args.tool:
        protonation_tool = args.tool
        protonation_enabled = True  # If tool is specified, enable protonation
    if args.force_field:
        force_field = args.force_field
    if args.chain:
        chain = args.chain
    if args.keep_water:
        remove_water = False
    if args.keep_hetatm:
        remove_hetatm = False
    if args.log_level:
        log_level = args.log_level

    if not output_dir:
        output_dir = "05_results/00b_receptor_preparation"

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))

    log_path = Path(output_dir) / "00b_receptor_preparation.log"
    setup_log_file(log_path, log_level)

    # =========================================================================
    # GATE: SKIP IF PROTONATION DISABLED
    # =========================================================================

    logger.info("=" * 60)
    logger.info("  REFERENCE_DOCKING - Module 00b: Receptor Preparation")
    logger.info("  Version 2.0 (ChimeraX backbone)")
    logger.info("=" * 60)
    logger.info(f"Campaign:     {campaign_id}")

    if not protonation_enabled:
        if prepared_mol2 and Path(prepared_mol2).exists():
            logger.info(f"Protonation DISABLED.")
            logger.info(f"Using pre-prepared mol2: {prepared_mol2}")
            validation = validate_prepared_mol2(prepared_mol2)
            logger.info(f"Validation: atoms={validation['n_atoms']}, "
                        f"residues={validation['n_residues']}, "
                        f"sybyl={validation['has_sybyl_types']}, "
                        f"valid={validation['valid']}")

            # Copy to output dir for downstream consistency
            out_dir = Path(output_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            dest = out_dir / "rec_charged.mol2"
            shutil.copy2(prepared_mol2, dest)
            logger.info(f"Copied to: {dest}")

            if not validation["valid"]:
                logger.warning("WARNING: mol2 validation issues detected. "
                              "Check your prepared mol2 file.")
            return 0

        else:
            logger.info("Protonation DISABLED and no prepared_mol2 provided.")
            logger.info("Options:")
            logger.info("  1. Set receptor.protonation.enabled: true")
            logger.info("  2. Set receptor.prepared_mol2: path/to/rec_charged.mol2")
            logger.info("SKIPPED")
            return 0

    # =========================================================================
    # EXECUTE
    # =========================================================================

    if not receptor_pdb or not Path(receptor_pdb).exists():
        logger.error(f"Receptor PDB not found: {receptor_pdb}")
        return 1

    result = run_receptor_preparation(
        receptor_pdb=receptor_pdb,
        output_dir=output_dir,
        docking_ph=docking_ph,
        protonation_tool=protonation_tool,
        force_field=force_field,
        chain=chain,
        remove_water=remove_water,
        remove_hetatm=remove_hetatm,
        remove_alt_conformations=remove_alt_conformations,
    )

    if not result.get("success"):
        logger.error(f"Receptor preparation failed: {result.get('error')}")
        return 1

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  rec_charged.mol2: {result['rec_charged_mol2']}")
    logger.info(f"  rec_noH.pdb:      {result['rec_noH_pdb']}")
    logger.info(f"{'=' * 60}")
    logger.info(f"Next: python 02_scripts/00d_binding_site_definition.py "
                f"--config 03_configs/00d_binding_site_definition.yaml "
                f"--campaigns {args.campaigns or '<campaign_config.yaml>'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())