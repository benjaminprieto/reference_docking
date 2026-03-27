#!/usr/bin/env python3
"""
00d Binding Site Definition - CLI
====================================
Identifica el binding site y recorta el receptor para grid generation.

Proteinas grandes hacen que sphgen falle. Este modulo recorta
rec_noH.pdb a una esfera alrededor del binding site.

Reads campaign_config.yaml:
    - grids.binding_site.method          -> reference_ligand | residues | coordinates
    - grids.binding_site.reference_mol2  -> ligando de referencia (method A)
    - grids.binding_site.residues        -> lista de residuos (method B)
    - grids.binding_site.center          -> coordenadas (method C)
    - grids.binding_site.radius          -> radio sphere_selector (no trim)
    - receptor.chain                     -> filtro de cadena

Reads module YAML (03_configs/00d_binding_site_definition.yaml):
    - contact_cutoff     -> distancia para contactos ligando-proteina
    - trim_radius        -> radio de recorte del PDB
    - keep_whole_residues

Input:  05_results/{campaign_id}/00b_receptor_preparation/rec_noH.pdb
Output: 05_results/{campaign_id}/00d_binding_site/
    - rec_noH_site.pdb            (PDB recortado para 01b)
    - binding_site_report.json
    - binding_site_summary.txt

Usage:
    python 02_scripts/00d_binding_site_definition.py \
        --config 03_configs/00d_binding_site_definition.yaml \
        --campaigns 04_data/campaigns/example_campaign/campaign_config.yaml

Project: reference_docking
Module: 00d — renumbered from 00e (2026-03-16)
Version: 1.1
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

from reference_docking.m00_preparation.binding_site_definition import run_binding_site_definition

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
        description="00d Binding Site Definition — trim receptor PDB for grid generation",
    )
    parser.add_argument("--config", "-c", type=str, default=None,
                        help="Module config YAML")
    parser.add_argument("--campaigns", type=str, default=None,
                        help="Campaign config YAML")

    # Direct mode overrides
    parser.add_argument("--receptor", type=str, default=None,
                        help="Path to rec_noH.pdb (overrides auto-detection)")
    parser.add_argument("--reference-mol2", type=str, default=None,
                        help="Path to reference ligand mol2")
    parser.add_argument("--residues", nargs="+", default=None,
                        help="Residue IDs, e.g. SER575 HIS335 ASP361 TRP392 GLU529 TRP555 CYS574")
    parser.add_argument("--center", nargs=3, type=float, default=None,
                        help="Binding site center: x y z")
    parser.add_argument("--output", "-o", type=str, default=None)
    parser.add_argument("--trim-radius", type=float, default=None)
    parser.add_argument("--contact-cutoff", type=float, default=None)
    parser.add_argument("--log-level", type=str, default=None,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================

    receptor_noH = None
    output_dir = None
    campaign_id = "direct"
    method = "reference_ligand"
    reference_mol2 = None
    residue_ids = None
    center = None
    chain = None
    contact_cutoff = 5.0
    trim_radius = 25.0
    keep_whole_residues = True
    log_level = "INFO"

    # --- Campaign config ---
    if args.campaigns:
        cc = load_yaml(args.campaigns)
        campaign_dir = Path(args.campaigns).parent
        campaign_id = cc.get("campaign_id", campaign_dir.name)

        # Receptor
        rec_config = cc.get("receptor", {})
        chain = rec_config.get("chain")

        # rec_noH.pdb from 00b
        receptor_noH = str(
            Path("05_results") / campaign_id / "00b_receptor_preparation" / "rec_noH.pdb"
        )

        # Binding site from grids section
        gc = cc.get("grids", {})
        bs = gc.get("binding_site", {})
        method = bs.get("method", "reference_ligand")

        if method == "reference_ligand":
            ref_path = bs.get("reference_mol2")
            if ref_path:
                reference_mol2 = str(
                    Path(ref_path) if Path(ref_path).is_absolute()
                    else campaign_dir / ref_path
                )

        elif method == "residues":
            residue_ids = bs.get("residues", [])

        elif method == "coordinates":
            center = bs.get("center")


        output_dir = str(Path("05_results") / campaign_id / "00d_binding_site")

    # --- Module config ---
    if args.config:
        mc = load_yaml(args.config)
        params = mc.get("parameters", {})
        contact_cutoff = params.get("contact_cutoff", contact_cutoff)
        trim_radius = params.get("trim_radius", trim_radius)
        keep_whole_residues = params.get("keep_whole_residues", keep_whole_residues)
        log_level = params.get("log_level", log_level)

    # --- CLI overrides ---
    if args.receptor:
        receptor_noH = args.receptor
    if args.reference_mol2:
        reference_mol2 = args.reference_mol2
        method = "reference_ligand"
    if args.residues:
        residue_ids = args.residues
        method = "residues"
    if args.center:
        center = args.center
        method = "coordinates"
    if args.output:
        output_dir = args.output
    if args.trim_radius:
        trim_radius = args.trim_radius
    if args.contact_cutoff:
        contact_cutoff = args.contact_cutoff
    if args.log_level:
        log_level = args.log_level

    # --- Validate ---
    if not receptor_noH:
        parser.error("Provide --campaigns or --receptor. Run 00b first.")
    if not output_dir:
        output_dir = "05_results/00d_binding_site"

    if not Path(receptor_noH).exists():
        logger.error(f"Receptor PDB not found: {receptor_noH}")
        logger.error("Run module 00b first.")
        return 1

    if method == "reference_ligand" and (not reference_mol2 or not Path(reference_mol2).exists()):
        logger.error(f"Reference ligand not found: {reference_mol2}")
        logger.error("Set grids.binding_site.reference_mol2 in campaign_config.yaml")
        return 1

    # =========================================================================
    # SETUP LOGGING
    # =========================================================================

    if log_level:
        logging.getLogger().setLevel(getattr(logging, log_level.upper(), logging.INFO))


    log_path = Path(output_dir) / "00d_binding_site.log"
    setup_log_file(log_path, log_level)

    # =========================================================================
    # EXECUTE
    # =========================================================================

    logger.info("=" * 60)
    logger.info("  REFERENCE_DOCKING - Module 00d: Binding Site Definition")
    logger.info("=" * 60)
    logger.info(f"Campaign:       {campaign_id}")
    logger.info(f"Receptor:       {receptor_noH}")
    logger.info(f"Method:         {method}")
    logger.info(f"Trim radius:    {trim_radius} A")
    logger.info(f"Contact cutoff: {contact_cutoff} A")

    result = run_binding_site_definition(
        receptor_noH_pdb=receptor_noH,
        output_dir=output_dir,
        method=method,
        reference_mol2=reference_mol2,
        residue_ids=residue_ids,
        center=center,
        chain=chain,
        contact_cutoff=contact_cutoff,
        trim_radius=trim_radius,
        keep_whole_residues=keep_whole_residues,
    )

    if not result.get("success"):
        logger.error(f"Error: {result.get('error')}")
        return 1

    logger.info("")

    logger.info(f"Next: python 02_scripts/01b_grid_generation.py "
                f"--config 03_configs/01b_grid_generation.yaml "
                f"--campaigns {args.campaigns or '<campaign_config.yaml>'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())