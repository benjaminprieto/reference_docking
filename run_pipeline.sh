#!/bin/bash
# =============================================================================
# reference_docking — Full Pipeline
# =============================================================================
# Rigid re-docking of crystallographic ligands with DOCK6, followed by
# footprint-based interaction analysis and PLIP crystal validation.
#
# Usage:
#   bash run_pipeline.sh <campaign_config.yaml>
#
# Example:
#   bash run_pipeline.sh 04_data/campaigns/SD1_reference_pH63/campaign_config.yaml
#
# Pipeline:
#   00b  Receptor preparation (chimera DockPrep → rec_charged.mol2)
#   00d  Binding site definition (spheres from crystal ligand)
#   01b  Grid generation (DMS → spheres → energy grids)
#   01c  DOCK6 rigid re-docking (crystal ligand → scored poses)
#   01d  Footprint rescore (re-score with fps_primary)
#   01e  Score collection (parse → Excel)
#   03a  PLIP interaction analysis (crystal complex → interaction JSON)
#   04b  Footprint analysis (per-residue energy + consensus)
# =============================================================================
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: bash run_pipeline.sh <campaign_config.yaml>"
    exit 1
fi

CAMPAIGN="$1"
CONFIGS="03_configs"
SCRIPTS="02_scripts"

echo "============================================================"
echo "  REFERENCE_DOCKING — Pipeline"
echo "  Campaign: $CAMPAIGN"
echo "============================================================"
echo ""

# --- 00b: Receptor Preparation ---
echo "[00b] Receptor Preparation"
python "$SCRIPTS/00b_receptor_preparation.py" \
    --config "$CONFIGS/00b_receptor_preparation.yaml" \
    --campaign "$CAMPAIGN"

# --- 00d: Binding Site Definition ---
echo "[00d] Binding Site Definition"
python "$SCRIPTS/00d_binding_site_definition.py" \
    --config "$CONFIGS/00d_binding_site_definition.yaml" \
    --campaign "$CAMPAIGN"

# --- 01b: Grid Generation ---
echo "[01b] Grid Generation"
python "$SCRIPTS/01b_grid_generation.py" \
    --config "$CONFIGS/01b_grid_generation.yaml" \
    --campaign "$CAMPAIGN"

# --- 01c: DOCK6 Rigid Re-Docking ---
echo "[01c] DOCK6 Rigid Re-Docking"
python "$SCRIPTS/01c_dock6_run.py" \
    --config "$CONFIGS/01c_dock6_run.yaml" \
    --campaign "$CAMPAIGN"

# --- 01d: Footprint Rescore ---
echo "[01d] Footprint Rescore"
python "$SCRIPTS/01d_footprint_rescore.py" \
    --config "$CONFIGS/01d_footprint_rescore.yaml" \
    --campaign "$CAMPAIGN"

# --- 01e: Score Collection ---
echo "[01e] Score Collection"
python "$SCRIPTS/01e_score_collection.py" \
    --config "$CONFIGS/01e_score_collection.yaml" \
    --campaign "$CAMPAIGN"

# --- 03a: PLIP Interaction Analysis ---
echo "[03a] PLIP Interaction Analysis"
python "$SCRIPTS/03a_plip_interaction_analysis.py" \
    --config "$CONFIGS/03a_plip_interaction_analysis.yaml" \
    --campaign "$CAMPAIGN"

# --- 04b: Footprint Analysis ---
echo "[04b] Footprint Analysis"
python "$SCRIPTS/04b_footprint_analysis.py" \
    --config "$CONFIGS/04b_footprint_analysis.yaml" \
    --campaign "$CAMPAIGN"

echo ""
echo "============================================================"
echo "  PIPELINE COMPLETE"
echo "  Results: 05_results/"
echo "============================================================"
