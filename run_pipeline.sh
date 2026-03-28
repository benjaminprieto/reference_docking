#!/bin/bash
# =============================================================================
# run_pipeline.sh — Reference Docking Pipeline
# =============================================================================
#
# Usage:
#   bash run_pipeline.sh <campaign_config.yaml> [start] [stop]
#
# Examples:
#   bash run_pipeline.sh 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml
#   bash run_pipeline.sh 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml 01c
#   bash run_pipeline.sh 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml 00a 00d
#   bash run_pipeline.sh 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml 01g 01g
#
# Pipeline:
#   00a  Ligand preparation    (crystal mol2 with validated coords)
#   00b  Receptor preparation  (ChimeraX DockPrep → rec_charged.mol2)
#   00d  Binding site definition (trimmed receptor from crystal ligand)
#   01b  Grid generation       (DMS → spheres → tight box 1.0Å → grids)
#   01c  DOCK6 docking         (grid_score primary, rigid/flex)
#   01d  Footprint rescore     (fps_primary, per-residue vdW+ES)
#   01f  GB/SA rescore         (gbsa_hawkins_primary, solvation)
#   01e  Score collection      (parse all scores → Excel/CSV)
#   01g  MMPBSA decomposition  (AMBER MMPBSA per-residue: vdW+ES+GB+SA)
#   01h  MMPBSA analysis       (parse decomp + compare with footprint)
#   03a  PLIP analysis         (crystal complex → interaction JSON)
#   04b  Footprint analysis    (per-residue consensus)
#   06a  Pharmit pharmacophore (footprint+PLIP → Pharmit queries)
#   06b  Pharmit zone selector (map features → binding site zones)
# =============================================================================
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: bash run_pipeline.sh <campaign_config.yaml> [start] [stop]"
    echo ""
    echo "Modules: 00a 00b 00d 01b 01c 01d 01f 01e 01g 01h 03a 04b 06a 06b"
    echo ""
    echo "Examples:"
    echo "  bash run_pipeline.sh config.yaml            # full pipeline"
    echo "  bash run_pipeline.sh config.yaml 01c        # start from docking"
    echo "  bash run_pipeline.sh config.yaml 00a 00d    # 00a through 00d"
    echo "  bash run_pipeline.sh config.yaml 01g 01g    # only 01g"
    exit 1
fi

CAMPAIGN="$1"
START="${2:-00a}"
STOP="${3:-}"
CONFIGS="03_configs"
SCRIPTS="02_scripts"

if [ ! -f "$CAMPAIGN" ]; then
    echo "ERROR: Campaign config not found: $CAMPAIGN"
    exit 1
fi

# Module order
MODULES=(00a 00b 00d 01b 01c 01d 01f 01e 01g 01h 03a 04b 06a 06b)

# Find start index
START_IDX=0
for i in "${!MODULES[@]}"; do
    if [ "${MODULES[$i]}" = "$START" ]; then
        START_IDX=$i
        break
    fi
done

# Find stop index (default: last module)
STOP_IDX=$(( ${#MODULES[@]} - 1 ))
if [ -n "$STOP" ]; then
    for i in "${!MODULES[@]}"; do
        if [ "${MODULES[$i]}" = "$STOP" ]; then
            STOP_IDX=$i
            break
        fi
    done
fi

echo "============================================================"
echo "  REFERENCE_DOCKING — Pipeline v2.1"
echo "  Campaign: $CAMPAIGN"
echo "  Start:    $START"
echo "  Stop:     ${STOP:-end}"
echo "  Started:  $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo ""

run_module() {
    local mod="$1"
    local script="$2"
    local config="$3"
    local label="$4"

    echo "[$mod] $label — $(date '+%H:%M:%S')"
    python "$SCRIPTS/$script" \
        --config "$CONFIGS/$config" \
        --campaigns "$CAMPAIGN"
    echo ""
}

for i in "${!MODULES[@]}"; do
    [ "$i" -lt "$START_IDX" ] && continue
    [ "$i" -gt "$STOP_IDX" ] && break

    case "${MODULES[$i]}" in
        00a) run_module 00a 00a_ligand_preparation.py 00a_ligand_preparation.yaml "Ligand Preparation" ;;
        00b) run_module 00b 00b_receptor_preparation.py 00b_receptor_preparation.yaml "Receptor Preparation" ;;
        00d) run_module 00d 00d_binding_site_definition.py 00d_binding_site_definition.yaml "Binding Site Definition" ;;
        01b) run_module 01b 01b_grid_generation.py 01b_grid_generation.yaml "Grid Generation" ;;
        01c) run_module 01c 01c_dock6_run.py 01c_dock6_run.yaml "DOCK6 Docking" ;;
        01d) run_module 01d 01d_footprint_rescore.py 01d_footprint_rescore.yaml "Footprint Rescore" ;;
        01f) run_module 01f 01f_gbsa_rescore.py 01f_gbsa_rescore.yaml "GB/SA Hawkins Rescore" ;;
        01e) run_module 01e 01e_score_collection.py 01e_score_collection.yaml "Score Collection" ;;
        01g) run_module 01g 01g_mmpbsa_decomp.py 01g_mmpbsa_decomp.yaml "MMPBSA Decomposition" ;;
        01h) run_module 01h 01h_mmpbsa_analysis.py 01h_mmpbsa_analysis.yaml "MMPBSA Analysis" ;;
        03a) run_module 03a 03a_plip_interaction_analysis.py 03a_plip_interaction_analysis.yaml "PLIP Analysis" ;;
        04b) run_module 04b 04b_footprint_analysis.py 04b_footprint_analysis.yaml "Footprint Analysis" ;;
        06a) run_module 06a 06a_pharmit_pharmacophore.py 06a_pharmit_pharmacophore.yaml "Pharmit Pharmacophore" ;;
        06b) run_module 06b 06b_pharmit_zone_selector.py 06b_pharmit_zone_selector.yaml "Pharmit Zone Selector" ;;
    esac
done

echo "============================================================"
echo "  PIPELINE COMPLETE"
echo "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Results:  05_results/"
echo "============================================================"N