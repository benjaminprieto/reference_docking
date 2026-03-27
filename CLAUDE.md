# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

A modular pipeline (Python) for **DOCK6 rigid re-docking** of crystallographic ligands with footprint-based interaction validation and GB/SA implicit solvation. Produces per-residue energy profiles for comparison against experimental crystal contacts. Written in Spanish (README, comments) by Benjamin Prieto / Environ Bio.

## Commands

```bash
# Setup
conda env create -f environment.yaml
conda activate reference_docking_env
pip install -e ".[dev]"
bash check_dependencies.sh          # verify DOCK6, ChimeraX, OpenBabel, etc.

# Run tests
pytest                               # all tests (uses pyproject.toml config)
pytest tests/test_pipeline.py::TestDock6Runner  # single test class

# Run a single module
python 02_scripts/00a_reference_ligand_preparation.py \
  --config 03_configs/00a_reference_ligand_preparation.yaml \
  --campaigns 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml

# Run full pipeline
bash run_pipeline.sh 04_data/campaigns/UDX_reference_pH63/campaign_config.yaml
```

Every script in `02_scripts/` follows the same CLI pattern: `--config <yaml> --campaign <campaign_config.yaml>`.

## Architecture

### Directory Convention

| Dir | Role |
|---|---|
| `01_src/reference_docking/` | Core library (logic only, no CLI) |
| `02_scripts/` | CLI wrappers (argparse + YAML → core functions) |
| `03_configs/` | One YAML per module (algorithm parameters) |
| `04_data/campaign/<id>/` | Campaign inputs (receptor PDB, grids, campaign_config.yaml) |
| `05_results/<id>/<module>/` | All pipeline outputs (gitignored) |

### Module Numbering = Execution Order

Pipeline: `00a → 00b → 00d → 01b → 01c → 01d → 01f → 01e → 03a → 04b`

- **00a** Reference ligand preparation (extract HETATM → protonate → Gasteiger charges)
- **00b** Receptor preparation (ChimeraX DockPrep → rec_charged.mol2)
- **00d** Binding site definition (spheres from crystal ligand)
- **01b** Grid generation (DMS → spheres → energy grids)
- **01c** DOCK6 docking (grid_score_primary — rigid or flex)
- **01d** Footprint rescore (fps_primary — per-residue vdW+ES in Cartesian space)
- **01f** GB/SA Hawkins rescore (gbsa_hawkins_primary — implicit solvation)
- **01e** Score collection (parse all scored mol2 → Excel)
- **03a** PLIP interaction analysis (crystal complex → interaction JSON)
- **04b** Footprint analysis (per-residue energy consensus)

### DOCK6.13 Scoring Architecture

DOCK6.13 silently ignores all `_secondary` scoring functions during flexible docking. Therefore:
- **01c** uses only `grid_score_primary` (docking step)
- **01d** uses `footprint_similarity_score_primary` (separate rescore)
- **01f** uses `gbsa_hawkins_score_primary` (separate rescore)

Each scoring function runs as an independent rigid rescore pass with `orient_ligand=no`.

### Two-Layer Design

Each module has exactly two files:
1. **Core** (`01_src/.../module.py`): Pure logic, receives paths/dicts, returns results. All the science lives here.
2. **Script** (`02_scripts/module.py`): Thin CLI wrapper that merges YAML config + campaign_config + CLI overrides, then calls the core function.

When modifying behavior, edit the core module. When changing CLI args or config loading, edit the script.

### Campaign Config is the Source of Truth

`campaign_config.yaml` defines everything about a docking experiment: receptor PDB, ligand ID, pH, binding site, grid paths. Module YAMLs in `03_configs/` control only algorithm parameters.

### Key Domain Constraints

- **DOCK6 80-char path limit**: Fortran programs truncate paths. The pipeline uses symlinks automatically — don't break this by hardcoding long absolute paths in DOCK6-related modules.
- **No antechamber for reference ligands**: antechamber recenters molecules to the origin, destroying crystallographic coordinates. Module 00a uses OpenBabel instead.
- **Receptor mol2 must have AMBER ff14SB charges** and Sybyl atom types for DOCK6 compatibility. The `receptor_preparation.py` module handles this via ChimeraX.
- **Reference ligand coordinates** must come from the crystal structure PDB, not regenerated from SMILES, to correctly define the binding site.

## Testing

Tests are smoke tests (imports + basic logic). They don't require DOCK6/ChimeraX installed. No fixtures or mocks for external tools — the pipeline depends on system-installed binaries that can't be meaningfully mocked.

## Note on README

The README is written in Spanish. Maintain Spanish for README updates. Code, docstrings, and this file use English.
