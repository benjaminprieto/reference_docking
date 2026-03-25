# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

A modular molecular docking pipeline (Python) with two engines: **DOCK6** and **GNINA** (Vina+CNN). Produces outputs compatible with **dock2profile**. Written in Spanish (README, comments) by Benjamin Prieto / Environ Bio.

## Commands

```bash
# Setup
conda env create -f environment.yaml
conda activate reference_docking_env
pip install -e ".[dev]"
bash check_dependencies.sh          # verify DOCK6, ChimeraX, RDKit, etc.

# Run tests
pytest                               # all tests (uses pyproject.toml config)
pytest tests/test_pipeline.py::TestGridGeneration  # single test class
pytest -k "test_parse_vina_output"   # single test by name

# Run a single module
python 02_scripts/00a_molecule_parser.py \
  --config 03_configs/00a_molecule_parser.yaml \
  --campaign 04_data/campaigns/<campaign_id>/campaign_config.yaml

# Run full pipeline
bash run_pipeline.sh <campaign_id> [dock6|gnina|both]
```

Every script in `02_scripts/` follows the same CLI pattern: `--config <yaml> --campaign <campaign_config.yaml>`.

## Architecture

### Directory Convention

| Dir | Role |
|---|---|
| `01_src/reference_docking/` | Core library (logic only, no CLI) |
| `02_scripts/` | CLI wrappers (argparse + YAML → core functions) |
| `03_configs/` | One YAML per module (algorithm parameters) |
| `04_data/campaigns/<id>/` | Campaign inputs (receptor, molecules, grids, campaign_config.yaml) |
| `05_results/<id>/<module>/` | All pipeline outputs (gitignored) |

### Module Numbering = Execution Order

Modules are numbered `00a`–`02c`. The prefix encodes dependencies:
- **00x** (shared preparation): must run sequentially, pipeline stops on failure
- **01x** (DOCK6 engine): depends on 00x outputs
- **02x** (GNINA engine): depends on 00x outputs, independent of 01x

### Two-Layer Design

Each module has exactly two files:
1. **Core** (`01_src/.../module.py`): Pure logic, receives paths/dicts, returns results. All the science lives here.
2. **Script** (`02_scripts/module.py`): Thin CLI wrapper that merges YAML config + campaign_config + CLI overrides, then calls the core function.

When modifying behavior, edit the core module. When changing CLI args or config loading, edit the script.

### Campaign Config is the Source of Truth

`campaign_config.yaml` defines everything about a docking experiment: receptor, molecules, pH, binding site, grid paths. Module YAMLs in `03_configs/` control only algorithm parameters (exhaustiveness, charge method, conformer strategy, etc.).

### Key Domain Constraints

- **DOCK6 80-char path limit**: Fortran programs truncate paths. The pipeline uses symlinks automatically — don't break this by hardcoding long absolute paths in DOCK6-related modules.
- **SYBYL atom types are mandatory** for DOCK6 flexible docking. GAFF types cause silent fallback to rigid docking. `antechamber_preparation.py` enforces this.
- **Receptor mol2 must have AMBER ff14SB charges** and Sybyl atom types for DOCK6 compatibility. The `receptor_preparation.py` module handles this via ChimeraX.
- **Reference ligand coordinates** must come from the crystal structure PDB, not regenerated from SMILES, to correctly define the binding site.

## Testing

Tests are smoke tests (imports + basic logic). They don't require DOCK6/ChimeraX/GNINA installed. No fixtures or mocks for external tools — the pipeline depends on system-installed binaries that can't be meaningfully mocked.

## Note on README

The README is written in Spanish. Maintain Spanish for README updates. Code, docstrings, and this file use English.
