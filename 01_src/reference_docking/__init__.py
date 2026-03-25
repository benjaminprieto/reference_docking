"""
reference_docking - Crystallographic Reference Docking Pipeline
=================================================================
DOCK6 rigid re-docking of co-crystallized ligands with footprint-based
interaction validation. Produces per-residue energy profiles for
comparison against experimental crystal contacts.

Modules:
    m00_preparation      - Prepare receptor, define binding site
    m01_docking          - DOCK6 engine (grids, rigid docking, footprint rescore, scores)
    m03_crystal_analysis - PLIP interaction analysis of crystal complexes
    m04_dock6_analysis   - Footprint analysis (per-residue energy decomposition)
"""
__version__ = "1.0.0"
