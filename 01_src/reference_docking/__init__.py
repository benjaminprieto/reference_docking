"""
reference_docking - Crystallographic Reference Docking Pipeline
=================================================================
DOCK6 rigid re-docking of co-crystallized ligands with footprint-based
interaction validation and GB/SA solvation scoring.

Produces per-residue energy profiles for comparison against experimental
crystal contacts. Output feeds into molecular_docking pharmacophore.

Modules:
    m00_preparation      - Ligand preparation (00a), receptor (00b), binding site (00d)
    m01_docking          - DOCK6 engine: grids (01b), docking (01c), footprint (01d),
                           GB/SA (01f), scores (01e)
    m03_crystal_analysis - PLIP interaction analysis (03a)
    m04_dock6_analysis   - Footprint analysis: per-residue energy consensus (04b)
"""
__version__ = "2.0.0"