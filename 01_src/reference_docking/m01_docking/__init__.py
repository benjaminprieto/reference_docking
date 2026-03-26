"""
m01_docking - DOCK6 Engine
============================
DOCK6 docking pipeline for reference complexes.

Modules:
    grid_generation    01b — DMS → spheres → grids
    dock6_runner       01c — dock6 docking (grid_score primary)
    footprint_rescore  01d — re-score with fps_primary (per-residue vdW+ES)
    gbsa_rescore       01f — re-score with gbsa_hawkins_primary (solvation)
    score_collector    01e — parse scores → Excel

DOCK6.13 note: all _secondary score params are ignored during flex docking.
Each scoring method runs as a separate rigid rescore step.
"""