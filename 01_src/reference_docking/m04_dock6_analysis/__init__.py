"""
m04_dock6_analysis - DOCK6 Post-Docking Analysis
===================================================
Residue-based analysis of DOCK6 docking results using footprint scoring.

Pipeline:
    01d footprint_rescore      — Re-score poses with fps_primary (in m01)
    04b footprint_analysis     — Per-residue energy parsing + consensus (reads 01d)
    04b footprint_rescoring    — Rescore with alternate references

Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56
"""
__version__ = "2.0.0"
