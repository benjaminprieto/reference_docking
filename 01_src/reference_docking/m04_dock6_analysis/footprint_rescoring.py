"""
Footprint Re-scoring - Core Module (04b rescore)
====================================================
Re-scores existing DOCK6 poses with footprint_similarity_score_primary
to obtain per-residue vdW + ES energy decomposition.

This is fast (seconds per molecule), no re-docking. It evaluates each
docked pose against a reference ligand (e.g., UDX) and the receptor,
writing per-residue FPS_vdw / FPS_es fields into the scored mol2 header.

Depends on m01 execution utilities only:
    - find_dock6_params()   → locates vdw_defn, flex_defn, flex_drive
    - run_dock6_single()    → executes dock6 binary
    - _create_symlink()     → 80-char path workaround

Reference: Balius et al. J Chem Inf Model 2011, 51(8):1942-56
Tutorial:  https://ringo.ams.stonybrook.edu/index.php/2023_DOCK_tutorial_3_with_PDBID_2P16

Location: 01_src/reference_docking/m04_dock6_analysis/footprint_rescoring.py
Project: reference_docking
Module: 04b (DOCK6 analysis — rescore phase)
Version: 1.0 (2026-03-22)
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional, Union

import pandas as pd

from reference_docking.m01_docking.dock6_runner import (
    find_dock6_params,
    run_dock6_single,
    _create_symlink,
)

logger = logging.getLogger(__name__)


# =============================================================================
# DOCK6 FOOTPRINT RE-SCORE TEMPLATE
# =============================================================================
# Rigid re-scoring only (orient_ligand=no, bump_filter=no, minimize=no).
# footprint_similarity_score_primary=yes triggers per-residue decomposition.
# =============================================================================

DOCK6_FOOTPRINT_RESCORE_TEMPLATE = """\
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             {ligand_mol2}
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_primary                           yes
fps_score_use_footprint_reference_mol2                       yes
fps_score_footprint_reference_mol2_filename                  {reference_mol2}
fps_score_receptor_filename                                  {receptor_mol2}
fps_score_foot_compare_type                                  Euclidean
fps_score_normalize_foot                                     no
fps_score_foot_comp_all_residue                              yes
fps_score_vdw_att_exp                                        6
fps_score_vdw_rep_exp                                        9
fps_score_vdw_rep_rad_scale                                  1
fps_score_use_distance_dependent_dielectric                  yes
fps_score_dielectric                                         4.0
fps_score_vdw_fp_scale                                       1
fps_score_es_fp_scale                                        1
fps_score_hb_fp_scale                                        0
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                {vdw_defn_file}
flex_defn_file                                               {flex_defn_file}
flex_drive_file                                              {flex_drive_file}
ligand_outfile_prefix                                        {output_prefix}
write_orientations                                           no
num_scored_conformers                                        {num_scored_conformers}
rank_ligands                                                 no
"""


# =============================================================================
# INPUT FILE GENERATION
# =============================================================================

def generate_footprint_rescore_input(
        ligand_mol2: str,
        reference_mol2: str,
        receptor_mol2: str,
        output_prefix: str,
        output_path: str,
        dock6_params: Optional[Dict[str, str]] = None,
        num_scored_conformers: int = 100,
) -> str:
    """
    Generate a DOCK6 input file for footprint re-scoring.

    This does NOT re-dock — it evaluates existing poses against a reference
    molecule using footprint_similarity_score as primary scoring.

    The output mol2 header contains per-residue vdW + ES decomposition
    for both the docked ligand AND the reference.

    Args:
        ligand_mol2:   Scored poses mol2 (from 01c docking)
        reference_mol2: Reference ligand mol2 (e.g., UDX best pose)
        receptor_mol2:  Receptor mol2 (rec_charged.mol2 from 00b)
        output_prefix:  Prefix for output scored mol2
        output_path:    Path to write the dock6.in file
        dock6_params:   DOCK6 parameter file paths
        num_scored_conformers: Max conformers to score (100 = all poses)

    Returns:
        Path to the generated input file
    """
    if dock6_params is None:
        dock6_params = find_dock6_params()

    template_vars = {
        "ligand_mol2": ligand_mol2,
        "reference_mol2": reference_mol2,
        "receptor_mol2": receptor_mol2,
        "output_prefix": output_prefix,
        "vdw_defn_file": dock6_params.get("vdw_defn_file", "vdw_AMBER_parm99.defn"),
        "flex_defn_file": dock6_params.get("flex_defn_file", "flex.defn"),
        "flex_drive_file": dock6_params.get("flex_drive_file", "flex_drive.tbl"),
        "num_scored_conformers": num_scored_conformers,
    }

    content = DOCK6_FOOTPRINT_RESCORE_TEMPLATE.format(**template_vars)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(content)

    return output_path


# =============================================================================
# FOOTPRINT RE-SCORING ORCHESTRATOR
# =============================================================================

def run_footprint_rescoring(
        docking_dir: Union[str, Path],
        output_dir: Union[str, Path],
        reference_mol2: Union[str, Path],
        receptor_mol2: Union[str, Path],
        dock6_home: Optional[str] = None,
        timeout: int = 120,
        num_scored_conformers: int = 100,
        dry_run: bool = False,
) -> Dict[str, Any]:
    """
    Re-score all docked poses with DOCK6 footprint similarity.

    Scans docking_dir for {name}/{name}_scored.mol2 files and runs
    rigid re-scoring with footprint_similarity_score_primary=yes
    against the reference molecule.

    This is fast — seconds per molecule. No re-docking.

    Args:
        docking_dir:   Path to 01c_dock6_run output
        output_dir:    Path for footprint re-scoring output
        reference_mol2: Reference ligand (e.g., UDX best pose from 01d)
        receptor_mol2:  Receptor mol2 (rec_charged.mol2 from 00b)
        dock6_home:    DOCK6 installation path (for param files)
        timeout:       Timeout per molecule in seconds
        num_scored_conformers: Max conformers to score
        dry_run:       Generate inputs only

    Returns:
        Dict with: success, n_molecules, n_rescored, results list
    """
    docking_dir = Path(docking_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_mol2 = Path(reference_mol2).resolve()
    receptor_mol2 = Path(receptor_mol2).resolve()

    if not reference_mol2.exists():
        return {"success": False, "error": f"Reference mol2 not found: {reference_mol2}"}
    if not receptor_mol2.exists():
        return {"success": False, "error": f"Receptor mol2 not found: {receptor_mol2}"}

    logger.info("=" * 60)
    logger.info("  DOCK6 Footprint Re-scoring v1.0")
    logger.info("=" * 60)
    logger.info(f"  Docking dir:  {docking_dir}")
    logger.info(f"  Reference:    {reference_mol2.name}")
    logger.info(f"  Receptor:     {receptor_mol2.name}")
    logger.info(f"  Output:       {output_dir}")

    # --- Find DOCK6 parameter files ---
    # find_dock6_params() discovers via env vars ($DOCK_HOME, $DOCK6_HOME)
    # If dock6_home is explicitly provided, set it before discovery.
    if dock6_home:
        import os
        os.environ["DOCK6_HOME"] = str(dock6_home)
    dock6_params = find_dock6_params()

    # --- Find scored mol2 files ---
    mol_dirs = sorted([
        d for d in docking_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    ])

    scored_files = []
    for d in mol_dirs:
        scored = d / f"{d.name}_scored.mol2"
        if scored.exists() and scored.stat().st_size > 0:
            scored_files.append((d.name, scored))

    logger.info(f"  Molecules found: {len(scored_files)}")

    if not scored_files:
        return {"success": False, "error": "No scored mol2 files found"}

    # --- Re-score each molecule ---
    results = []
    total_time = 0

    for idx, (name, scored_mol2) in enumerate(scored_files):
        mol_out = output_dir / name
        mol_out.mkdir(parents=True, exist_ok=True)

        logger.info(f"  [{idx + 1}/{len(scored_files)}] {name}")

        # --- Symlinks for 80-char path compliance ---
        lig_link = mol_out / "scored_poses.mol2"
        _create_symlink(str(scored_mol2), lig_link)

        ref_link = mol_out / "reference.mol2"
        _create_symlink(str(reference_mol2), ref_link)

        rec_link = mol_out / "receptor.mol2"
        _create_symlink(str(receptor_mol2), rec_link)

        # Symlink parameter files
        short_params = {}
        for key in ["vdw_defn_file", "flex_defn_file", "flex_drive_file"]:
            src = dock6_params.get(key, "")
            if src and (("/" in src) or ("\\" in src)):
                link_name = Path(src).name
                _create_symlink(src, mol_out / link_name)
                short_params[key] = link_name
            else:
                short_params[key] = src

        # --- Generate footprint rescore input ---
        dock_input = "dock6_fps.in"
        dock_output = "dock6_fps.out"
        output_prefix = f"{name}_footprint"
        fps_scored_mol2 = str(mol_out / f"{name}_footprint_scored.mol2")

        generate_footprint_rescore_input(
            ligand_mol2="scored_poses.mol2",
            reference_mol2="reference.mol2",
            receptor_mol2="receptor.mol2",
            output_prefix=output_prefix,
            output_path=str(mol_out / dock_input),
            dock6_params=short_params,
            num_scored_conformers=num_scored_conformers,
        )

        if dry_run:
            results.append({
                "Name": name,
                "status": "DRY_RUN",
                "fps_scored_mol2": None,
                "runtime_sec": 0,
                "error": None,
            })
            logger.info(f"    -> DRY_RUN (input generated)")
            continue

        # --- Execute DOCK6 re-scoring with cwd=mol_out ---
        run_result = run_dock6_single(
            dock_input, dock_output,
            timeout=timeout,
            cwd=str(mol_out),
        )

        has_output = Path(fps_scored_mol2).exists() and Path(fps_scored_mol2).stat().st_size > 0

        status = "OK" if run_result["success"] and has_output else "FAILED"
        error = None
        if not run_result["success"]:
            error = run_result.get("error", "Unknown error")
        elif not has_output:
            error = "dock6 succeeded but no footprint scored mol2 produced"

        results.append({
            "Name": name,
            "status": status,
            "fps_scored_mol2": fps_scored_mol2 if has_output else None,
            "runtime_sec": run_result["runtime_sec"],
            "error": error,
        })

        total_time += run_result["runtime_sec"]

        if status == "OK":
            logger.info(f"    -> OK ({run_result['runtime_sec']:.1f}s)")
        else:
            logger.warning(f"    -> FAILED: {error}")

    # --- Summary ---
    n_ok = sum(1 for r in results if r["status"] == "OK")
    n_fail = sum(1 for r in results if r["status"] == "FAILED")
    n_dry = sum(1 for r in results if r["status"] == "DRY_RUN")

    # Save status CSV
    df_status = pd.DataFrame(results)
    status_csv = output_dir / "footprint_rescore_status.csv"
    df_status.to_csv(status_csv, index=False, encoding="utf-8")

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  FOOTPRINT RE-SCORING: {n_ok}/{len(results)} completed "
                f"({total_time:.0f}s total)")
    if n_fail > 0:
        logger.info(f"  FAILED: {n_fail}")
    logger.info(f"{'=' * 60}")

    return {
        "success": True,
        "n_total": len(results),
        "n_ok": n_ok,
        "n_failed": n_fail,
        "n_dry_run": n_dry,
        "total_runtime_sec": round(total_time, 1),
        "results": results,
        "status_csv": str(status_csv),
        "output_dir": str(output_dir),
    }
