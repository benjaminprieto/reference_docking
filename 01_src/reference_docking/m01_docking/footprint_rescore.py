"""
DOCK6 Footprint Re-scoring — Core Module (01d)
=================================================
Re-scores existing DOCK6 poses with per-residue energy decomposition.

The footprint score calculates VDW + ES interaction energy per receptor
residue in Cartesian space (not grid space). This is a POST-DOCKING step:
it reads scored mol2 from 01c and generates new mol2 with FPS columns.

Based on Rizzo Lab protocol (run.009) and DOCK6 manual §2.11.9:
  - orient_ligand = no  (poses already positioned)
  - grid_score_primary = no
  - footprint_similarity_score_primary = yes
  - fps_score_receptor_filename = receptor mol2
  - fps_score_footprint_reference_mol2_filename = reference ligand

v2.0 (2026-03-25): Added GB/SA Hawkins implicit solvation as secondary
  score. Corrects in-vacuo electrostatic artifacts (e.g., ASP494 repulsion
  caused by unscreened charge-charge interactions without explicit water).
  GB/SA gives a global solvation-corrected score per pose; the per-residue
  footprint remains in-vacuo (DOCK6 limitation).

Input:  01c_dock6_run/{name}/{name}_scored.mol2
Output: 01d_footprint_rescore/{name}/{name}_fps.mol2

Location: 01_src/reference_docking/m01_docking/footprint_rescore.py
Project: reference_docking
Module: 01d (core)
Version: 2.0
"""

import logging
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# DOCK6 FOOTPRINT RE-SCORE TEMPLATE
# =============================================================================
# {gbsa_hawkins_block} is replaced at runtime with either:
#   - "gbsa_hawkins_score_secondary  no"  (default, in-vacuo only)
#   - Full GB/SA Hawkins block          (implicit solvation enabled)
# =============================================================================

FPS_RESCORE_TEMPLATE = """\
conformer_search_type                            rigid
use_internal_energy                              no
ligand_atom_file                                 {ligand_mol2}
limit_max_ligands                                no
skip_molecule                                    no
read_mol_solvation                               no
calculate_rmsd                                   no
use_database_filter                              no
orient_ligand                                    no
bump_filter                                      no
score_molecules                                  yes
contact_score_primary                            no
contact_score_secondary                          no
grid_score_primary                               no
grid_score_secondary                             no
multigrid_score_primary                          no
multigrid_score_secondary                        no
dock3.5_score_primary                            no
dock3.5_score_secondary                          no
continuous_score_primary                         no
continuous_score_secondary                       no
footprint_similarity_score_primary               yes
footprint_similarity_score_secondary             no
fps_score_use_footprint_reference_mol2           yes
fps_score_footprint_reference_mol2_filename      {reference_mol2}
fps_score_foot_compare_type                      Euclidean
fps_score_normalize_foot                         no
fps_score_foot_comp_all_residue                  yes
fps_score_receptor_filename                      {receptor_mol2}
fps_score_vdw_att_exp                            6
fps_score_vdw_rep_exp                            9
fps_score_vdw_rep_rad_scale                      1
fps_score_use_distance_dependent_dielectric      yes
fps_score_dielectric                             4.0
fps_score_vdw_fp_scale                           1
fps_score_es_fp_scale                            1
fps_score_hb_fp_scale                            0
pharmacophore_score_secondary                    no
descriptor_score_secondary                       no
gbsa_zou_score_secondary                         no
{gbsa_hawkins_block}
SASA_score_secondary                             no
amber_score_secondary                            no
minimize_ligand                                  no
atom_model                                       all
vdw_defn_file                                    {vdw_defn_file}
flex_defn_file                                   {flex_defn_file}
flex_drive_file                                  {flex_drive_file}
ligand_outfile_prefix                            {output_prefix}
write_footprints                                 yes
write_hbonds                                     yes
write_orientations                               no
num_scored_conformers                            {num_scored_conformers}
rank_ligands                                     no
"""


# =============================================================================
# GB/SA HAWKINS BLOCK GENERATORS
# =============================================================================

def _gbsa_hawkins_block_off() -> str:
    """Return disabled GB/SA line."""
    return "gbsa_hawkins_score_secondary                     no"


def _gbsa_hawkins_block_on(
        solvent_dielectric: float = 78.5,
        salt_concentration: float = 0.15,
        gb_offset: float = 0.09,
) -> str:
    """
    Return full GB/SA Hawkins block for DOCK6.

    Uses continuous vdW+ES (Cartesian space, consistent with footprint).
    Solvent dielectric 78.5 = water at 25°C.
    Salt concentration 0.15 M = physiological (~150 mM NaCl).

    Reference: Hawkins, Cramer, Truhlar. J Phys Chem 1996.
    DOCK6 manual §2.11.7 (Hawkins GB/SA Score).
    """
    return f"""\
gbsa_hawkins_score_secondary                     yes
gbsa_hawkins_score_rec_filename                  receptor.mol2
gbsa_hawkins_score_solvent_dielectric            {solvent_dielectric}
gbsa_hawkins_score_salt_conc                     {salt_concentration}
gbsa_hawkins_score_gb_offset                     {gb_offset}
gbsa_hawkins_score_cont_vdw_and_es               yes
gbsa_hawkins_score_vdw_att_exp                   6
gbsa_hawkins_score_vdw_rep_exp                   12"""


# =============================================================================
# SYMLINK HELPER (reuse pattern from dock6_runner)
# =============================================================================

def _create_symlink(target: str, link_path: Path) -> bool:
    """Create a symlink, removing any existing link/file."""
    target_path = Path(target).resolve()
    if not target_path.exists():
        logger.debug(f"    Symlink target not found: {target_path}")
        return False
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()
    link_path.symlink_to(target_path)
    return True


def _find_dock6_params(dock6_home: str = "/opt/dock6") -> Dict[str, str]:
    """Find DOCK6 parameter files."""
    params_dir = Path(dock6_home) / "parameters"
    return {
        "vdw_defn_file": str(params_dir / "vdw_AMBER_parm99.defn"),
        "flex_defn_file": str(params_dir / "flex.defn"),
        "flex_drive_file": str(params_dir / "flex_drive.tbl"),
    }


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_footprint_rescore(
        docking_dir: Union[str, Path],
        output_dir: Union[str, Path],
        receptor_mol2: str,
        reference_mol2: str,
        dock6_home: str = "/opt/dock6",
        timeout_per_molecule: int = 300,
        molecule_filter: Optional[List[str]] = None,
        gbsa_hawkins: bool = False,
        solvent_dielectric: float = 78.5,
        salt_concentration: float = 0.15,
        gb_offset: float = 0.09,
) -> Dict[str, Any]:
    """
    Re-score all DOCK6 poses with footprint decomposition.

    Args:
        docking_dir:    01c_dock6_run directory (contains {name}/ subdirs)
        output_dir:     01d_footprint_rescore output directory
        receptor_mol2:  Path to receptor mol2 (from 00b)
        reference_mol2: Path to reference ligand mol2 (crystallographic)
        dock6_home:     DOCK6 installation path
        timeout_per_molecule: Timeout in seconds per molecule
        molecule_filter: Optional list of molecule names to process
        gbsa_hawkins:   Enable GB/SA Hawkins implicit solvation (secondary)
        solvent_dielectric: Solvent dielectric constant (78.5 = water)
        salt_concentration: Salt concentration in M (0.15 = physiological)
        gb_offset:      GB radius offset (default 0.09)

    Returns:
        Dict with n_total, n_ok, n_failed, results
    """
    docking_dir = Path(docking_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate inputs
    if not docking_dir.exists():
        return {"success": False, "error": f"Docking dir not found: {docking_dir}"}
    if not Path(receptor_mol2).exists():
        return {"success": False, "error": f"Receptor mol2 not found: {receptor_mol2}"}
    if not Path(reference_mol2).exists():
        return {"success": False, "error": f"Reference mol2 not found: {reference_mol2}"}

    # Find scored mol2 files from 01c
    mol_dirs = sorted([d for d in docking_dir.iterdir() if d.is_dir()])
    molecules = []
    for d in mol_dirs:
        name = d.name
        scored = d / f"{name}_scored.mol2"
        if scored.exists() and scored.stat().st_size > 0:
            molecules.append((name, scored))

    if molecule_filter:
        molecules = [(n, s) for n, s in molecules if n in molecule_filter]

    if not molecules:
        return {"success": False, "error": "No scored mol2 files found in docking dir"}

    # Find DOCK6 parameter files
    dock6_params = _find_dock6_params(dock6_home)

    # Build GB/SA block
    if gbsa_hawkins:
        gbsa_block = _gbsa_hawkins_block_on(
            solvent_dielectric=solvent_dielectric,
            salt_concentration=salt_concentration,
            gb_offset=gb_offset,
        )
    else:
        gbsa_block = _gbsa_hawkins_block_off()

    logger.info("=" * 60)
    logger.info("  DOCK6 Footprint Re-scoring")
    logger.info("=" * 60)
    logger.info(f"  Docking dir:  {docking_dir}")
    logger.info(f"  Receptor:     {Path(receptor_mol2).name}")
    logger.info(f"  Reference:    {Path(reference_mol2).name}")
    logger.info(f"  Molecules:    {len(molecules)}")
    logger.info(f"  Timeout:      {timeout_per_molecule}s per molecule")
    logger.info(f"  GB/SA Hawkins: {'YES (dielectric={}, salt={}M)'.format(solvent_dielectric, salt_concentration) if gbsa_hawkins else 'no'}")

    results = []
    total_time = 0

    for i, (name, scored_mol2) in enumerate(molecules, 1):
        logger.info(f"  [{i}/{len(molecules)}] {name}")

        mol_out = output_dir / name
        mol_out.mkdir(parents=True, exist_ok=True)

        # Create symlinks (80-char path fix)
        _create_symlink(str(scored_mol2), mol_out / "poses.mol2")
        _create_symlink(receptor_mol2, mol_out / "receptor.mol2")
        _create_symlink(reference_mol2, mol_out / "reference.mol2")

        for key in ["vdw_defn_file", "flex_defn_file", "flex_drive_file"]:
            src = dock6_params[key]
            _create_symlink(src, mol_out / Path(src).name)

        # Count poses in scored mol2
        n_poses = 0
        try:
            text = scored_mol2.read_text()
            n_poses = text.count("@<TRIPOS>MOLECULE")
        except Exception:
            n_poses = 100  # safe fallback

        # Generate dock6_fps.in
        fps_in = FPS_RESCORE_TEMPLATE.format(
            ligand_mol2="poses.mol2",
            receptor_mol2="receptor.mol2",
            reference_mol2="reference.mol2",
            vdw_defn_file=Path(dock6_params["vdw_defn_file"]).name,
            flex_defn_file=Path(dock6_params["flex_defn_file"]).name,
            flex_drive_file=Path(dock6_params["flex_drive_file"]).name,
            output_prefix=f"{name}_fps",
            num_scored_conformers=max(n_poses, 1),
            gbsa_hawkins_block=gbsa_block,
        )

        fps_in_path = mol_out / "dock6_fps.in"
        fps_in_path.write_text(fps_in)

        # Run dock6
        t0 = time.time()
        try:
            proc = subprocess.run(
                ["dock6", "-i", "dock6_fps.in", "-o", "dock6_fps.out"],
                cwd=str(mol_out),
                capture_output=True, text=True,
                timeout=timeout_per_molecule,
            )
            runtime = time.time() - t0
            fps_mol2 = mol_out / f"{name}_fps_scored.mol2"
            success = proc.returncode == 0 and fps_mol2.exists() and fps_mol2.stat().st_size > 0

            results.append({
                "Name": name,
                "status": "OK" if success else "FAILED",
                "fps_mol2": str(fps_mol2) if success else None,
                "runtime_sec": round(runtime, 1),
                "error": proc.stderr[:200] if not success else None,
            })

            total_time += runtime

            if success:
                logger.info(f"    -> OK ({runtime:.1f}s, {n_poses} poses)")
            else:
                logger.warning(f"    -> FAILED ({runtime:.1f}s)")

        except subprocess.TimeoutExpired:
            runtime = time.time() - t0
            total_time += runtime
            results.append({
                "Name": name,
                "status": "TIMEOUT",
                "fps_mol2": None,
                "runtime_sec": round(runtime, 1),
                "error": f"Timeout after {timeout_per_molecule}s",
            })
            logger.warning(f"    -> TIMEOUT ({timeout_per_molecule}s)")

        except Exception as e:
            runtime = time.time() - t0
            total_time += runtime
            results.append({
                "Name": name,
                "status": "FAILED",
                "fps_mol2": None,
                "runtime_sec": round(runtime, 1),
                "error": str(e),
            })
            logger.warning(f"    -> ERROR: {e}")

    # Save status CSV
    df_status = pd.DataFrame(results)
    status_csv = output_dir / "fps_status.csv"
    df_status.to_csv(status_csv, index=False, encoding="utf-8")

    n_ok = sum(1 for r in results if r["status"] == "OK")
    n_fail = len(results) - n_ok

    logger.info("")
    logger.info(f"{'=' * 60}")
    logger.info(f"  FOOTPRINT: {n_ok}/{len(results)} completed ({total_time:.0f}s)")
    if gbsa_hawkins:
        logger.info(f"  GB/SA Hawkins: enabled (scores in mol2 header)")
    logger.info(f"{'=' * 60}")

    return {
        "success": True,
        "n_total": len(results),
        "n_ok": n_ok,
        "n_failed": n_fail,
        "total_runtime_sec": round(total_time, 1),
        "gbsa_hawkins": gbsa_hawkins,
        "results": results,
    }