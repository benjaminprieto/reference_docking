"""
GB/SA Hawkins Re-scoring - Core Module (01f)
===============================================
Re-scores existing DOCK6 poses with GB/SA Hawkins implicit solvation
as PRIMARY scoring function.

DOCK6.13 ignores all _secondary score parameters during flex docking.
This module runs GB/SA as a separate rigid rescore step with
gbsa_hawkins_score_primary=yes.

GB/SA calculates:
  - MM (vdW + ES in Cartesian space, unscaled)
  - deltaGBSA = GBSA_complex - (GBSA_receptor + GBSA_ligand)
  - Total = MM + deltaGBSA

This corrects in-vacuo electrostatic artifacts where charged residues
(e.g., ASP494 COO⁻ near PO₄²⁻) show artificial repulsion without
solvent screening.

Reference:
  Hawkins, Cramer, Truhlar. J Phys Chem 1996.
  Srinivasan et al. J Am Chem Soc 1998 (MM-GBSA).
  DOCK6 manual §2.11.7 (Hawkins GB/SA Score).

Input:  01c_dock6_run/{name}/{name}_scored.mol2
Output: 01f_gbsa_rescore/{name}/{name}_gbsa_scored.mol2

Location: 01_src/reference_docking/m01_docking/gbsa_rescore.py
Project: reference_docking
Module: 01f (core)
Version: 1.0 (2026-03-25)
"""

import logging
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# DOCK6 GB/SA RESCORE TEMPLATE
# =============================================================================
# gbsa_hawkins_score_primary = yes (NOT secondary!)
# DOCK6.13 ignores secondary params. Primary works for rigid rescore.
#
# orient_ligand = no (poses already positioned from 01c)
# minimize_ligand configurable (recommended by DOCK6 docs for GB/SA)
# =============================================================================

GBSA_RESCORE_TEMPLATE = """\
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
grid_score_primary                               no
multigrid_score_primary                          no
dock3.5_score_primary                            no
continuous_score_primary                         no
footprint_similarity_score_primary               no
pharmacophore_score_primary                      no
descriptor_score_primary                         no
gbsa_zou_score_primary                           no
gbsa_hawkins_score_primary                       yes
gbsa_hawkins_score_rec_filename                  {receptor_mol2}
gbsa_hawkins_score_solvent_dielectric            {solvent_dielectric}
gbsa_hawkins_score_salt_conc                     {salt_concentration}
gbsa_hawkins_score_gb_offset                     {gb_offset}
gbsa_hawkins_score_cont_vdw_and_es               yes
gbsa_hawkins_score_vdw_att_exp                   6
gbsa_hawkins_score_vdw_rep_exp                   12
SASA_score_primary                               no
amber_score_primary                              no
minimize_ligand                                  {minimize}
{minimize_params}
atom_model                                       all
vdw_defn_file                                    {vdw_defn_file}
flex_defn_file                                   {flex_defn_file}
flex_drive_file                                  {flex_drive_file}
ligand_outfile_prefix                            {output_prefix}
write_orientations                               no
num_scored_conformers                            {num_scored_conformers}
rank_ligands                                     no
"""

# Minimization params (only included when minimize=yes)
MINIMIZE_PARAMS = """\
simplex_max_iterations                           1000
simplex_max_cycles                               1
simplex_score_converge                           0.1
simplex_cycle_converge                           1.0
simplex_trans_step                               1.0
simplex_rot_step                                 0.1
simplex_tors_step                                10.0
simplex_random_seed                              0
simplex_restraint_min                            no"""


# =============================================================================
# HELPERS
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

def run_gbsa_rescore(
        docking_dir: Union[str, Path],
        output_dir: Union[str, Path],
        receptor_mol2: str,
        dock6_home: str = "/opt/dock6",
        timeout_per_molecule: int = 600,
        molecule_filter: Optional[List[str]] = None,
        minimize: bool = False,
        solvent_dielectric: float = 78.5,
        salt_concentration: float = 0.15,
        gb_offset: float = 0.09,
) -> Dict[str, Any]:
    """
    Re-score all DOCK6 poses with GB/SA Hawkins implicit solvation.

    Args:
        docking_dir:    01c_dock6_run directory
        output_dir:     01f_gbsa_rescore output directory
        receptor_mol2:  Path to receptor mol2 (from 00b)
        dock6_home:     DOCK6 installation path
        timeout_per_molecule: Timeout in seconds
        molecule_filter: Optional molecule name filter
        minimize:       Run simplex minimization before scoring
        solvent_dielectric: Solvent dielectric (78.5 = water)
        salt_concentration: Salt concentration in M (0.15 = physiological)
        gb_offset:      GB radius offset

    Returns:
        Dict with n_total, n_ok, n_failed, results
    """
    docking_dir = Path(docking_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate
    if not docking_dir.exists():
        return {"success": False, "error": f"Docking dir not found: {docking_dir}"}
    if not Path(receptor_mol2).exists():
        return {"success": False, "error": f"Receptor mol2 not found: {receptor_mol2}"}

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
        return {"success": False, "error": "No scored mol2 files found"}

    dock6_params = _find_dock6_params(dock6_home)

    logger.info("=" * 60)
    logger.info("  GB/SA Hawkins Re-scoring (Primary)")
    logger.info("=" * 60)
    logger.info(f"  Docking dir:  {docking_dir}")
    logger.info(f"  Receptor:     {Path(receptor_mol2).name}")
    logger.info(f"  Molecules:    {len(molecules)}")
    logger.info(f"  Minimize:     {'yes' if minimize else 'no'}")
    logger.info(f"  Dielectric:   {solvent_dielectric}")
    logger.info(f"  Salt:         {salt_concentration} M")

    results = []
    total_time = 0

    for i, (name, scored_mol2) in enumerate(molecules, 1):
        logger.info(f"  [{i}/{len(molecules)}] {name}")

        mol_out = output_dir / name
        mol_out.mkdir(parents=True, exist_ok=True)

        # Create symlinks
        _create_symlink(str(scored_mol2), mol_out / "poses.mol2")
        _create_symlink(receptor_mol2, mol_out / "receptor.mol2")

        for key in ["vdw_defn_file", "flex_defn_file", "flex_drive_file"]:
            src = dock6_params[key]
            _create_symlink(src, mol_out / Path(src).name)

        # Count poses
        n_poses = 0
        try:
            text = scored_mol2.read_text()
            n_poses = text.count("@<TRIPOS>MOLECULE")
        except Exception:
            n_poses = 100

        # Generate dock6_gbsa.in
        minimize_str = "yes" if minimize else "no"
        minimize_params = MINIMIZE_PARAMS if minimize else ""

        gbsa_in = GBSA_RESCORE_TEMPLATE.format(
            ligand_mol2="poses.mol2",
            receptor_mol2="receptor.mol2",
            vdw_defn_file=Path(dock6_params["vdw_defn_file"]).name,
            flex_defn_file=Path(dock6_params["flex_defn_file"]).name,
            flex_drive_file=Path(dock6_params["flex_drive_file"]).name,
            output_prefix=f"{name}_gbsa",
            num_scored_conformers=max(n_poses, 1),
            minimize=minimize_str,
            minimize_params=minimize_params,
            solvent_dielectric=solvent_dielectric,
            salt_concentration=salt_concentration,
            gb_offset=gb_offset,
        )

        gbsa_in_path = mol_out / "dock6_gbsa.in"
        gbsa_in_path.write_text(gbsa_in)

        # Run dock6
        t0 = time.time()
        try:
            proc = subprocess.run(
                ["dock6", "-i", "dock6_gbsa.in", "-o", "dock6_gbsa.out"],
                cwd=str(mol_out),
                capture_output=True, text=True,
                timeout=timeout_per_molecule,
            )
            runtime = time.time() - t0
            gbsa_mol2 = mol_out / f"{name}_gbsa_scored.mol2"
            success = proc.returncode == 0 and gbsa_mol2.exists() and gbsa_mol2.stat().st_size > 0

            results.append({
                "Name": name,
                "status": "OK" if success else "FAILED",
                "gbsa_mol2": str(gbsa_mol2) if success else None,
                "runtime_sec": round(runtime, 1),
                "error": proc.stderr[:200] if not success else None,
            })

            total_time += runtime
            if success:
                logger.info(f"    -> OK ({runtime:.1f}s, {n_poses} poses)")
            else:
                logger.warning(f"    -> FAILED ({runtime:.1f}s)")
                if proc.stderr:
                    logger.debug(f"    stderr: {proc.stderr[:300]}")

        except subprocess.TimeoutExpired:
            runtime = time.time() - t0
            total_time += runtime
            results.append({
                "Name": name, "status": "TIMEOUT", "gbsa_mol2": None,
                "runtime_sec": round(runtime, 1),
                "error": f"Timeout after {timeout_per_molecule}s",
            })
            logger.warning(f"    -> TIMEOUT")

        except Exception as e:
            runtime = time.time() - t0
            total_time += runtime
            results.append({
                "Name": name, "status": "FAILED", "gbsa_mol2": None,
                "runtime_sec": round(runtime, 1), "error": str(e),
            })
            logger.warning(f"    -> ERROR: {e}")

    # Save status
    df = pd.DataFrame(results)
    df.to_csv(output_dir / "gbsa_status.csv", index=False)

    n_ok = sum(1 for r in results if r["status"] == "OK")
    logger.info("")
    logger.info(f"  GB/SA RESCORE: {n_ok}/{len(results)} completed ({total_time:.0f}s)")

    return {
        "success": True,
        "n_total": len(results),
        "n_ok": n_ok,
        "n_failed": len(results) - n_ok,
        "total_runtime_sec": round(total_time, 1),
        "results": results,
    }
