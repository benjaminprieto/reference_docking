"""
DOCK6 Runner - Core Module (01c)
==================================
Genera input files de DOCK6 y ejecuta rigid re-docking para ligandos
cristalográficos. Soporta también flexible docking (anchor-and-grow).

Soporta:
  - Rigid ligand docking (default for reference)
  - Flexible ligand docking (anchor-and-grow)
  - Grid-based scoring (vdW + electrostatic)
  - Footprint scoring (per-residue energy decomposition)
  - GB/SA Hawkins implicit solvation (secondary score)
  - Simplex minimization

DOCK6 path limit:
  DOCK6 programs (Fortran legacy) truncate file paths at ~80 characters.
  This module creates symlinks in each molecule's output directory and
  uses short filenames in dock6.in. dock6 is executed with cwd=mol_out
  so all paths are relative and short.

Pipeline por molecula:
    1. Crear symlinks en mol_out/ para spheres, grids, param files, ligand
    2. Generar dock6.in con filenames cortos (no paths absolutos)
    3. Ejecutar dock6 -i dock6.in -o dock6.out con cwd=mol_out
    4. Verificar output ({name}_scored.mol2)

Input:
  - mol2 individuales (ligandos cristalográficos con coordenadas originales)
  - Grids DOCK6 (de 01b grid generation)
  - spheres_ligand.sph

Output (por molecula):
  - {name}/dock6.in
  - {name}/dock6.out
  - {name}/{name}_scored.mol2

Location: 01_src/reference_docking/m01_docking/dock6_runner.py
Project: reference_docking
Module: 01c (core)
Version: 3.0 — adapted for reference docking (2026-03-25)
"""

import logging
import os
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# DOCK6 INPUT TEMPLATES
# =============================================================================

DOCK6_FLEX_TEMPLATE = """\
conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              {min_anchor_size}
pruning_use_clustering                                       yes
pruning_max_orients                                          {pruning_max_orients}
pruning_clustering_cutoff                                    {pruning_clustering_cutoff}
pruning_conformer_score_cutoff                               {pruning_conformer_score_cutoff}
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            no
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             {ligand_mol2}
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           {spheres_file}
max_orientations                                             {max_orientations}
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       {grid_prefix}
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
{footprint_block}
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
{gbsa_hawkins_block}
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              {minimize}
simplex_max_iterations                                       {simplex_max_iterations}
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           {simplex_max_cycles}
simplex_score_converge                                       {simplex_score_converge}
simplex_cycle_converge                                       {simplex_cycle_converge}
simplex_trans_step                                           {simplex_trans_step}
simplex_rot_step                                             {simplex_rot_step}
simplex_tors_step                                            {simplex_tors_step}
simplex_random_seed                                          {simplex_random_seed}
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                {vdw_defn_file}
flex_defn_file                                               {flex_defn_file}
flex_drive_file                                              {flex_drive_file}
ligand_outfile_prefix                                        {output_prefix}
write_orientations                                           {write_orientations}
num_scored_conformers                                        {num_scored_conformers}
num_final_scored_poses                                       {num_final_scored_poses}
num_preclustered_conformers                                  {num_preclustered_conformers}
rank_ligands                                                 no
"""

DOCK6_RIGID_TEMPLATE = """\
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             {ligand_mol2}
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           {spheres_file}
max_orientations                                             {max_orientations}
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       {grid_prefix}
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
{footprint_block}
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
{gbsa_hawkins_block}
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              {minimize}
simplex_max_iterations                                       {simplex_max_iterations}
simplex_max_cycles                                           {simplex_max_cycles}
simplex_score_converge                                       {simplex_score_converge}
simplex_cycle_converge                                       {simplex_cycle_converge}
simplex_trans_step                                           {simplex_trans_step}
simplex_rot_step                                             {simplex_rot_step}
simplex_tors_step                                            {simplex_tors_step}
simplex_random_seed                                          {simplex_random_seed}
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                {vdw_defn_file}
flex_defn_file                                               {flex_defn_file}
flex_drive_file                                              {flex_drive_file}
ligand_outfile_prefix                                        {output_prefix}
write_orientations                                           {write_orientations}
num_scored_conformers                                        {num_scored_conformers}
num_final_scored_poses                                       {num_final_scored_poses}
num_preclustered_conformers                                  {num_preclustered_conformers}
rank_ligands                                                 no
"""


# =============================================================================
# DOCK6 PARAMETER FILE DISCOVERY
# =============================================================================

def find_dock6_params() -> Dict[str, str]:
    """
    Find DOCK6 parameter files (vdw_AMBER_parm99.defn, flex.defn, flex_drive.tbl).

    Search order:
      1. $DOCK_HOME/parameters/
      2. $DOCK6_HOME/parameters/
      3. $DOCK_BASE/parameters/
      4. Relative to dock6 binary: dirname(which dock6)/../parameters/
      5. Common install paths

    Returns:
        Dict with keys: vdw_defn_file, flex_defn_file, flex_drive_file
        Values are absolute paths if found, bare filenames otherwise.
    """
    param_files = {
        "vdw_defn_file": "vdw_AMBER_parm99.defn",
        "flex_defn_file": "flex.defn",
        "flex_drive_file": "flex_drive.tbl",
    }

    search_paths = []

    # Environment variables
    for var in ["DOCK_HOME", "DOCK6_HOME", "DOCK_BASE"]:
        val = os.environ.get(var)
        if val:
            search_paths.append(Path(val) / "parameters")

    # Relative to dock6 binary
    try:
        dock6_which = subprocess.run(
            ["which", "dock6"], capture_output=True, text=True, timeout=5,
        ).stdout.strip()
        if dock6_which:
            # dock6 at /path/to/dock6/bin/dock6 -> parameters at /path/to/dock6/parameters/
            search_paths.append(Path(dock6_which).resolve().parent.parent / "parameters")
    except Exception:
        pass

    # Common install paths
    search_paths.extend([
        Path("/opt/dock6/parameters"),
        Path("/usr/local/dock6/parameters"),
        Path.home() / "dock6" / "parameters",
        Path.home() / "software" / "dock6" / "parameters",
        Path.home() / "programs" / "dock6" / "parameters",
    ])

    result = {}
    for key, filename in param_files.items():
        result[key] = filename  # Default: bare name (hope it's findable)
        for search_dir in search_paths:
            candidate = search_dir / filename
            if candidate.exists():
                result[key] = str(candidate.resolve())
                break

    # Log what we found
    found = sum(1 for v in result.values() if "/" in v or "\\" in v)
    logger.info(f"  DOCK6 parameters: {found}/{len(param_files)} found with full paths")
    for key, path in result.items():
        logger.debug(f"    {key}: {path}")

    return result


# =============================================================================
# GRID PREFIX RESOLUTION
# =============================================================================

def resolve_grid_prefix(grid_dir: str, energy_grid: str) -> str:
    """
    Resolve the grid prefix for DOCK6.

    DOCK6 expects a prefix like "/path/to/grids/grid" and appends
    .nrg and .bmp itself.

    Args:
        grid_dir: Directory containing grid files
        energy_grid: Energy grid filename (e.g., "ligand.nrg")

    Returns:
        Full path prefix (without extension)
    """
    prefix = energy_grid.replace(".nrg", "").replace(".NRG", "")
    return str(Path(grid_dir) / prefix)


# =============================================================================
# INPUT FILE GENERATION
# =============================================================================

def generate_dock6_input(
        ligand_mol2: str,
        spheres_file: str,
        grid_prefix: str,
        output_prefix: str,
        output_path: str,
        search_method: str = "flex",
        dock6_params: Optional[Dict[str, str]] = None,
        **kwargs,
) -> str:
    """
    Generate a DOCK6 input file.

    Args:
        ligand_mol2: Path to ligand mol2 file
        spheres_file: Path to spheres_ligand.sph
        grid_prefix: Grid prefix (without .nrg/.bmp)
        output_prefix: Prefix for output scored mol2
        output_path: Path to write the dock6.in file
        search_method: "flex" or "rigid"
        dock6_params: DOCK6 parameter file paths (from find_dock6_params)
        **kwargs: Override any template variable

    Returns:
        Path to the generated input file
    """
    if dock6_params is None:
        dock6_params = find_dock6_params()

    template = DOCK6_FLEX_TEMPLATE if search_method == "flex" else DOCK6_RIGID_TEMPLATE

    # Build template variables with defaults
    template_vars = {
        # Paths
        "ligand_mol2": ligand_mol2,
        "spheres_file": spheres_file,
        "grid_prefix": grid_prefix,
        "output_prefix": output_prefix,
        # DOCK6 parameter files
        "vdw_defn_file": dock6_params.get("vdw_defn_file", "vdw_AMBER_parm99.defn"),
        "flex_defn_file": dock6_params.get("flex_defn_file", "flex.defn"),
        "flex_drive_file": dock6_params.get("flex_drive_file", "flex_drive.tbl"),
        # Flex-specific
        "min_anchor_size": 5,
        "pruning_max_orients": 1000,
        "pruning_clustering_cutoff": 100,
        "pruning_conformer_score_cutoff": 100.0,
        # Orientation
        "max_orientations": 1000,
        # Minimization
        "minimize": "yes",
        "simplex_max_iterations": 500,
        "simplex_max_cycles": 1,
        "simplex_score_converge": 0.1,
        "simplex_cycle_converge": 1.0,
        "simplex_trans_step": 1.0,
        "simplex_rot_step": 0.1,
        "simplex_tors_step": 10.0,
        "simplex_random_seed": 0,
        # Footprint (per-residue energy decomposition)
        "footprint_block": "footprint_similarity_score_secondary                         no",
        # GB/SA Hawkins implicit solvation (secondary score)
        "gbsa_hawkins_block": "gbsa_hawkins_score_secondary                                 no",
        # Output
        "num_scored_conformers": 20,
        "num_final_scored_poses": 100,
        "num_preclustered_conformers": 500,
        "write_orientations": "no",
    }

    # Apply kwargs overrides
    for key, val in kwargs.items():
        if key in template_vars:
            # Convert booleans to DOCK6 yes/no
            if isinstance(val, bool):
                val = "yes" if val else "no"
            template_vars[key] = val

    # Handle footprint scoring
    receptor_mol2 = kwargs.pop("receptor_mol2", None)
    reference_mol2 = kwargs.pop("reference_mol2", None)
    compute_footprint = kwargs.pop("compute_footprint_score", False)
    if compute_footprint and receptor_mol2 and reference_mol2:
        template_vars["footprint_block"] = (
            "footprint_similarity_score_secondary                         yes\n"
            f"fps_score_use_footprint_reference_mol2                       yes\n"
            f"fps_score_footprint_reference_mol2_filename                  {reference_mol2}\n"
            "fps_score_foot_compare_type                                  Euclidean\n"
            "fps_score_normalize_foot                                     no\n"
            "fps_score_foot_comp_all_residue                              yes\n"
            f"fps_score_receptor_filename                                  {receptor_mol2}\n"
            "fps_score_vdw_att_exp                                        6\n"
            "fps_score_vdw_rep_exp                                        9\n"
            "fps_score_vdw_rep_rad_scale                                  1\n"
            "fps_score_use_distance_dependent_dielectric                  yes\n"
            "fps_score_dielectric                                         4.0\n"
            "fps_score_vdw_fp_scale                                       1\n"
            "fps_score_es_fp_scale                                        1\n"
            "fps_score_hb_fp_scale                                        0"
        )

    # Handle GB/SA Hawkins implicit solvation (secondary score)
    gbsa_hawkins = kwargs.pop("gbsa_hawkins", False)
    gbsa_solvent_dielectric = kwargs.pop("solvent_dielectric", 78.5)
    gbsa_salt_concentration = kwargs.pop("salt_concentration", 0.15)
    gbsa_gb_offset = kwargs.pop("gb_offset", 0.09)
    if gbsa_hawkins and receptor_mol2:
        template_vars["gbsa_hawkins_block"] = (
            "gbsa_hawkins_score_secondary                                 yes\n"
            f"gbsa_hawkins_score_rec_filename                              {receptor_mol2}\n"
            f"gbsa_hawkins_score_solvent_dielectric                        {gbsa_solvent_dielectric}\n"
            f"gbsa_hawkins_score_salt_conc                                 {gbsa_salt_concentration}\n"
            f"gbsa_hawkins_score_gb_offset                                 {gbsa_gb_offset}\n"
            "gbsa_hawkins_score_cont_vdw_and_es                           yes\n"
            "gbsa_hawkins_score_vdw_att_exp                               6\n"
            "gbsa_hawkins_score_vdw_rep_exp                               12"
        )

    # Handle 'minimize' bool -> string
    if isinstance(template_vars["minimize"], bool):
        template_vars["minimize"] = "yes" if template_vars["minimize"] else "no"
    elif template_vars["minimize"] is True:
        template_vars["minimize"] = "yes"
    elif template_vars["minimize"] is False:
        template_vars["minimize"] = "no"

    if isinstance(template_vars["write_orientations"], bool):
        template_vars["write_orientations"] = "yes" if template_vars["write_orientations"] else "no"

    content = template.format(**template_vars)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(content)

    return output_path


# =============================================================================
# SYMLINK HELPER (80-char path workaround)
# =============================================================================

def _create_symlink(target: str, link_path: Path) -> bool:
    """
    Create a symlink, removing any existing link/file.

    Args:
        target: Absolute path to the real file
        link_path: Path where the symlink should be created

    Returns:
        True if symlink was created successfully
    """
    target_path = Path(target).resolve()
    if not target_path.exists():
        logger.debug(f"    Symlink target not found: {target_path}")
        return False

    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()

    link_path.symlink_to(target_path)
    return True


def _setup_mol_symlinks(
        mol_out: Path,
        mol2_file: Path,
        spheres_file: str,
        grid_prefix: str,
        dock6_params: Dict[str, str],
        receptor_mol2: Optional[str] = None,
        reference_mol2: Optional[str] = None,
) -> Dict[str, str]:
    """
    Create symlinks in the molecule output directory for all files
    referenced in dock6.in. Returns a dict of short filenames to use
    in the dock6.in template.

    This is the fix for DOCK6's 80-character path limit.
    All symlinks point to resolved absolute paths, and dock6.in
    uses only the short filenames with cwd=mol_out.
    """
    short_names = {}

    # Ligand mol2
    lig_link = mol_out / mol2_file.name
    _create_symlink(str(mol2_file), lig_link)
    short_names["ligand_mol2"] = mol2_file.name

    # Spheres
    _create_symlink(spheres_file, mol_out / "spheres_ligand.sph")
    short_names["spheres_file"] = "spheres_ligand.sph"

    # Grid files (prefix.nrg, prefix.bmp)
    grid_nrg = grid_prefix + ".nrg"
    grid_bmp = grid_prefix + ".bmp"
    _create_symlink(grid_nrg, mol_out / "ligand.nrg")
    _create_symlink(grid_bmp, mol_out / "ligand.bmp")
    short_names["grid_prefix"] = "ligand"

    # DOCK6 parameter files (vdw_defn, flex_defn, flex_drive)
    for key in ["vdw_defn_file", "flex_defn_file", "flex_drive_file"]:
        src = dock6_params.get(key, "")
        if "/" in src or "\\" in src:
            link_name = Path(src).name
            _create_symlink(src, mol_out / link_name)
            short_names[key] = link_name
        else:
            short_names[key] = src

    # Receptor mol2 (for footprint scoring)
    if receptor_mol2 and Path(receptor_mol2).exists():
        rec_name = Path(receptor_mol2).name
        _create_symlink(receptor_mol2, mol_out / rec_name)
        short_names["receptor_mol2"] = rec_name

    # Reference mol2 (for footprint comparison)
    if reference_mol2 and Path(reference_mol2).exists():
        ref_name = "fps_reference.mol2"
        _create_symlink(reference_mol2, mol_out / ref_name)
        short_names["reference_mol2"] = ref_name

    return short_names

def run_dock6_single(
        dock_input: str,
        dock_output: str,
        timeout: int = 600,
        cwd: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Execute dock6 for a single input file.

    Args:
        dock_input: Path to dock6.in (filename if cwd is set)
        dock_output: Path for dock6.out (filename if cwd is set)
        timeout: Maximum seconds to wait
        cwd: Working directory for dock6 execution (for short paths)

    Returns:
        Dict with: success, returncode, runtime_sec, error
    """
    start = time.time()

    try:
        result = subprocess.run(
            ["dock6", "-i", dock_input, "-o", dock_output],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=cwd,
        )
        elapsed = time.time() - start

        return {
            "success": result.returncode == 0,
            "returncode": result.returncode,
            "runtime_sec": round(elapsed, 1),
            "error": result.stderr[:500] if result.returncode != 0 and result.stderr else None,
        }

    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "returncode": -1,
            "runtime_sec": round(time.time() - start, 1),
            "error": f"Timeout after {timeout}s",
        }
    except FileNotFoundError:
        return {
            "success": False,
            "returncode": -2,
            "runtime_sec": 0,
            "error": "dock6 not found in PATH",
        }


# =============================================================================
# VALIDATION
# =============================================================================

def validate_grids(spheres_file: str, grid_prefix: str) -> List[str]:
    """
    Validate that all required grid files exist.

    Returns:
        List of error messages (empty = all OK)
    """
    errors = []

    if not Path(spheres_file).exists():
        errors.append(f"Spheres not found: {spheres_file}")

    nrg = grid_prefix + ".nrg"
    if not Path(nrg).exists():
        errors.append(f"Energy grid not found: {nrg}")

    bmp = grid_prefix + ".bmp"
    if not Path(bmp).exists():
        errors.append(f"Bump grid not found: {bmp}")

    return errors


def validate_dock6_available() -> bool:
    """Check if dock6 binary is available."""
    try:
        result = subprocess.run(
            ["which", "dock6"], capture_output=True, text=True, timeout=5,
        )
        return result.returncode == 0
    except Exception:
        return False


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

def run_dock6_batch(
        ligand_mol2_dir: Union[str, Path],
        spheres_file: Union[str, Path],
        grid_prefix: str,
        output_dir: Union[str, Path],
        search_method: str = "flex",
        max_orientations: int = 1000,
        num_scored_conformers: int = 20,
        minimize: bool = True,
        simplex_max_iterations: int = 500,
        timeout_per_molecule: int = 600,
        molecule_filter: Optional[List[str]] = None,
        receptor_mol2: Optional[str] = None,
        reference_mol2: Optional[str] = None,
        dry_run: bool = False,
        **kwargs,
) -> Dict[str, Any]:
    """
    Run DOCK6 for all prepared ligands in a directory.

    Scans ligand_mol2_dir for *.mol2 files. For each, generates a
    dock6.in and executes dock6.

    Args:
        ligand_mol2_dir: Directory with prepared mol2 files (from 00c)
        spheres_file: Path to spheres_ligand.sph
        grid_prefix: Grid prefix (without .nrg/.bmp)
        output_dir: Directory for docking output
        search_method: "flex" | "rigid"
        max_orientations: Maximum orientations
        num_scored_conformers: Legacy param (DOCK6 <6.13). Also set via kwargs:
            num_final_scored_poses (default 100) and
            num_preclustered_conformers (default 500) for DOCK6 6.13+
        minimize: Run simplex minimization
        simplex_max_iterations: Max minimization iterations
        timeout_per_molecule: Timeout in seconds per molecule
        molecule_filter: Only dock these molecules (None = all)
        dry_run: Generate input files only, don't execute dock6
        **kwargs: Additional DOCK6 parameters (min_anchor_size, etc.)

    Returns:
        Dict with per-molecule status and summary
    """
    ligand_dir = Path(ligand_mol2_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    spheres_file = str(spheres_file)

    # --- Validate prerequisites ---
    if not ligand_dir.exists():
        logger.error(f"Ligand directory not found: {ligand_dir}")        logger.error("Check ligand_mol2 path in campaign_config.yaml.")
        return {"n_total": 0, "n_ok": 0, "n_failed": 0, "error": "Ligand dir not found"}

    grid_errors = validate_grids(spheres_file, grid_prefix)
    if grid_errors:
        for err in grid_errors:
            logger.error(err)
        logger.error("Check grids.grid_dir in campaign_config.yaml")
        return {"n_total": 0, "n_ok": 0, "n_failed": 0, "error": "; ".join(grid_errors)}

    if not dry_run and not validate_dock6_available():
        logger.error("dock6 binary not found in PATH")
        return {"n_total": 0, "n_ok": 0, "n_failed": 0, "error": "dock6 not in PATH"}

    # --- Find ligand mol2 files ---
    mol2_files = sorted(ligand_dir.glob("*.mol2"))

    # Skip empty or broken mol2 files
    valid_mol2 = []
    for f in mol2_files:
        if f.stat().st_size == 0:
            logger.warning(f"  Skipping empty mol2: {f.name}")
        else:
            valid_mol2.append(f)
    mol2_files = valid_mol2

    if molecule_filter:
        filter_set = set(molecule_filter)
        mol2_files = [f for f in mol2_files if f.stem in filter_set]

    if not mol2_files:
        logger.error(f"No mol2 files found in {ligand_dir}")
        if molecule_filter:
            logger.error(f"  Filter: {molecule_filter}")
        return {"n_total": 0, "n_ok": 0, "n_failed": 0, "error": "No mol2 files"}

    logger.info(f"Found {len(mol2_files)} ligands to dock")
    logger.info(f"Grid prefix: {grid_prefix}")
    logger.info(f"Spheres:     {spheres_file}")
    logger.info(f"Method:      {search_method}")
    logger.info(f"Orientations:{max_orientations}")
    if dry_run:
        logger.info("*** DRY RUN — generating input files only ***")

    # --- Find DOCK6 parameter files ---
    dock6_params = find_dock6_params()

    # --- Process each molecule ---
    results = []
    total_time = 0

    for idx, mol2_file in enumerate(mol2_files):
        name = mol2_file.stem
        mol_out = output_dir / name
        mol_out.mkdir(parents=True, exist_ok=True)

        logger.info(f"  [{idx + 1}/{len(mol2_files)}] {name}")

        # --- Create symlinks for short paths (80-char fix) ---
        short = _setup_mol_symlinks(
            mol_out=mol_out,
            mol2_file=mol2_file,
            spheres_file=spheres_file,
            grid_prefix=grid_prefix,
            dock6_params=dock6_params,
            receptor_mol2=receptor_mol2,
            reference_mol2=reference_mol2,
        )

        # --- Generate input file with short filenames ---
        dock_input = "dock6.in"
        dock_output = "dock6.out"
        output_prefix = name  # Just the name, no path (cwd=mol_out)
        scored_mol2 = str(mol_out / f"{name}_scored.mol2")

        generate_dock6_input(
            ligand_mol2=short["ligand_mol2"],
            spheres_file=short["spheres_file"],
            grid_prefix=short["grid_prefix"],
            output_prefix=output_prefix,
            output_path=str(mol_out / dock_input),
            search_method=search_method,
            dock6_params=short,  # short filenames for param files
            max_orientations=max_orientations,
            num_scored_conformers=num_scored_conformers,
            minimize=minimize,
            simplex_max_iterations=simplex_max_iterations,
            receptor_mol2=short.get("receptor_mol2"),
            reference_mol2=short.get("reference_mol2"),
            **kwargs,
        )

        if dry_run:
            results.append({
                "Name": name,
                "status": "DRY_RUN",
                "dock_input": str(mol_out / dock_input),
                "scored_mol2": None,
                "runtime_sec": 0,
                "error": None,
            })
            logger.info(f"    -> DRY_RUN (input generated)")
            continue

        # --- Execute DOCK6 with cwd=mol_out (short paths) ---
        run_result = run_dock6_single(
            dock_input, dock_output,
            timeout=timeout_per_molecule,
            cwd=str(mol_out),
        )
        has_output = Path(scored_mol2).exists() and Path(scored_mol2).stat().st_size > 0

        status = "OK" if run_result["success"] and has_output else "FAILED"
        error = None
        if not run_result["success"]:
            error = run_result.get("error", "Unknown error")
        elif not has_output:
            error = "dock6 succeeded but no scored mol2 produced"

        results.append({
            "Name": name,
            "status": status,
            "dock_input": str(mol_out / dock_input),
            "dock_output": str(mol_out / dock_output),
            "scored_mol2": scored_mol2 if has_output else None,
            "runtime_sec": run_result["runtime_sec"],
            "error": error,
        })

        total_time += run_result["runtime_sec"]

        if status == "OK":
            logger.info(f"    -> OK ({run_result['runtime_sec']}s)")
        else:
            logger.warning(f"    -> FAILED: {error}")

    # --- Save status CSV ---
    df_status = pd.DataFrame(results)
    status_csv = output_dir / "docking_status.csv"
    df_status.to_csv(status_csv, index=False, encoding="utf-8")
    logger.info(f"  Saved: {status_csv}")

    # --- Save summary TXT ---
    n_ok = sum(1 for r in results if r["status"] == "OK")
    n_dry = sum(1 for r in results if r["status"] == "DRY_RUN")
    n_fail = sum(1 for r in results if r["status"] == "FAILED")

    summary_path = output_dir / "docking_summary.txt"
    w = 70
    lines = [
        "=" * w,
        "01b DOCK6 RUN - SUMMARY",
        "=" * w,
        "",
        f"Date:              {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Molecules:         {len(results)}",
        f"Successful:        {n_ok}",
        f"Failed:            {n_fail}",
        f"Dry run:           {n_dry}",
        f"Total time:        {total_time:.0f}s ({total_time / 60:.1f} min)",
        f"Avg per molecule:  {total_time / max(n_ok, 1):.1f}s",
        "",
        f"Method:            {search_method}",
        f"Orientations:      {max_orientations}",
        f"Minimize:          {minimize}",
        f"Max iterations:    {simplex_max_iterations}",
        "",
        "-" * w,
        f"{'Name':<30} {'Status':>8} {'Time(s)':>8} {'Scored mol2':<20}",
        "-" * w,
    ]

    for r in results:
        scored = "yes" if r.get("scored_mol2") else "—"
        lines.append(
            f"{r['Name']:<30} {r['status']:>8} {r['runtime_sec']:>8.1f} {scored:<20}"
        )

    if n_fail > 0:
        lines.extend(["", "FAILURES:"])
        for r in results:
            if r["status"] == "FAILED":
                lines.append(f"  {r['Name']}: {r.get('error', 'unknown')}")

    lines.extend(["", "=" * w])

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    logger.info("")
    logger.info(f"{'=' * 60}")
    if dry_run:
        logger.info(f"  DRY RUN: {n_dry} input files generated")
    else:
        logger.info(f"  RESULTS: {n_ok}/{len(results)} dockings completed "
                     f"({total_time:.0f}s total)")
        if n_fail > 0:
            logger.info(f"  FAILED:  {n_fail}")
    logger.info(f"{'=' * 60}")

    return {
        "n_total": len(results),
        "n_ok": n_ok,
        "n_failed": n_fail,
        "n_dry_run": n_dry,
        "total_runtime_sec": round(total_time, 1),
        "results": results,
        "status_csv": str(status_csv),
        "summary_txt": str(summary_path),
        "output_dir": str(output_dir),
    }