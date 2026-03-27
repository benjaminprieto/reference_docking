"""
Grid Generation - Core Module (01a)
=====================================
Genera grids DOCK6 a partir del receptor preparado.

OPCIONAL: Si grids pre-existentes son validos en campaign_config.grids.grid_dir,
este modulo se SALTA automaticamente.

DOCK6 grid generation pipeline:
    1. DMS surface   (dms)             -> receptor.ms
    2. Spheres       (sphgen)          -> all_spheres.sph
    3. Selection     (sphere_selector) -> spheres_ligand.sph
    4. Box           (showbox)         -> spheres_ligand_box.pdb
    5. Grid          (grid)            -> ligand.nrg, ligand.bmp

Binding site definition (step 3) supports 4 methods:
    A. reference_ligand  -> sphere_selector uses ligand mol2/pdb
    B. residues          -> centroid of specified residues
    C. coordinates       -> explicit (x, y, z) center
    D. (skip)            -> grids already exist

Requires DOCK6 accessory programs in PATH:
    dms, sphgen, sphere_selector, showbox, grid

Input:
  - rec_noH.pdb (receptor without H, for DMS)
  - rec_charged.mol2 (receptor with charges, for grid)
  - Binding site reference (ligand, residues, or coordinates)

Output (all visualizable in ChimeraX):
  - receptor.ms               DMS surface
  - all_spheres.sph           All spheres (sphgen)
  - INSPH, OUTSPH             sphgen input/output logs
  - spheres_ligand.sph        Selected spheres near binding site
  - box.in                    showbox input (reproducibility)
  - spheres_ligand_box.pdb    Grid box as PDB (open in ChimeraX)
  - grid.in                   grid input
  - grid.out                  grid output (diagnostics)
  - ligand.nrg, ligand.bmp    Energy and bump grids
  - grid_generation_report.json

Location: 01_src/reference_docking/m01_docking/grid_generation.py
Project: reference_docking
Module: 01a (core)
Version: 2.2 — tutorial-compatible output names (2026-03-13)
"""

import json
import logging
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple

logger = logging.getLogger(__name__)


# =============================================================================
# VALIDATION
# =============================================================================

def validate_existing_grids(
        grid_dir: Union[str, Path],
        spheres_file: str = "spheres_ligand.sph",
        energy_grid: str = "ligand.nrg",
        bump_grid: str = "ligand.bmp",
) -> bool:
    """
    Check if pre-existing grids are valid and complete.

    Returns True if all required files exist and are non-empty.
    """
    grid_dir = Path(grid_dir)
    required = [
        grid_dir / spheres_file,
        grid_dir / energy_grid,
        grid_dir / bump_grid,
    ]
    return all(f.exists() and f.stat().st_size > 0 for f in required)


def check_dock6_tools() -> Dict[str, bool]:
    """Check which DOCK6 accessory programs are available."""
    tools = {}
    for tool in ["dms", "sphgen", "sphere_selector", "showbox", "grid"]:
        try:
            result = subprocess.run(
                ["which", tool], capture_output=True, text=True, timeout=5,
            )
            tools[tool] = result.returncode == 0
        except Exception:
            tools[tool] = False
    return tools


# =============================================================================
# BINDING SITE: CENTROID FROM RESIDUES
# =============================================================================

def compute_residue_centroid(
        pdb_path: str,
        residues: List[str],
        chain: Optional[str] = None,
) -> Tuple[float, float, float]:
    """
    Compute the geometric centroid of specified residues from a PDB file.

    Reads CA (alpha carbon) atoms for each residue. Falls back to any
    atom if CA not found.

    Args:
        pdb_path: Path to PDB file
        residues: Residue identifiers, e.g. ["GLU529", "TRP555", "CYS574"]
                  Format: 3-letter-code + residue number
        chain: Chain ID filter (None = any chain)

    Returns:
        (x, y, z) centroid coordinates

    Raises:
        ValueError: If no matching atoms found
    """
    import re

    # Parse residue specs: "GLU529" -> ("GLU", 529)
    residue_specs = []
    for r in residues:
        match = re.match(r"([A-Z]{3})(\d+)", r)
        if match:
            residue_specs.append((match.group(1), int(match.group(2))))
        else:
            logger.warning(f"  Cannot parse residue: {r} (expected format: GLU529)")

    if not residue_specs:
        raise ValueError(f"No valid residues parsed from: {residues}")

    # Read PDB and collect coordinates (CA preferred, any atom as fallback)
    ca_coords: Dict[str, Tuple[float, float, float]] = {}
    fallback_coords: Dict[str, Tuple[float, float, float]] = {}

    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            res_num = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if chain and chain_id != chain:
                continue

            for spec_name, spec_num in residue_specs:
                if res_name == spec_name and res_num == spec_num:
                    res_key = f"{res_name}{res_num}"
                    if atom_name == "CA":
                        ca_coords[res_key] = (x, y, z)
                    elif res_key not in fallback_coords:
                        fallback_coords[res_key] = (x, y, z)

    # Merge: prefer CA, fallback to first atom
    coords = []
    found_residues = set()
    for spec_name, spec_num in residue_specs:
        res_key = f"{spec_name}{spec_num}"
        if res_key in ca_coords:
            coords.append(ca_coords[res_key])
            found_residues.add(res_key)
        elif res_key in fallback_coords:
            coords.append(fallback_coords[res_key])
            found_residues.add(res_key)

    if not coords:
        raise ValueError(
            f"No atoms found for residues {residues} in {pdb_path}. "
            f"Check residue names/numbers and chain ID."
        )

    # Compute centroid
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)

    logger.info(f"  Binding site centroid from {len(found_residues)} residues: "
                f"({cx:.3f}, {cy:.3f}, {cz:.3f})")

    return (cx, cy, cz)


# =============================================================================
# STEP 1: DMS SURFACE
# =============================================================================

def generate_dms_surface(
        receptor_noH_pdb: str,
        output_dir: str,
        probe_radius: float = 1.4,
) -> Optional[str]:
    """
    Generate molecular surface using dms.

    Args:
        receptor_noH_pdb: Receptor PDB without hydrogens
        output_dir: Output directory
        probe_radius: Probe radius in Angstroms

    Returns:
        Path to receptor.ms, or None if failed
    """
    output_path = Path(output_dir) / "receptor.ms"

    cmd = [
        "dms", receptor_noH_pdb,
        "-n",  # use normals
        "-w", str(probe_radius),
        "-o", str(output_path),
    ]

    logger.info(f"  Step 1: DMS surface (probe={probe_radius} A)")
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300,
        )
        if result.returncode == 0 and output_path.exists():
            logger.info(f"    -> {output_path.name} ({output_path.stat().st_size} bytes)")
            return str(output_path)
        else:
            logger.error(f"    dms failed (rc={result.returncode})")
            if result.stderr:
                logger.error(f"    {result.stderr[:300]}")
            return None
    except FileNotFoundError:
        logger.error("    dms not found in PATH")
        return None
    except subprocess.TimeoutExpired:
        logger.error("    dms timed out")
        return None


# =============================================================================
# STEP 2: SPHERE GENERATION
# =============================================================================

def generate_spheres(
        dms_path: str,
        output_dir: str,
        steric_close_distance: float = 0.0,
        min_radius: float = 1.4,
        max_radius: float = 4.0,
) -> Optional[str]:
    """
    Generate spheres using sphgen.

    sphgen reads from a file called INSPH in the current working directory.
    It creates temporary files (temp1.ms, temp2.sph, temp3.atc) and writes
    the sphere output to the filename specified in INSPH (last line).
    OUTSPH is the log/info file.

    CRITICAL: sphgen (Fortran) will fail if OUTSPH or temp files already
    exist. They must be deleted before each run.

    Returns:
        Path to all_spheres.sph, or None if failed
    """
    output_path = Path(output_dir) / "all_spheres.sph"
    sph_filename = "all_spheres.sph"
    dms_filename = Path(dms_path).name

    # --- Clean up files that sphgen refuses to overwrite ---
    cleanup_files = [
        "OUTSPH", "temp1.ms", "temp2.sph", "temp3.atc",
        sph_filename,
    ]
    for fname in cleanup_files:
        fpath = Path(output_dir) / fname
        if fpath.exists():
            fpath.unlink()
            logger.debug(f"    Cleaned up: {fname}")

    # --- Write INSPH (format from DOCK6 documentation) ---
    # Line 1: DMS surface file
    # Line 2: R (receptor) or L (ligand)
    # Line 3: X (all surface points) or a specific subset
    # Line 4: steric close distance (0.0 = no filter)
    # Line 5: maximum sphere radius
    # Line 6: minimum sphere radius
    # Line 7: output sphere file
    insph_content = f"""{dms_filename}
R
X
{steric_close_distance}
{max_radius}
{min_radius}
{sph_filename}
"""
    insph_path = Path(output_dir) / "INSPH"
    insph_path.write_text(insph_content)

    logger.info(f"  Step 2: Sphere generation (sphgen)")
    logger.debug(f"    INSPH: dms={dms_filename}, R={max_radius}, r={min_radius}, out={sph_filename}")

    try:
        result = subprocess.run(
            ["sphgen", "-i", "INSPH", "-o", "OUTSPH"],
            capture_output=True, text=True, timeout=600,
            cwd=output_dir,
        )

        # Log any output for debugging
        if result.stdout:
            logger.debug(f"    sphgen stdout: {result.stdout[:500]}")
        if result.stderr:
            logger.warning(f"    sphgen stderr: {result.stderr[:500]}")
        if result.returncode != 0:
            logger.warning(f"    sphgen returncode: {result.returncode}")

        # Check if sphere file was created
        if output_path.exists() and output_path.stat().st_size > 0:
            # Count clusters from OUTSPH
            outsph = Path(output_dir) / "OUTSPH"
            if outsph.exists():
                outsph_text = outsph.read_text()
                logger.debug(f"    OUTSPH: {outsph_text[:300]}")
            logger.info(f"    -> {output_path.name} ({output_path.stat().st_size} bytes)")
            return str(output_path)
        else:
            # Log OUTSPH for diagnostics
            outsph = Path(output_dir) / "OUTSPH"
            if outsph.exists() and outsph.stat().st_size > 0:
                logger.error(f"    OUTSPH content: {outsph.read_text()[:500]}")
            logger.error("    sphgen produced no sphere output")
            return None

    except FileNotFoundError:
        logger.error("    sphgen not found in PATH")
        return None
    except subprocess.TimeoutExpired:
        logger.error("    sphgen timed out (>600s)")
        return None


# =============================================================================
# STEP 3: SPHERE SELECTION
# =============================================================================

def select_spheres_by_ligand(
        all_spheres: str,
        ligand_file: str,
        output_path: str,
        radius: float = 10.0,
) -> Optional[str]:
    """
    Method A: Select spheres near a reference ligand.

    Uses sphere_selector: sphere_selector all_spheres.sph ligand.mol2 radius
    """
    logger.info(f"  Step 3a: Sphere selection by reference ligand (radius={radius} A)")

    # Use absolute paths for sphere_selector arguments
    all_spheres_abs = str(Path(all_spheres).resolve())
    ligand_abs = str(Path(ligand_file).resolve())
    target = Path(output_path)
    target_dir = str(target.parent)

    try:
        result = subprocess.run(
            ["sphere_selector", all_spheres_abs, ligand_abs, str(radius)],
            capture_output=True, text=True, timeout=60,
            cwd=target_dir,
        )

        # sphere_selector always writes selected_spheres.sph to cwd
        # We rename it to spheres_ligand.sph (our target name)
        cwd_output = Path(target_dir) / "selected_spheres.sph"
        target = Path(output_path)

        if cwd_output.exists() and cwd_output != target:
            cwd_output.rename(target)

        if target.exists() and target.stat().st_size > 0:
            logger.info(f"    -> {target.name}")
            return str(target)
        else:
            logger.error("    sphere_selector produced no output")
            return None
    except FileNotFoundError:
        logger.error("    sphere_selector not found in PATH")
        return None


def select_spheres_by_center(
        all_spheres: str,
        center: Tuple[float, float, float],
        output_path: str,
        radius: float = 10.0,
) -> Optional[str]:
    """
    Method B/C: Select spheres within radius of a center point.

    sphere_selector needs a mol2 file as reference, so we create a
    dummy single-atom mol2 at the specified coordinates.
    """
    logger.info(f"  Step 3b: Sphere selection by center "
                f"({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}), radius={radius} A")

    # Create a dummy mol2 with a single atom at the center
    dummy_mol2 = f"""\
@<TRIPOS>MOLECULE
DUMMY
 1 0 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 DU        {center[0]:10.4f}{center[1]:10.4f}{center[2]:10.4f} Du      1 DUM         0.0000
@<TRIPOS>BOND
"""
    with tempfile.NamedTemporaryFile(
        suffix=".mol2", mode="w", delete=False, dir=Path(output_path).parent,
    ) as tmp:
        tmp.write(dummy_mol2)
        dummy_path = tmp.name

    try:
        result = select_spheres_by_ligand(all_spheres, dummy_path, output_path, radius)
        return result
    finally:
        Path(dummy_path).unlink(missing_ok=True)


# =============================================================================
# STEP 4: BOX GENERATION
# =============================================================================

def generate_box(
        selected_spheres: str,
        output_dir: str,
        margin: float = 10.0,
) -> Optional[str]:
    """
    Generate grid box using showbox.

    showbox reads from an input file interactively. We create the input
    and pipe it.
    """
    box_path = Path(output_dir) / "spheres_ligand_box.pdb"

    # showbox input: use filenames only (runs with cwd=output_dir)
    spheres_filename = Path(selected_spheres).name
    box_filename = box_path.name

    # showbox input format (5 lines):
    # 1. Y/N  - use spheres to define box
    # 2. margin - extra margin around spheres (Angstroms)
    # 3. sphere file
    # 4. cluster number (1 = largest)
    # 5. output box file
    showbox_input = f"Y\n{margin}\n{spheres_filename}\n1\n{box_filename}\n"

    # Save box.in for reproducibility
    box_in_path = Path(output_dir) / "box.in"
    box_in_path.write_text(showbox_input)

    logger.info(f"  Step 4: Box generation (margin={margin} A)")
    try:
        result = subprocess.run(
            ["showbox"],
            input=showbox_input,
            capture_output=True, text=True, timeout=60,
            cwd=output_dir,
        )

        if result.stderr:
            logger.warning(f"    showbox stderr: {result.stderr[:300]}")

        if box_path.exists() and box_path.stat().st_size > 0:
            logger.info(f"    -> {box_path.name}")
            return str(box_path)
        else:
            logger.error("    showbox produced no output")
            return None
    except FileNotFoundError:
        logger.error("    showbox not found in PATH")
        return None


# =============================================================================
# STEP 5: GRID CALCULATION
# =============================================================================

def generate_grid(
        receptor_charged_mol2: str,
        box_path: str,
        output_dir: str,
        grid_prefix: str = "grid",
        grid_spacing: float = 0.3,
        energy_cutoff_distance: float = 9999.0,
        attractive_exponent: int = 6,
        repulsive_exponent: int = 12,
        dielectric_factor: int = 4,
        bump_overlap: float = 0.75,
        vdw_defn_file: Optional[str] = None,
        dock6_home: Optional[str] = None,
) -> Optional[Dict[str, str]]:
    """
    Generate energy and bump grids using grid program.

    This is the most time-consuming step (~5-45 min depending on box size).

    Args:
        dock6_home: Path to DOCK6 installation (e.g. /opt/dock6).
                    Used to find parameter files. Takes priority over
                    environment variables.
    """
    import os

    # Find vdw definition file — search order:
    # 1. Explicit dock6_home from config YAML
    # 2. Environment variables (DOCK_HOME, DOCK6_HOME, DOCK_BASE)
    # 3. Relative to grid binary
    # 4. Common install paths
    if vdw_defn_file is None:
        search_paths = []

        # 1. From config
        if dock6_home:
            search_paths.append(Path(dock6_home) / "parameters")

        # 2. Environment variables
        for var in ["DOCK_HOME", "DOCK6_HOME", "DOCK_BASE"]:
            val = os.environ.get(var)
            if val:
                search_paths.append(Path(val) / "parameters")

        # 3. Relative to grid binary
        try:
            grid_which = subprocess.run(
                ["which", "grid"], capture_output=True, text=True, timeout=5,
            ).stdout.strip()
            if grid_which:
                search_paths.append(Path(grid_which).resolve().parent.parent / "parameters")
        except Exception:
            pass

        # 4. Common install paths
        search_paths.extend([
            Path("/opt/dock6/parameters"),
            Path("/usr/local/dock6/parameters"),
            Path.home() / "dock6" / "parameters",
        ])

        for search_dir in search_paths:
            candidate = search_dir / "vdw_AMBER_parm99.defn"
            if candidate.exists():
                vdw_defn_file = str(candidate.resolve())
                logger.info(f"    vdw_defn_file: {vdw_defn_file}")
                break

        if vdw_defn_file is None:
            logger.error("    vdw_AMBER_parm99.defn not found.")
            logger.error("    Set dock6_home in 03_configs/01b_grid_generation.yaml")
            logger.error(f"    Searched: {[str(p) for p in search_paths]}")
            return None

    # --- Create symlinks in output_dir for short paths ---
    # ALL DOCK6 programs (including grid) have a Fortran 80-char path limit.
    # We symlink all input files into the output directory and use filenames only.
    out = Path(output_dir)

    receptor_link = out / "receptor.mol2"
    box_file = out / "spheres_ligand_box.pdb"  # already here from step 4
    vdw_link = out / "vdw_AMBER_parm99.defn"

    # Symlink receptor mol2
    if not receptor_link.exists():
        receptor_link.symlink_to(Path(receptor_charged_mol2).resolve())
        logger.debug(f"    Symlinked: receptor.mol2 -> {receptor_charged_mol2}")

    # Symlink vdw definition file
    if not vdw_link.exists():
        vdw_link.symlink_to(Path(vdw_defn_file).resolve())
        logger.debug(f"    Symlinked: vdw_AMBER_parm99.defn -> {vdw_defn_file}")

    # grid.in content — ALL filenames only (no paths > 80 chars)
    grid_in_content = f"""\
compute_grids                  yes
grid_spacing                   {grid_spacing}
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         {energy_cutoff_distance}
atom_model                     a
attractive_exponent            {attractive_exponent}
repulsive_exponent             {repulsive_exponent}
distance_dielectric            yes
dielectric_factor              {dielectric_factor}
bump_filter                    yes
bump_overlap                   {bump_overlap}
allow_non_integral_charges     yes
receptor_file                  receptor.mol2
box_file                       spheres_ligand_box.pdb
vdw_definition_file            vdw_AMBER_parm99.defn
score_grid_prefix              ligand
"""

    grid_in_path = out / "grid.in"
    grid_in_path.write_text(grid_in_content)

    logger.info(f"  Step 5: Grid calculation (spacing={grid_spacing} A)")
    logger.info(f"    This may take 5-45 minutes...")

    try:
        result = subprocess.run(
            ["grid", "-i", "grid.in", "-o", "grid.out"],
            capture_output=True, text=True,
            timeout=3600,  # 1 hour max
            cwd=str(out),
        )

        nrg_path = out / "ligand.nrg"
        bmp_path = out / "ligand.bmp"

        if nrg_path.exists() and bmp_path.exists():
            logger.info(f"    -> {nrg_path.name} ({nrg_path.stat().st_size} bytes)")
            logger.info(f"    -> {bmp_path.name} ({bmp_path.stat().st_size} bytes)")
            return {
                "energy_grid": str(nrg_path),
                "bump_grid": str(bmp_path),
                "grid_in": str(grid_in_path),
            }
        else:
            logger.error("    grid calculation produced no output")
            grid_out = out / "grid.out"
            if grid_out.exists():
                logger.error(f"    grid.out tail: {grid_out.read_text()[-500:]}")
            if result.stderr:
                logger.error(f"    stderr: {result.stderr[:500]}")
            return None
    except FileNotFoundError:
        logger.error("    grid program not found in PATH")
        return None
    except subprocess.TimeoutExpired:
        logger.error("    grid calculation timed out (>1 hour)")
        return None


# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

def run_grid_generation(
        receptor_noH_pdb: Union[str, Path],
        receptor_charged_mol2: Union[str, Path],
        output_dir: Union[str, Path],
        binding_site_method: str = "reference_ligand",
        reference_mol2: Optional[str] = None,
        residues: Optional[List[str]] = None,
        center: Optional[List[float]] = None,
        receptor_pdb_for_residues: Optional[str] = None,
        chain: Optional[str] = None,
        radius: float = 10.0,
        probe_radius: float = 1.4,
        max_spheres: int = 50,
        box_margin: float = 1.0,
        grid_spacing: float = 0.3,
        energy_cutoff_distance: float = 9999.0,
        attractive_exponent: int = 6,
        repulsive_exponent: int = 12,
        dielectric_factor: int = 4,
        bump_overlap: float = 0.75,
        vdw_defn_file: Optional[str] = None,
        dock6_home: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run the complete DOCK6 grid generation pipeline.

    Args:
        receptor_noH_pdb: Receptor PDB without hydrogens (for DMS)
        receptor_charged_mol2: Receptor mol2 with charges (for grid)
        output_dir: Directory for all output files
        binding_site_method: "reference_ligand" | "residues" | "coordinates"
        reference_mol2: Path to reference ligand mol2/pdb (method A)
        residues: List of residue IDs, e.g. ["GLU529", "TRP555"] (method B)
        center: Explicit coordinates [x, y, z] (method C)
        receptor_pdb_for_residues: PDB to read residue coords from (method B)
        chain: Chain ID for residue selection
        radius: Radius for sphere selection (Angstroms)
        probe_radius: Probe radius for DMS surface
        max_spheres: Max spheres to keep (not currently enforced by sphgen)
        box_margin: Margin around spheres for grid box
        grid_spacing: Grid spacing in Angstroms
        energy_cutoff_distance: Energy cutoff distance
        attractive_exponent: VdW attractive exponent
        repulsive_exponent: VdW repulsive exponent
        dielectric_factor: Dielectric constant factor
        bump_overlap: Bump filter overlap threshold
        vdw_defn_file: Path to vdw_AMBER_parm99.defn (auto-detected if None)

    Returns:
        Dict with output paths and generation report
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    receptor_noH_pdb = str(receptor_noH_pdb)
    receptor_charged_mol2 = str(receptor_charged_mol2)

    logger.info("=" * 60)
    logger.info("  Grid Generation Pipeline")
    logger.info("=" * 60)
    logger.info(f"  Receptor (noH):    {Path(receptor_noH_pdb).name}")
    logger.info(f"  Receptor (mol2):   {Path(receptor_charged_mol2).name}")
    logger.info(f"  Binding site:      {binding_site_method}")
    logger.info(f"  Output:            {output_dir}")

    # --- Check tools ---
    tools = check_dock6_tools()
    missing = [t for t, ok in tools.items() if not ok]
    if missing:
        logger.error(f"  Missing DOCK6 tools: {missing}")
        return {"success": False, "error": f"Missing tools: {missing}"}

    report = {
        "start_time": datetime.now().isoformat(),
        "binding_site_method": binding_site_method,
        "steps": {},
    }

    # --- Step 1: DMS surface ---
    dms_path = generate_dms_surface(receptor_noH_pdb, str(output_dir), probe_radius)
    if not dms_path:
        return {"success": False, "error": "DMS surface generation failed", "report": report}
    report["steps"]["dms"] = dms_path

    # --- Step 2: Sphere generation ---
    sph_path = generate_spheres(dms_path, str(output_dir))
    if not sph_path:
        return {"success": False, "error": "Sphere generation failed", "report": report}
    report["steps"]["spheres"] = sph_path

    # --- Step 3: Sphere selection ---
    selected_path = str(output_dir / "spheres_ligand.sph")

    if binding_site_method == "reference_ligand":
        if not reference_mol2 or not Path(reference_mol2).exists():
            return {"success": False,
                    "error": f"Reference ligand not found: {reference_mol2}"}
        result = select_spheres_by_ligand(sph_path, reference_mol2, selected_path, radius)

    elif binding_site_method == "residues":
        if not residues:
            return {"success": False, "error": "No residues specified"}
        pdb_for_centroid = receptor_pdb_for_residues or receptor_noH_pdb
        try:
            centroid = compute_residue_centroid(pdb_for_centroid, residues, chain)
        except ValueError as e:
            return {"success": False, "error": str(e)}
        result = select_spheres_by_center(sph_path, centroid, selected_path, radius)
        report["binding_site_centroid"] = list(centroid)
        report["binding_site_residues"] = residues

    elif binding_site_method == "coordinates":
        if not center or len(center) != 3:
            return {"success": False, "error": f"Invalid center coordinates: {center}"}
        result = select_spheres_by_center(
            sph_path, tuple(center), selected_path, radius,
        )
        report["binding_site_center"] = center

    else:
        return {"success": False,
                "error": f"Unknown binding_site_method: {binding_site_method}"}

    if not result:
        return {"success": False, "error": "Sphere selection failed", "report": report}
    report["steps"]["selected_spheres"] = selected_path

    # --- Step 4: Box generation ---
    box_path = generate_box(selected_path, str(output_dir), box_margin)
    if not box_path:
        return {"success": False, "error": "Box generation failed", "report": report}
    report["steps"]["box"] = box_path

    # --- Step 5: Grid calculation ---
    grid_result = generate_grid(
        receptor_charged_mol2=receptor_charged_mol2,
        box_path=box_path,
        output_dir=str(output_dir),
        grid_spacing=grid_spacing,
        energy_cutoff_distance=energy_cutoff_distance,
        attractive_exponent=attractive_exponent,
        repulsive_exponent=repulsive_exponent,
        dielectric_factor=dielectric_factor,
        bump_overlap=bump_overlap,
        vdw_defn_file=vdw_defn_file,
        dock6_home=dock6_home,
    )
    if not grid_result:
        return {"success": False, "error": "Grid calculation failed", "report": report}
    report["steps"]["grid"] = grid_result

    # --- Save report ---
    report["end_time"] = datetime.now().isoformat()
    report["success"] = True

    report_path = output_dir / "grid_generation_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    # --- Summary ---
    logger.info("")
    logger.info(f"  Grid generation complete!")
    logger.info(f"  spheres_ligand.sph:     {selected_path}")
    logger.info(f"  spheres_ligand_box.pdb: {box_path}")
    logger.info(f"  ligand.nrg:             {grid_result['energy_grid']}")
    logger.info(f"  ligand.bmp:             {grid_result['bump_grid']}")

    return {
        "success": True,
        "spheres_file": selected_path,
        "energy_grid": grid_result["energy_grid"],
        "bump_grid": grid_result["bump_grid"],
        "box_file": box_path,
        "report": report,
        "report_path": str(report_path),
    }