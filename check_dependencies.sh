#!/bin/bash
# =============================================================================
# reference_docking — Dependency Checker
# =============================================================================
# Verifica que todas las dependencias estan instaladas.
#
# Usage:
#   bash check_dependencies.sh
#
# =============================================================================

echo "============================================================"
echo "  MOLECULAR_DOCKING v2.0 — Dependency Check"
echo "============================================================"
echo ""

PASS=0; FAIL=0; WARN=0

check_cmd() {
    local name="$1" cmd="$2" req="$3"
    if eval "$cmd" > /dev/null 2>&1; then
        echo "  [OK]   $name"
        ((PASS++))
    elif [ "$req" = "required" ]; then
        echo "  [FAIL] $name: NOT FOUND"
        ((FAIL++))
    else
        echo "  [WARN] $name: not found (optional)"
        ((WARN++))
    fi
}

# --- Python ---
echo "--- Python ---"
check_cmd "python3"     "python3 --version"         "required"
echo ""

# --- Python Libraries ---
echo "--- Python Libraries ---"
for lib in rdkit pandas numpy openpyxl yaml tqdm; do
    mod=$lib
    if python3 -c "import $mod" 2>/dev/null; then
        echo "  [OK]   $lib"; ((PASS++))
    else
        echo "  [FAIL] $lib"; ((FAIL++))
    fi
done

# Optional Python libs
for lib in pdb2pqr openbabel; do
    if python3 -c "import $lib" 2>/dev/null; then
        echo "  [OK]   $lib"; ((PASS++))
    else
        echo "  [WARN] $lib (optional)"; ((WARN++))
    fi
done
echo ""

# --- DOCK6 ---
echo "--- DOCK6 ---"
check_cmd "dock6"            "which dock6"            "required"
check_cmd "grid"             "which grid"             "required"
check_cmd "sphgen"           "which sphgen"           "required"
check_cmd "sphere_selector"  "which sphere_selector"  "required"
check_cmd "showbox"          "which showbox"          "required"
check_cmd "dms"              "which dms"              "required"

# Check DOCK6 parameter files
if [ -f "/opt/dock6/parameters/vdw_AMBER_parm99.defn" ]; then
    echo "  [OK]   vdw_AMBER_parm99.defn (/opt/dock6/parameters/)"
    ((PASS++))
elif [ -n "$DOCK_HOME" ] && [ -f "$DOCK_HOME/parameters/vdw_AMBER_parm99.defn" ]; then
    echo "  [OK]   vdw_AMBER_parm99.defn (\$DOCK_HOME/parameters/)"
    ((PASS++))
else
    echo "  [WARN] vdw_AMBER_parm99.defn: not found (set DOCK_HOME or install to /opt/dock6/)"
    ((WARN++))
fi
echo ""

# --- ChimeraX ---
echo "--- ChimeraX (receptor preparation) ---"
if [ -x "/usr/bin/chimerax-daily" ]; then
    echo "  [OK]   chimerax-daily (/usr/bin/chimerax-daily)"
    ((PASS++))
elif which chimerax-daily > /dev/null 2>&1; then
    echo "  [OK]   chimerax-daily ($(which chimerax-daily))"
    ((PASS++))
elif which chimerax > /dev/null 2>&1; then
    echo "  [OK]   chimerax ($(which chimerax))"
    ((PASS++))
else
    echo "  [WARN] ChimeraX: not found (needed for 00b receptor preparation)"
    echo "         Install from: https://www.cgl.ucsf.edu/chimerax/download.html"
    ((WARN++))
fi
echo ""

# --- AmberTools ---
echo "--- AmberTools ---"
check_cmd "antechamber"  "which antechamber"  "optional"
check_cmd "parmchk2"     "which parmchk2"     "optional"
check_cmd "reduce"       "which reduce"       "optional"
echo ""

# --- OpenBabel ---
echo "--- OpenBabel ---"
check_cmd "obabel"  "obabel -V"  "required"
echo ""

# --- Summary ---
echo "============================================================"
echo "  SUMMARY: $PASS passed, $FAIL failed, $WARN warnings"
echo "============================================================"
if [ $FAIL -gt 0 ]; then
    echo "  Fix FAIL items before running pipeline."
    echo ""
    echo "  Quick fix:"
    echo "    conda activate reference_docking_env"
    echo "    conda env update -f environment.yaml --prune"
    echo "    pip install -e '.[dev]'"
    exit 1
fi
echo "  Ready!"
exit 0