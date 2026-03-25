"""Tests for reference_docking pipeline."""
import pytest


class TestBindingSiteDefinition:
    def test_import(self):
        from reference_docking.m00_preparation.binding_site_definition import run_binding_site_definition
        assert callable(run_binding_site_definition)


class TestGridGeneration:
    def test_import(self):
        from reference_docking.m01_docking.grid_generation import validate_existing_grids
        assert callable(validate_existing_grids)

    def test_missing_grids(self, tmp_path):
        from reference_docking.m01_docking.grid_generation import validate_existing_grids
        assert validate_existing_grids(str(tmp_path), "s.sph", "g.nrg", "g.bmp") is False


class TestDock6Runner:
    def test_import(self):
        from reference_docking.m01_docking.dock6_runner import run_dock6_batch
        assert callable(run_dock6_batch)

    def test_resolve_grid_prefix(self):
        from reference_docking.m01_docking.dock6_runner import resolve_grid_prefix
        result = resolve_grid_prefix("/path/to/grids", "grid.nrg")
        assert result == "/path/to/grids/grid"

    def test_rigid_template_exists(self):
        from reference_docking.m01_docking.dock6_runner import DOCK6_RIGID_TEMPLATE
        assert "conformer_search_type" in DOCK6_RIGID_TEMPLATE
        assert "rigid" in DOCK6_RIGID_TEMPLATE


class TestScoreCollector:
    def test_import(self):
        from reference_docking.m01_docking.score_collector import run_score_collection
        assert callable(run_score_collection)
