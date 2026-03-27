"""Tests for reference_docking pipeline."""
import pytest


class TestLigandPreparation:
    def test_import(self):
        from reference_docking.m00_preparation.ligand_preparation import run_ligand_preparation
        assert callable(run_ligand_preparation)


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

    def test_rigid_template_no_secondary(self):
        from reference_docking.m01_docking.dock6_runner import DOCK6_RIGID_TEMPLATE
        assert "conformer_search_type" in DOCK6_RIGID_TEMPLATE
        assert "rigid" in DOCK6_RIGID_TEMPLATE
        assert "_secondary" not in DOCK6_RIGID_TEMPLATE

    def test_flex_template_no_secondary(self):
        from reference_docking.m01_docking.dock6_runner import DOCK6_FLEX_TEMPLATE
        assert "_secondary" not in DOCK6_FLEX_TEMPLATE


class TestFootprintRescore:
    def test_import(self):
        from reference_docking.m01_docking.footprint_rescore import run_footprint_rescore
        assert callable(run_footprint_rescore)

    def test_template_no_secondary(self):
        from reference_docking.m01_docking.footprint_rescore import FPS_RESCORE_TEMPLATE
        assert "footprint_similarity_score_primary" in FPS_RESCORE_TEMPLATE
        assert "_secondary" not in FPS_RESCORE_TEMPLATE
        assert "gbsa_hawkins" not in FPS_RESCORE_TEMPLATE


class TestGbsaRescore:
    def test_import(self):
        from reference_docking.m01_docking.gbsa_rescore import run_gbsa_rescore
        assert callable(run_gbsa_rescore)

    def test_template_primary(self):
        from reference_docking.m01_docking.gbsa_rescore import GBSA_RESCORE_TEMPLATE
        assert "gbsa_hawkins_score_primary               yes" in GBSA_RESCORE_TEMPLATE
        assert "orient_ligand                            no" in GBSA_RESCORE_TEMPLATE
        assert "grid_score_primary                       no" in GBSA_RESCORE_TEMPLATE


class TestScoreCollector:
    def test_import(self):
        from reference_docking.m01_docking.score_collector import run_score_collection
        assert callable(run_score_collection)