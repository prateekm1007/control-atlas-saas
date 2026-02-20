"""
Regression tests for canonical math migration (Module 1.2.1).

Verifies that migrating from local _dist/_angle/_dihedral to
coord_math.distance/angle_deg/dihedral_deg produces identical
results for all affected laws.
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.coord_math import distance, angle_deg, dihedral_deg
from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Atom, Structure, ConfidenceProfile


# ─────────────────────────────────────────────────────────────
# Part 1: Verify coord_math matches old local implementations
# ─────────────────────────────────────────────────────────────

class TestDistanceMigration:
    """Verify coord_math.distance matches old _dist."""

    def test_peptide_bond_length(self):
        c_pos = np.array([2.5, 1.0, 0.0])
        n_pos = np.array([3.8, 1.0, 0.0])
        assert abs(distance(c_pos, n_pos) - 1.3) < 0.01

    def test_ca_ca_spacing(self):
        ca1 = np.array([1.47, 0.0, 0.0])
        ca2 = np.array([5.27, 1.0, 0.0])
        d = distance(ca1, ca2)
        assert 3.8 < d < 4.0

    def test_zero_distance(self):
        p = np.array([1.0, 2.0, 3.0])
        assert distance(p, p) == 0.0


class TestAngleMigration:
    """Verify coord_math.angle_deg matches old _angle behavior."""

    def test_right_angle(self):
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 0.0, 0.0])
        c = np.array([0.0, 1.0, 0.0])
        assert abs(angle_deg(a, b, c) - 90.0) < 1e-6

    def test_straight_angle(self):
        a = np.array([0.0, 0.0, 0.0])
        b = np.array([1.0, 0.0, 0.0])
        c = np.array([2.0, 0.0, 0.0])
        assert abs(angle_deg(a, b, c) - 180.0) < 1e-6

    def test_backbone_nca_c_angle(self):
        """
        The fixture's synthetic geometry gives N-CA-C = ~136°.
        This is NOT ideal protein geometry (~111°) but the test verifies
        that angle_deg computes correctly, not that the fixture is ideal.
        """
        n = np.array([0.0, 0.0, 0.0])
        ca = np.array([1.47, 0.0, 0.0])
        c = np.array([2.5, 1.0, 0.0])
        ang = angle_deg(n, ca, c)
        # Verify it returns a valid angle (the exact value is ~135.8°)
        assert 130.0 < ang < 140.0

    def test_degenerate_returns_zero(self):
        p = np.array([1.0, 1.0, 1.0])
        assert angle_deg(p, p, p) == 0.0


class TestDihedralMigration:
    """Verify coord_math.dihedral_deg matches old _dihedral behavior."""

    def test_trans_configuration(self):
        a = np.array([1.0, 1.0, 0.0])
        b = np.array([0.0, 0.0, 0.0])
        c = np.array([1.0, 0.0, 0.0])
        d = np.array([1.0, -1.0, 0.0])
        dih = dihedral_deg(a, b, c, d)
        assert abs(abs(dih) - 180.0) < 1.0

    def test_cis_configuration(self):
        a = np.array([0.0, 1.0, 0.0])
        b = np.array([0.0, 0.0, 0.0])
        c = np.array([1.0, 0.0, 0.0])
        d = np.array([1.0, 1.0, 0.0])
        dih = dihedral_deg(a, b, c, d)
        assert abs(dih) < 1.0

    def test_range(self):
        np.random.seed(42)
        for _ in range(50):
            pts = [np.random.randn(3) for _ in range(4)]
            dih = dihedral_deg(*pts)
            assert -180.0 <= dih <= 180.0

    def test_degenerate_collinear(self):
        a = np.array([0.0, 0.0, 0.0])
        b = np.array([1.0, 0.0, 0.0])
        c = np.array([2.0, 0.0, 0.0])
        d = np.array([3.0, 1.0, 0.0])
        dih = dihedral_deg(a, b, c, d)
        assert dih == 0.0


# ─────────────────────────────────────────────────────────────
# Part 2: Full audit regression — verify law outputs
# ─────────────────────────────────────────────────────────────

class TestFullAuditRegression:
    """Run full audit on fixture structures and verify key laws."""

    def test_simple_backbone_all_pass(self, simple_backbone):
        """With b_iso=90.0, core is populated and all laws compute."""
        struct = Structure(
            atoms=simple_backbone, audit_id="REGTEST",
            confidence=ConfidenceProfile([], False, "predicted")
        )
        results, coverage, _ = Tier1Measurements.run_full_audit(struct)

        # LAW-100: With b_iso=90 (above pLDDT 70 threshold), core is populated
        assert "LAW-100" in results
        assert results["LAW-100"]["sample"] > 0

        # LAW-110: No backbone gaps in simple_backbone
        assert results["LAW-110"]["observed"] == 0
        assert results["LAW-110"]["status"] == "PASS"

        # LAW-160: CA-CA within threshold
        assert results["LAW-160"]["status"] == "PASS"
        assert 3.7 <= results["LAW-160"]["observed"] <= 4.0

        # LAW-120: Should have real angle data now
        assert results["LAW-120"]["sample"] > 0
        assert results["LAW-120"]["observed"] >= 0.0

        # All 15 laws present
        for lid in ["LAW-100", "LAW-105", "LAW-110", "LAW-120", "LAW-125",
                     "LAW-130", "LAW-135", "LAW-145", "LAW-150", "LAW-155",
                     "LAW-160", "LAW-170", "LAW-182", "LAW-195", "LAW-200"]:
            assert lid in results, f"{lid} missing from results"

    def test_broken_backbone_veto(self, broken_backbone):
        struct = Structure(
            atoms=broken_backbone, audit_id="REGTEST",
            confidence=ConfidenceProfile([], False, "predicted")
        )
        results, coverage, _ = Tier1Measurements.run_full_audit(struct)
        assert results["LAW-160"]["status"] == "VETO"
        assert results["LAW-160"]["observed"] > 9.0

    def test_law_160_sample_count_is_pairs(self, simple_backbone):
        """Verify LAW-160 reports CA pairs checked, not total residues."""
        struct = Structure(
            atoms=simple_backbone, audit_id="REGTEST",
            confidence=ConfidenceProfile([], False, "predicted")
        )
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        assert results["LAW-160"]["sample"] == 2

    def test_law_155_present(self, simple_backbone):
        """Verify LAW-155 is now explicitly in results."""
        struct = Structure(
            atoms=simple_backbone, audit_id="REGTEST",
            confidence=ConfidenceProfile([], False, "predicted")
        )
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        assert "LAW-155" in results
        assert results["LAW-155"]["status"] == "PASS"

    def test_no_local_math_functions(self):
        """Verify the local _dist, _angle, _dihedral are truly gone."""
        import inspect
        import tos.engine.tier1_measurements as mod
        source = inspect.getsource(mod)
        assert "def _dist(" not in source, "_dist still exists"
        assert "def _angle(" not in source, "_angle still exists"
        assert "def _dihedral(" not in source, "_dihedral still exists"
