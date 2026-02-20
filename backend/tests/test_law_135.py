"""
Unit tests for LAW-135 (Omega Planarity).

Tests verify that the omega dihedral measurement correctly identifies
non-planar peptide bonds (deviation > 30° from 0° or 180°).

Test atom positions are analytically chosen with verified dihedral angles:
  Trans (180°): CA(i) and CA(i+1) on opposite sides of C-N bond
  Cis (0°):     CA(i) and CA(i+1) on same side of C-N bond
  Twisted (90°): CA(i+1) perpendicular to the CA(i)-C-N plane

Geometry verification (all analytically confirmed):
  dihedral((0,1.5,0), (1,0,0), (2,0,0), (2,-1,0))  = 180°  (trans)
  dihedral((0,1.5,0), (1,0,0), (2,0,0), (2,1,0))   =   0°  (cis)
  dihedral((0,1.5,0), (1,0,0), (2,0,0), (2,0,-1))   = -90°  (twisted)
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.engine.coord_math import dihedral_deg
from tos.ingestion.processor import Atom, Structure, ConfidenceProfile


def _make_atom(atom_name, element, pos, res_name="ALA", res_seq=1,
               chain_id="A", insertion_code="", b_iso=90.0):
    return Atom(
        atom_name=atom_name, element=element,
        pos=np.array(pos, dtype=float),
        res_name=res_name, res_seq=res_seq,
        chain_id=chain_id, insertion_code=insertion_code,
        b_iso=b_iso,
    )


def _make_struct(atoms):
    return Structure(
        atoms=atoms, audit_id="TEST",
        confidence=ConfidenceProfile([], False, "predicted")
    )


class TestLaw135:
    def test_trans_peptide_pass(self):
        """All omega angles ~180° (trans) → 0% outliers → PASS.

        Geometry: CA(i) above C-N line, CA(i+1) below → opposite sides → 180°
        """
        atoms = [
            # Res 1: CA above
            _make_atom("N",  "N", [0.0, 0.5, 0.0], res_seq=1),
            _make_atom("CA", "C", [0.0, 1.5, 0.0], res_seq=1),
            _make_atom("C",  "C", [1.0, 0.0, 0.0], res_seq=1),
            _make_atom("O",  "O", [1.0, -1.0, 0.0], res_seq=1),
            # Res 2: CA below (opposite side → trans)
            _make_atom("N",  "N", [2.0, 0.0, 0.0], res_seq=2),
            _make_atom("CA", "C", [2.0, -1.0, 0.0], res_seq=2),
            _make_atom("C",  "C", [3.0, 0.0, 0.0], res_seq=2),
            _make_atom("O",  "O", [3.0, 1.0, 0.0], res_seq=2),
            # Res 3: CA above (opposite from CA2 → trans)
            _make_atom("N",  "N", [4.0, 0.0, 0.0], res_seq=3),
            _make_atom("CA", "C", [4.0, 1.0, 0.0], res_seq=3),
            _make_atom("C",  "C", [5.0, 0.0, 0.0], res_seq=3),
            _make_atom("O",  "O", [5.0, -1.0, 0.0], res_seq=3),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-135"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0.0
        assert res["sample"] == 2

    def test_twisted_peptide_veto(self):
        """All omega angles ~90° (twisted) → 100% outliers → VETO.

        Geometry: CA(i+1) perpendicular to the CA(i)-C-N plane (z-offset).
        Verified: dev_from_trans = 90°, dev_from_cis = 90°, min = 90° > 30°.
        """
        atoms = [
            # Res 1: CA in xy-plane
            _make_atom("N",  "N", [0.0, 0.5, 0.0], res_seq=1),
            _make_atom("CA", "C", [0.0, 1.5, 0.0], res_seq=1),
            _make_atom("C",  "C", [1.0, 0.0, 0.0], res_seq=1),
            _make_atom("O",  "O", [1.0, -1.0, 0.0], res_seq=1),
            # Res 2: CA out of plane (z = -1 → omega = -90°)
            _make_atom("N",  "N", [2.0, 0.0, 0.0], res_seq=2),
            _make_atom("CA", "C", [2.0, 0.0, -1.0], res_seq=2),
            _make_atom("C",  "C", [3.0, 0.0, 0.0], res_seq=2),
            _make_atom("O",  "O", [3.0, 1.0, 0.0], res_seq=2),
            # Res 3: CA perpendicular to CA2 (y = -1 → omega = -90°)
            _make_atom("N",  "N", [4.0, 0.0, 0.0], res_seq=3),
            _make_atom("CA", "C", [4.0, -1.0, 0.0], res_seq=3),
            _make_atom("C",  "C", [5.0, 0.0, 0.0], res_seq=3),
            _make_atom("O",  "O", [5.0, -1.0, 0.0], res_seq=3),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-135"]
        assert res["status"] == "VETO"
        assert res["observed"] == 100.0
        assert res["sample"] == 2

    def test_cis_peptide_not_outlier(self):
        """Omega ~0° (cis) is planar — NOT an outlier.

        Geometry: CA(i) and CA(i+1) on same side of C-N bond → cis → 0°.
        Cis peptides are rare but structurally valid (especially cis-Pro).
        LAW-135 measures planarity, not cis/trans classification.
        """
        atoms = [
            # Res 1: CA above
            _make_atom("N",  "N", [0.0, 0.5, 0.0], res_seq=1),
            _make_atom("CA", "C", [0.0, 1.5, 0.0], res_seq=1),
            _make_atom("C",  "C", [1.0, 0.0, 0.0], res_seq=1),
            _make_atom("O",  "O", [1.0, -1.0, 0.0], res_seq=1),
            # Res 2: CA above (same side → cis)
            _make_atom("N",  "N", [2.0, 0.0, 0.0], res_seq=2),
            _make_atom("CA", "C", [2.0, 1.0, 0.0], res_seq=2),
            _make_atom("C",  "C", [3.0, 0.0, 0.0], res_seq=2),
            _make_atom("O",  "O", [3.0, -1.0, 0.0], res_seq=2),
            # Res 3: CA above (same side → cis)
            _make_atom("N",  "N", [4.0, 0.0, 0.0], res_seq=3),
            _make_atom("CA", "C", [4.0, 1.0, 0.0], res_seq=3),
            _make_atom("C",  "C", [5.0, 0.0, 0.0], res_seq=3),
            _make_atom("O",  "O", [5.0, -1.0, 0.0], res_seq=3),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-135"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0.0
        assert res["sample"] == 2

    def test_sample_equals_sequential_pairs(self, simple_backbone):
        """Sample should equal number of sequential pairs with required atoms."""
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(simple_backbone))
        res = results["LAW-135"]
        # simple_backbone: 3 residues in chain A = 2 sequential pairs
        assert res["sample"] == 2

    def test_uses_station_sop_threshold(self):
        """Verify threshold comes from station_sop.py, not hardcoded."""
        from tos.governance.station_sop import LAW_CANON
        assert LAW_CANON["LAW-135"]["threshold"] == 3.0
        assert LAW_CANON["LAW-135"]["unit"] == "%"
        assert LAW_CANON["LAW-135"]["type"] == "Percentage"

    def test_omega_math_verification(self):
        """Directly verify dihedral computation for all test geometries.

        This test catches any discrepancy between our analytical calculations
        and coord_math.dihedral_deg before the full-audit tests run.
        """
        # Trans: opposite sides of C-N bond → 180°
        omega_trans = dihedral_deg(
            np.array([0.0, 1.5, 0.0]),   # CA_i
            np.array([1.0, 0.0, 0.0]),   # C_i
            np.array([2.0, 0.0, 0.0]),   # N_i+1
            np.array([2.0, -1.0, 0.0])   # CA_i+1
        )
        assert abs(abs(omega_trans) - 180.0) < 0.1, f"Expected ~180°, got {omega_trans}"

        # Cis: same side of C-N bond → 0°
        omega_cis = dihedral_deg(
            np.array([0.0, 1.5, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
            np.array([2.0, 1.0, 0.0])
        )
        assert abs(omega_cis) < 0.1, f"Expected ~0°, got {omega_cis}"

        # Twisted: CA(i+1) out of plane → ±90°
        omega_twist = dihedral_deg(
            np.array([0.0, 1.5, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
            np.array([2.0, 0.0, -1.0])
        )
        assert abs(abs(omega_twist) - 90.0) < 0.1, f"Expected ~90°, got {omega_twist}"

    def test_simple_backbone_no_crash(self, simple_backbone):
        """Verify LAW-135 runs on simple_backbone without errors."""
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(simple_backbone))
        assert "LAW-135" in results
        assert results["LAW-135"]["observed"] >= 0.0
        assert results["LAW-135"]["sample"] >= 0
