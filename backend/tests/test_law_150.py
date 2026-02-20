"""
Unit tests for LAW-150 (Rotamer Audit).

Tests verify that chi1 dihedral angles are correctly measured and
classified against the three canonical rotamer wells:
  gauche+ (g+):  60° ± 30°  → [30°, 90°]
  gauche- (g-): -60° ± 30°  → [-90°, -30°]
  trans   (t):  180° ± 30°  → |chi1| >= 150°

Outlier: chi1 outside all three wells.

Geometry construction (verified empirically):
  Backbone: N(-0.5,0.8,0), CA(0,0,0), CB(0,-1.53,0)
  CG at distance 1.5 from CB, rotated in xz-plane.
  Formula: rotation = 270 - chi1_target (degrees)
  CG = CB + [sin(rotation)*1.5, 0, cos(rotation)*1.5]

  Verified mapping (all errors = 0.0°):
    target=  60° → chi1=  60°  (g+)
    target= -60° → chi1= -60°  (g-)
    target= 180° → chi1=-180°  (trans, equivalent)
    target= 120° → chi1= 120°  (outlier)
    target=   0° → chi1=   0°  (outlier)
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.engine.coord_math import dihedral_deg, CHI1_ATOMS, NO_CHI1_RESIDUES
from tos.ingestion.processor import Atom, Structure, ConfidenceProfile


def _make_atom(atom_name, element, pos, res_name="LEU", res_seq=1,
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


def _leu_with_chi1(chi1_target_deg, res_seq=1):
    """Create a LEU residue with a specific chi1 angle.

    Backbone: N(-0.5, 0.8, 0), CA(0, 0, 0) — non-collinear
    CB at (0, -1.53, 0) — along -y axis from CA
    CG position computed via: rotation = 270 - chi1_target
    CG = CB + [sin(rotation)*1.5, 0, cos(rotation)*1.5]

    Verified: all targets produce exact chi1 values (error = 0.0°).
    """
    rotation = 270.0 - chi1_target_deg
    rad = np.radians(rotation)
    cb = np.array([0.0, -1.53, 0.0])
    cg = cb + np.array([np.sin(rad) * 1.5, 0.0, np.cos(rad) * 1.5])

    return [
        _make_atom("N",  "N", [-0.5, 0.8, 0.0], res_name="LEU", res_seq=res_seq),
        _make_atom("CA", "C", [0.0, 0.0, 0.0],  res_name="LEU", res_seq=res_seq),
        _make_atom("C",  "C", [1.2, 0.7, 0.0],  res_name="LEU", res_seq=res_seq),
        _make_atom("O",  "O", [1.5, 1.5, 0.0],  res_name="LEU", res_seq=res_seq),
        _make_atom("CB", "C", cb.tolist(),        res_name="LEU", res_seq=res_seq),
        _make_atom("CG", "C", cg.tolist(),        res_name="LEU", res_seq=res_seq),
    ]


class TestLaw150:
    def test_verify_chi1_geometry(self):
        """Verify that _leu_with_chi1 produces the intended chi1 angles."""
        for target in [60.0, -60.0, 180.0, 120.0, 0.0]:
            atoms = _leu_with_chi1(target)
            n_pos = atoms[0].pos
            ca_pos = atoms[1].pos
            cb_pos = atoms[4].pos
            cg_pos = atoms[5].pos
            measured = dihedral_deg(n_pos, ca_pos, cb_pos, cg_pos)
            diff = abs(measured - target)
            if diff > 180:
                diff = 360 - diff
            assert diff < 1.0, f"Target {target}°, got {measured}°"

    def test_trans_rotamer_pass(self):
        """Chi1 ≈ 180° (trans) → in well → not outlier → PASS."""
        atoms = _leu_with_chi1(180.0)
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0.0
        assert res["sample"] == 1

    def test_gauche_plus_pass(self):
        """Chi1 ≈ 60° (g+) → in well → not outlier → PASS."""
        atoms = _leu_with_chi1(60.0)
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0.0
        assert res["sample"] == 1

    def test_gauche_minus_pass(self):
        """Chi1 ≈ -60° (g-) → in well → not outlier → PASS."""
        atoms = _leu_with_chi1(-60.0)
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0.0
        assert res["sample"] == 1

    def test_outlier_chi1_detected(self):
        """Chi1 ≈ 120° → outside all wells → outlier → 100% → VETO."""
        atoms = _leu_with_chi1(120.0)
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["status"] == "VETO"
        assert res["observed"] == 100.0
        assert res["sample"] == 1

    def test_zero_chi1_is_outlier(self):
        """Chi1 ≈ 0° → outside all wells → outlier."""
        atoms = _leu_with_chi1(0.0)
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["observed"] == 100.0
        assert res["sample"] == 1

    def test_gly_ala_excluded(self):
        """GLY and ALA have no chi1 → excluded from audit."""
        atoms = [
            _make_atom("N",  "N", [-0.5, 0.8, 0.0], res_name="GLY", res_seq=1),
            _make_atom("CA", "C", [0.0, 0.0, 0.0],  res_name="GLY", res_seq=1),
            _make_atom("C",  "C", [1.2, 0.7, 0.0],  res_name="GLY", res_seq=1),
            _make_atom("N",  "N", [3.0, 0.8, 0.0],  res_name="ALA", res_seq=2),
            _make_atom("CA", "C", [3.5, 0.0, 0.0],  res_name="ALA", res_seq=2),
            _make_atom("C",  "C", [4.7, 0.7, 0.0],  res_name="ALA", res_seq=2),
            _make_atom("CB", "C", [3.5, -1.0, -1.0], res_name="ALA", res_seq=2),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["sample"] == 0
        assert res["status"] == "PASS"

    def test_missing_chi1_atoms_skipped(self):
        """LEU without CG → cannot compute chi1 → skipped."""
        atoms = [
            _make_atom("N",  "N", [-0.5, 0.8, 0.0], res_name="LEU", res_seq=1),
            _make_atom("CA", "C", [0.0, 0.0, 0.0],  res_name="LEU", res_seq=1),
            _make_atom("C",  "C", [1.2, 0.7, 0.0],  res_name="LEU", res_seq=1),
            _make_atom("CB", "C", [0.0, -1.53, 0.0], res_name="LEU", res_seq=1),
            # No CG atom
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["sample"] == 0
        assert res["status"] == "PASS"

    def test_mixed_good_and_bad(self):
        """3 residues: 2 in wells + 1 outlier → 33.3% → VETO (> 20%)."""
        atoms = (
            _leu_with_chi1(60.0, res_seq=1) +
            _leu_with_chi1(-60.0, res_seq=2) +
            _leu_with_chi1(120.0, res_seq=3)
        )
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["status"] == "VETO"
        assert abs(res["observed"] - 33.33) < 1.0
        assert res["sample"] == 3

    def test_under_threshold_pass(self):
        """5 residues: 4 in wells + 1 outlier → 20% → PASS (exactly at threshold)."""
        atoms = (
            _leu_with_chi1(60.0, res_seq=1) +
            _leu_with_chi1(-60.0, res_seq=2) +
            _leu_with_chi1(180.0, res_seq=3) +
            _leu_with_chi1(60.0, res_seq=4) +
            _leu_with_chi1(120.0, res_seq=5)
        )
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-150"]
        assert res["status"] == "PASS"
        assert res["observed"] == 20.0
        assert res["sample"] == 5

    def test_uses_station_sop_threshold(self):
        """Verify threshold comes from station_sop.py."""
        from tos.governance.station_sop import LAW_CANON
        assert LAW_CANON["LAW-150"]["threshold"] == 20.0
        assert LAW_CANON["LAW-150"]["unit"] == "%"
        assert LAW_CANON["LAW-150"]["type"] == "Percentage"

    def test_chi1_atoms_table_coverage(self):
        """Verify all standard residues are covered by CHI1_ATOMS + NO_CHI1."""
        from tos.governance.station_sop import STANDARD_RESIDUES
        covered = set(CHI1_ATOMS.keys()) | NO_CHI1_RESIDUES
        assert covered == STANDARD_RESIDUES
