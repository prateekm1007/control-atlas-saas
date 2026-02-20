"""
Unit tests for LAW-120 (Bond Angle RMSD).

Tests verify that the angle measurement engine correctly computes
RMSD of backbone angles against Engh-Huber ideals.
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Atom, Structure, ConfidenceProfile
from tos.governance.station_sop import IDEAL_TABLE


def _make_atom(atom_name, element, pos, res_name="ALA", res_seq=1,
               chain_id="A", insertion_code="", b_iso=90.0):
    """Helper to create Atom with b_iso=90 (passes pLDDT coverage filter)."""
    return Atom(
        atom_name=atom_name, element=element,
        pos=np.array(pos, dtype=float),
        res_name=res_name, res_seq=res_seq,
        chain_id=chain_id, insertion_code=insertion_code,
        b_iso=b_iso,
    )


class TestLaw120:
    def test_measures_angles_on_simple_backbone(self, simple_backbone):
        """simple_backbone produces non-zero angle samples."""
        struct = Structure(
            atoms=simple_backbone, audit_id="TEST",
            confidence=ConfidenceProfile([], False, "predicted")
        )
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        res = results["LAW-120"]
        assert res["sample"] > 0
        assert res["observed"] >= 0.0

    def test_sample_counts_angles(self, simple_backbone):
        """Verify sample reflects actual angle measurements taken.
        3 residues with N,CA,C each:
          - 3 intra-residue N-CA-C angles
          - 2 intra-residue O-C-N angles (res 1 and 3 have O; res 2 has O too)
          - 2 inter-residue CA-C-N angles
        Total depends on which atoms are present and which pairs are sequential.
        """
        struct = Structure(
            atoms=simple_backbone, audit_id="TEST",
            confidence=ConfidenceProfile([], False, "predicted")
        )
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        res = results["LAW-120"]
        # At minimum: 3 N-CA-C + 2 CA-C-N = 5 (O-C-N depends on atom presence)
        assert res["sample"] >= 5

    def test_veto_on_distorted_angles(self):
        """Create geometry with badly distorted angles — should VETO.
        N-CA-C nearly straight (~170°) vs ideal 111° = ~59° deviation.
        With multiple such angles, RMSD >> 10°.
        """
        atoms = [
            # Residue 1: N-CA-C nearly collinear (angle ~170°)
            _make_atom("N",  "N", [0.0, 0.0, 0.0], res_seq=1),
            _make_atom("CA", "C", [1.47, 0.0, 0.0], res_seq=1),
            _make_atom("C",  "C", [2.94, 0.1, 0.0], res_seq=1),
            _make_atom("O",  "O", [2.94, 0.1, 1.2], res_seq=1),
            # Residue 2: Also distorted
            _make_atom("N",  "N", [4.2, 0.1, 0.0], res_seq=2),
            _make_atom("CA", "C", [5.67, 0.1, 0.0], res_seq=2),
            _make_atom("C",  "C", [7.14, 0.2, 0.0], res_seq=2),
            _make_atom("O",  "O", [7.14, 0.2, 1.2], res_seq=2),
            # Residue 3: Also distorted
            _make_atom("N",  "N", [8.4, 0.2, 0.0], res_seq=3),
            _make_atom("CA", "C", [9.87, 0.2, 0.0], res_seq=3),
            _make_atom("C",  "C", [11.34, 0.3, 0.0], res_seq=3),
            _make_atom("O",  "O", [11.34, 0.3, 1.2], res_seq=3),
        ]
        struct = Structure(
            atoms=atoms, audit_id="DISTORTED",
            confidence=ConfidenceProfile([], False, "predicted")
        )
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        res = results["LAW-120"]
        assert res["sample"] > 0, f"No angles measured — check b_iso and coverage filter"
        # N-CA-C ~170° vs ideal 111° = ~59° deviation per angle
        # RMSD of multiple 59° deviations >> 10°
        assert res["observed"] > 10.0, f"Expected RMSD > 10°, got {res['observed']}"
        assert res["status"] == "VETO"

    def test_uses_station_sop_threshold(self):
        """Verify LAW-120 threshold comes from station_sop, not hardcoded."""
        from tos.governance.station_sop import LAW_CANON
        assert LAW_CANON["LAW-120"]["threshold"] == 18.0
        assert LAW_CANON["LAW-120"]["unit"] == "°"
        assert LAW_CANON["LAW-120"]["type"] == "RMSD"

    def test_empty_structure(self):
        """Empty structure should not crash."""
        struct = Structure(
            atoms=[], audit_id="EMPTY",
            confidence=ConfidenceProfile([], False, "predicted")
        )
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        assert results == {} or "LAW-120" not in results or results.get("LAW-120", {}).get("sample", 0) == 0
