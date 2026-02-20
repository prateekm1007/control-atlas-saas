"""
Unit tests for LAW-195 (Disulfide Geometry).

Tests verify that the disulfide bond measurement correctly identifies
CYS pairs with SG-SG distance < 3.0Å and measures deviation from
the ideal S-S bond length of 2.033Å.

Geometry notes:
  Ideal S-S: 2.033Å (Engh & Huber, 1991)
  Threshold: 0.20Å maximum deviation
  Detection cutoff: SG-SG < 3.0Å (putative disulfide)
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


class TestLaw195:
    def test_ideal_disulfide_pass(self, disulfide_pair):
        """Fixture has SG-SG exactly at 2.033Å → deviation 0.0 → PASS."""
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(disulfide_pair))
        res = results["LAW-195"]
        assert res["status"] == "PASS"
        assert res["observed"] <= 0.01  # Essentially zero deviation
        assert res["sample"] == 1

    def test_stretched_disulfide_veto(self):
        """SG-SG at 2.5Å → deviation 0.467Å → exceeds 0.20Å → VETO."""
        atoms = [
            _make_atom("N",  "N", [0.0, 0.0, 0.0],   res_name="CYS", res_seq=5),
            _make_atom("CA", "C", [1.47, 0.0, 0.0],   res_name="CYS", res_seq=5),
            _make_atom("CB", "C", [1.47, 1.53, 0.0],  res_name="CYS", res_seq=5),
            _make_atom("SG", "S", [0.0, 0.0, 0.0],    res_name="CYS", res_seq=5),
            _make_atom("N",  "N", [10.0, 0.0, 0.0],   res_name="CYS", res_seq=20),
            _make_atom("CA", "C", [11.47, 0.0, 0.0],  res_name="CYS", res_seq=20),
            _make_atom("CB", "C", [11.47, 1.53, 0.0], res_name="CYS", res_seq=20),
            _make_atom("SG", "S", [2.5, 0.0, 0.0],    res_name="CYS", res_seq=20),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-195"]
        assert res["status"] == "VETO"
        assert res["observed"] > 0.20
        assert res["sample"] == 1

    def test_no_cys_residues(self):
        """Structure with no CYS → sample=0, observed=0, PASS."""
        atoms = [
            _make_atom("N",  "N", [0.0, 0.0, 0.0],  res_name="ALA", res_seq=1),
            _make_atom("CA", "C", [1.47, 0.0, 0.0],  res_name="ALA", res_seq=1),
            _make_atom("C",  "C", [2.5, 1.0, 0.0],   res_name="ALA", res_seq=1),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-195"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0.0
        assert res["sample"] == 0

    def test_distant_cys_not_disulfide(self):
        """Two CYS with SG-SG > 3.0Å → not a disulfide pair → sample=0."""
        atoms = [
            _make_atom("N",  "N", [0.0, 0.0, 0.0],   res_name="CYS", res_seq=5),
            _make_atom("CA", "C", [1.47, 0.0, 0.0],   res_name="CYS", res_seq=5),
            _make_atom("CB", "C", [1.47, 1.53, 0.0],  res_name="CYS", res_seq=5),
            _make_atom("SG", "S", [0.0, 0.0, 0.0],    res_name="CYS", res_seq=5),
            _make_atom("N",  "N", [20.0, 0.0, 0.0],   res_name="CYS", res_seq=20),
            _make_atom("CA", "C", [21.47, 0.0, 0.0],  res_name="CYS", res_seq=20),
            _make_atom("CB", "C", [21.47, 1.53, 0.0], res_name="CYS", res_seq=20),
            _make_atom("SG", "S", [20.0, 0.0, 0.0],   res_name="CYS", res_seq=20),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-195"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0.0
        assert res["sample"] == 0

    def test_multiple_disulfides(self):
        """Two disulfide pairs: one ideal, one stretched → reports max deviation."""
        atoms = [
            # Pair 1: ideal (2.033Å)
            _make_atom("N",  "N", [0.0, 0.0, 0.0],     res_name="CYS", res_seq=1),
            _make_atom("CA", "C", [1.47, 0.0, 0.0],     res_name="CYS", res_seq=1),
            _make_atom("SG", "S", [0.0, 0.0, 0.0],      res_name="CYS", res_seq=1),
            _make_atom("N",  "N", [10.0, 0.0, 0.0],     res_name="CYS", res_seq=2),
            _make_atom("CA", "C", [11.47, 0.0, 0.0],    res_name="CYS", res_seq=2),
            _make_atom("SG", "S", [2.033, 0.0, 0.0],    res_name="CYS", res_seq=2),
            # Pair 2: stretched (2.4Å, deviation 0.367Å)
            _make_atom("N",  "N", [30.0, 0.0, 0.0],     res_name="CYS", res_seq=10),
            _make_atom("CA", "C", [31.47, 0.0, 0.0],    res_name="CYS", res_seq=10),
            _make_atom("SG", "S", [30.0, 0.0, 0.0],     res_name="CYS", res_seq=10),
            _make_atom("N",  "N", [40.0, 0.0, 0.0],     res_name="CYS", res_seq=11),
            _make_atom("CA", "C", [41.47, 0.0, 0.0],    res_name="CYS", res_seq=11),
            _make_atom("SG", "S", [32.4, 0.0, 0.0],     res_name="CYS", res_seq=11),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-195"]
        # Max deviation is from pair 2: |2.4 - 2.033| = 0.367
        assert res["status"] == "VETO"
        assert abs(res["observed"] - 0.367) < 0.01
        assert res["sample"] == 2

    def test_slightly_off_pass(self):
        """SG-SG at 2.15Å → deviation 0.117Å → under 0.20Å → PASS."""
        atoms = [
            _make_atom("N",  "N", [0.0, 0.0, 0.0],   res_name="CYS", res_seq=5),
            _make_atom("CA", "C", [1.47, 0.0, 0.0],   res_name="CYS", res_seq=5),
            _make_atom("SG", "S", [0.0, 0.0, 0.0],    res_name="CYS", res_seq=5),
            _make_atom("N",  "N", [10.0, 0.0, 0.0],   res_name="CYS", res_seq=20),
            _make_atom("CA", "C", [11.47, 0.0, 0.0],  res_name="CYS", res_seq=20),
            _make_atom("SG", "S", [2.15, 0.0, 0.0],   res_name="CYS", res_seq=20),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-195"]
        assert res["status"] == "PASS"
        assert abs(res["observed"] - 0.117) < 0.01
        assert res["sample"] == 1

    def test_uses_station_sop_threshold(self):
        """Verify threshold and ideal come from station_sop.py."""
        from tos.governance.station_sop import LAW_CANON
        assert LAW_CANON["LAW-195"]["threshold"] == 0.20
        assert LAW_CANON["LAW-195"]["unit"] == "Å"
        assert LAW_CANON["LAW-195"]["type"] == "Scalar"
        assert IDEAL_TABLE["S-S"] == 2.033

    def test_cys_without_sg_skipped(self):
        """CYS residue without SG atom → not included in search."""
        atoms = [
            _make_atom("N",  "N", [0.0, 0.0, 0.0],   res_name="CYS", res_seq=5),
            _make_atom("CA", "C", [1.47, 0.0, 0.0],   res_name="CYS", res_seq=5),
            # No SG atom
            _make_atom("N",  "N", [10.0, 0.0, 0.0],   res_name="CYS", res_seq=20),
            _make_atom("CA", "C", [11.47, 0.0, 0.0],  res_name="CYS", res_seq=20),
            _make_atom("SG", "S", [10.0, 2.0, 0.0],   res_name="CYS", res_seq=20),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-195"]
        assert res["status"] == "PASS"
        assert res["sample"] == 0
