"""
Unit tests for LAW-182 (Hydrophobic Burial).

Tests verify that hydrophobic core compactness is correctly measured
per chain, using the centroid-radius method.

Algorithm:
  1. For each chain, collect CA positions of hydrophobic residues
  2. Compute centroid of those positions
  3. Count fraction within 15Å radius
  4. Report minimum ratio across all chains

Hydrophobic residues: ALA, ILE, LEU, MET, PHE, PRO, TRP, VAL
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
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


class TestLaw182:
    def test_compact_hydrophobic_core_pass(self):
        """All hydrophobic CAs within 5Å of each other → ratio 1.0 → PASS."""
        atoms = [
            _make_atom("CA", "C", [0, 0, 0], res_name="LEU", res_seq=1),
            _make_atom("CA", "C", [2, 0, 0], res_name="ILE", res_seq=2),
            _make_atom("CA", "C", [0, 2, 0], res_name="VAL", res_seq=3),
            _make_atom("CA", "C", [2, 2, 0], res_name="PHE", res_seq=4),
            _make_atom("CA", "C", [1, 1, 0], res_name="ALA", res_seq=5),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-182"]
        assert res["status"] == "PASS"
        assert res["observed"] == 1.0
        assert res["sample"] == 5

    def test_scattered_hydrophobic_veto(self):
        """Hydrophobic CAs spread > 30Å apart → some outside 15Å radius → low ratio.
        With 5 residues spread along x from 0 to 100:
          centroid at x=50. Only the middle one (x=50) is within 15Å.
          Ratio = 1/5 = 0.2 < 0.3 → VETO.
        """
        atoms = [
            _make_atom("CA", "C", [0, 0, 0],   res_name="LEU", res_seq=1),
            _make_atom("CA", "C", [25, 0, 0],  res_name="ILE", res_seq=2),
            _make_atom("CA", "C", [50, 0, 0],  res_name="VAL", res_seq=3),
            _make_atom("CA", "C", [75, 0, 0],  res_name="PHE", res_seq=4),
            _make_atom("CA", "C", [100, 0, 0], res_name="ALA", res_seq=5),
        ]
        # centroid = (50, 0, 0)
        # distances: 50, 25, 0, 25, 50
        # within 15Å: only res 3 (d=0) → ratio = 1/5 = 0.2
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-182"]
        assert res["status"] == "VETO"
        assert res["observed"] < 0.3
        assert res["sample"] == 5

    def test_no_hydrophobic_residues(self):
        """Structure with only GLY/SER → no hydrophobic → ratio 1.0 → PASS."""
        atoms = [
            _make_atom("N",  "N", [0, 0, 0],  res_name="GLY", res_seq=1),
            _make_atom("CA", "C", [1.47, 0, 0], res_name="GLY", res_seq=1),
            _make_atom("C",  "C", [2.5, 1, 0], res_name="GLY", res_seq=1),
            _make_atom("N",  "N", [3.8, 1, 0], res_name="SER", res_seq=2),
            _make_atom("CA", "C", [5.27, 1, 0], res_name="SER", res_seq=2),
            _make_atom("C",  "C", [6.3, 2, 0], res_name="SER", res_seq=2),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-182"]
        assert res["status"] == "PASS"
        assert res["observed"] == 1.0
        assert res["sample"] == 0

    def test_multichain_worst_governs(self):
        """Chain A: compact (ratio 1.0). Chain B: scattered (ratio < 0.3).
        Reports minimum → VETO from chain B."""
        atoms = [
            # Chain A: compact cluster
            _make_atom("CA", "C", [0, 0, 0],  res_name="LEU", res_seq=1, chain_id="A"),
            _make_atom("CA", "C", [2, 0, 0],  res_name="ILE", res_seq=2, chain_id="A"),
            _make_atom("CA", "C", [1, 1, 0],  res_name="VAL", res_seq=3, chain_id="A"),
            # Chain B: scattered along x
            _make_atom("CA", "C", [0, 0, 0],   res_name="LEU", res_seq=1, chain_id="B"),
            _make_atom("CA", "C", [25, 0, 0],  res_name="ILE", res_seq=2, chain_id="B"),
            _make_atom("CA", "C", [50, 0, 0],  res_name="VAL", res_seq=3, chain_id="B"),
            _make_atom("CA", "C", [75, 0, 0],  res_name="PHE", res_seq=4, chain_id="B"),
            _make_atom("CA", "C", [100, 0, 0], res_name="ALA", res_seq=5, chain_id="B"),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-182"]
        assert res["status"] == "VETO"
        assert res["observed"] < 0.3
        assert res["sample"] == 8  # 3 from A + 5 from B

    def test_mixed_hydrophobic_and_polar(self):
        """Non-hydrophobic residues are ignored in the calculation."""
        atoms = [
            # Hydrophobic: compact
            _make_atom("CA", "C", [0, 0, 0], res_name="LEU", res_seq=1),
            _make_atom("CA", "C", [2, 0, 0], res_name="ILE", res_seq=2),
            _make_atom("CA", "C", [1, 1, 0], res_name="VAL", res_seq=3),
            # Polar: far away (should NOT affect hydrophobic burial)
            _make_atom("N",  "N", [100, 0, 0], res_name="GLY", res_seq=4),
            _make_atom("CA", "C", [101, 0, 0], res_name="GLY", res_seq=4),
            _make_atom("C",  "C", [102, 1, 0], res_name="GLY", res_seq=4),
            _make_atom("N",  "N", [103, 1, 0], res_name="SER", res_seq=5),
            _make_atom("CA", "C", [104, 1, 0], res_name="SER", res_seq=5),
            _make_atom("C",  "C", [105, 2, 0], res_name="SER", res_seq=5),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-182"]
        assert res["status"] == "PASS"
        assert res["observed"] == 1.0
        assert res["sample"] == 3  # Only hydrophobic counted

    def test_uses_station_sop_threshold(self):
        """Verify threshold and operator from station_sop.py."""
        from tos.governance.station_sop import LAW_CANON
        assert LAW_CANON["LAW-182"]["threshold"] == 0.3
        assert LAW_CANON["LAW-182"]["operator"] == ">="
        assert LAW_CANON["LAW-182"]["unit"] == "ratio"
        assert LAW_CANON["LAW-182"]["type"] == "Rate"

    def test_simple_backbone_pass(self, simple_backbone):
        """simple_backbone has 2 ALA residues (hydrophobic) close together → PASS."""
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(simple_backbone))
        res = results["LAW-182"]
        assert res["status"] == "PASS"
        assert res["observed"] >= 0.3
