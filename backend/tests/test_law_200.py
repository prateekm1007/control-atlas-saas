"""
Unit tests for LAW-200 (Packing Quality).

Tests verify that bounding box volume per atom is correctly computed.
Well-packed proteins have ~10-20 Å³/atom. Threshold is 300 Å³/atom.

Geometry:
  Volume = (xmax-xmin) * (ymax-ymin) * (zmax-zmin)
  vol_per_atom = volume / num_atoms
  Dimensions clamped to minimum 1.0Å to handle planar structures.
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
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


class TestLaw200:
    def test_compact_structure_pass(self):
        """10 atoms in a 5x5x5 box = 125ų / 10 = 12.5 ų/atom → PASS."""
        atoms = []
        for i in range(10):
            atoms.append(_make_atom(
                "CA", "C",
                [float(i % 5), float((i // 5) % 5), float(i // 25)],
                res_seq=i + 1
            ))
        # Manually verify: x range 0-4, y range 0-1, z range 0
        # With 1.0 clamp: dims = (4, 1, 1) → vol = 4. 4/10 = 0.4
        # Actually let me make a proper 3D spread
        atoms = [
            _make_atom("CA", "C", [0, 0, 0], res_seq=1),
            _make_atom("CA", "C", [5, 0, 0], res_seq=2),
            _make_atom("CA", "C", [0, 5, 0], res_seq=3),
            _make_atom("CA", "C", [5, 5, 0], res_seq=4),
            _make_atom("CA", "C", [0, 0, 5], res_seq=5),
            _make_atom("CA", "C", [5, 0, 5], res_seq=6),
            _make_atom("CA", "C", [0, 5, 5], res_seq=7),
            _make_atom("CA", "C", [5, 5, 5], res_seq=8),
            _make_atom("CA", "C", [2.5, 2.5, 2.5], res_seq=9),
            _make_atom("CA", "C", [1, 1, 1], res_seq=10),
        ]
        # bbox: 5x5x5 = 125. 125/10 = 12.5 Å³/atom
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-200"]
        assert res["status"] == "PASS"
        assert abs(res["observed"] - 12.5) < 0.1
        assert res["sample"] == 10

    def test_sparse_structure_veto(self):
        """2 atoms 100Å apart = 100x1x1 / 2 = 50 per atom.
        But let's make it really sparse: 2 atoms, bbox 1000x1x1."""
        atoms = [
            _make_atom("CA", "C", [0, 0, 0], res_seq=1),
            _make_atom("CA", "C", [1000, 0, 0], res_seq=2),
        ]
        # bbox: 1000 x 1(clamped) x 1(clamped) = 1000. 1000/2 = 500 Å³/atom
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-200"]
        assert res["status"] == "VETO"
        assert res["observed"] == 500.0
        assert res["sample"] == 2

    def test_simple_backbone_pass(self, simple_backbone):
        """simple_backbone has 14 atoms in a small box → well under 300."""
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(simple_backbone))
        res = results["LAW-200"]
        assert res["status"] == "PASS"
        assert res["observed"] < 300.0
        assert res["sample"] == 14

    def test_planar_structure_clamped(self):
        """All atoms in z=0 plane → z-dimension clamped to 1.0Å."""
        atoms = [
            _make_atom("CA", "C", [0, 0, 0], res_seq=1),
            _make_atom("CA", "C", [10, 0, 0], res_seq=2),
            _make_atom("CA", "C", [0, 10, 0], res_seq=3),
            _make_atom("CA", "C", [10, 10, 0], res_seq=4),
        ]
        # bbox: 10 x 10 x 1(clamped) = 100. 100/4 = 25 Å³/atom
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-200"]
        assert abs(res["observed"] - 25.0) < 0.1
        assert res["sample"] == 4

    def test_single_atom(self):
        """Single atom → bbox 1x1x1 (all clamped) = 1.0 Å³/atom."""
        atoms = [
            _make_atom("CA", "C", [5, 5, 5], res_seq=1),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-200"]
        assert res["status"] == "PASS"
        assert abs(res["observed"] - 1.0) < 0.1
        assert res["sample"] == 1

    def test_uses_station_sop_threshold(self):
        """Verify threshold comes from station_sop.py."""
        from tos.governance.station_sop import LAW_CANON
        assert LAW_CANON["LAW-200"]["threshold"] == 300.0
        assert LAW_CANON["LAW-200"]["unit"] == "Å³/atom"
        assert LAW_CANON["LAW-200"]["type"] == "Rate"

    def test_exactly_at_threshold(self):
        """Volume/atom exactly 300.0 → PASS (operator is <=)."""
        # Need: vol / n_atoms = 300
        # 3 atoms, bbox = 900. dims could be 30x30x1 = 900. 900/3 = 300.
        atoms = [
            _make_atom("CA", "C", [0, 0, 0], res_seq=1),
            _make_atom("CA", "C", [30, 0, 0], res_seq=2),
            _make_atom("CA", "C", [0, 30, 0], res_seq=3),
        ]
        # bbox: 30 x 30 x 1(clamped) = 900. 900/3 = 300.0
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-200"]
        assert res["status"] == "PASS"
        assert abs(res["observed"] - 300.0) < 0.1
