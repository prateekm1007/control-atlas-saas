"""
Unit tests for LAW-155 (Voxel Occupancy).

Tests verify that the 2Å voxel grid correctly counts occupied voxels
and reports voxels per residue.

Algorithm:
  1. Compute bounding box minimum
  2. For each atom, compute voxel index: floor((pos - bbox_min) / 2.0)
  3. Count unique voxel indices
  4. Report: occupied_voxels / num_residues

Threshold: 2.0 voxels/residue (operator >=)
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


class TestLaw155:
    def test_well_resolved_pass(self):
        """5 atoms in different 2Å voxels, 1 residue → 5 voxels/res → PASS."""
        atoms = [
            _make_atom("N",  "N", [0, 0, 0], res_seq=1),
            _make_atom("CA", "C", [3, 0, 0], res_seq=1),   # different voxel (3/2=1)
            _make_atom("C",  "C", [0, 3, 0], res_seq=1),   # different voxel
            _make_atom("O",  "O", [0, 0, 3], res_seq=1),   # different voxel
            _make_atom("CB", "C", [3, 3, 0], res_seq=1),   # different voxel
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-155"]
        assert res["status"] == "PASS"
        assert res["observed"] == 5.0
        assert res["sample"] == 5

    def test_sparse_structure_veto(self):
        """1 atom per residue, all in same voxel → 1 voxel / 3 residues < 2.0 → VETO.
        Actually: atoms in same 2Å cube, 3 residues → 1/3 = 0.33.
        """
        atoms = [
            _make_atom("CA", "C", [0.0, 0.0, 0.0], res_seq=1),
            _make_atom("CA", "C", [0.5, 0.5, 0.5], res_seq=2),  # same voxel
            _make_atom("CA", "C", [1.0, 1.0, 1.0], res_seq=3),  # same voxel (1/2=0)
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-155"]
        assert res["status"] == "VETO"
        assert res["observed"] < 2.0
        assert res["sample"] == 3

    def test_simple_backbone_pass(self, simple_backbone):
        """simple_backbone: 14 atoms spread across ~10Å, 3 residues.
        Many different voxels → high voxels/residue → PASS."""
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(simple_backbone))
        res = results["LAW-155"]
        assert res["status"] == "PASS"
        assert res["observed"] >= 2.0
        assert res["sample"] == 14

    def test_atoms_in_same_voxel_counted_once(self):
        """Multiple atoms in the same 2Å cube count as 1 voxel."""
        atoms = [
            _make_atom("N",  "N", [0.0, 0.0, 0.0], res_seq=1),
            _make_atom("CA", "C", [0.1, 0.1, 0.1], res_seq=1),  # same voxel
            _make_atom("C",  "C", [0.2, 0.2, 0.2], res_seq=1),  # same voxel
            _make_atom("O",  "O", [0.3, 0.3, 0.3], res_seq=1),  # same voxel
            _make_atom("CB", "C", [0.4, 0.4, 0.4], res_seq=1),  # same voxel
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-155"]
        # All 5 atoms in voxel (0,0,0) → 1 voxel / 1 residue = 1.0
        assert res["observed"] == 1.0
        assert res["status"] == "VETO"  # 1.0 < 2.0

    def test_spread_atoms_many_voxels(self):
        """Atoms spread across many voxels in 1 residue → high ratio."""
        atoms = [
            _make_atom("N",  "N", [0, 0, 0], res_seq=1),
            _make_atom("CA", "C", [3, 0, 0], res_seq=1),
            _make_atom("C",  "C", [6, 0, 0], res_seq=1),
            _make_atom("O",  "O", [9, 0, 0], res_seq=1),
            _make_atom("CB", "C", [12, 0, 0], res_seq=1),
            _make_atom("CG", "C", [15, 0, 0], res_seq=1),
        ]
        # Voxels: 0, 1, 3, 4, 6, 7 → 6 unique voxels / 1 residue = 6.0
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-155"]
        assert res["status"] == "PASS"
        assert res["observed"] >= 2.0

    def test_exactly_at_threshold(self):
        """Exactly 2.0 voxels/residue → PASS (operator >=)."""
        # 2 residues, 4 atoms in 4 different voxels → 4/2 = 2.0
        atoms = [
            _make_atom("N",  "N", [0, 0, 0], res_seq=1),
            _make_atom("CA", "C", [3, 0, 0], res_seq=1),
            _make_atom("N",  "N", [6, 0, 0], res_seq=2),
            _make_atom("CA", "C", [9, 0, 0], res_seq=2),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-155"]
        assert res["status"] == "PASS"
        assert abs(res["observed"] - 2.0) < 0.01

    def test_uses_station_sop_threshold(self):
        """Verify threshold and operator from station_sop.py."""
        from tos.governance.station_sop import LAW_CANON
        assert LAW_CANON["LAW-155"]["threshold"] == 2.0
        assert LAW_CANON["LAW-155"]["operator"] == ">="
        assert LAW_CANON["LAW-155"]["unit"] == "V"
        assert LAW_CANON["LAW-155"]["type"] == "Scalar"
