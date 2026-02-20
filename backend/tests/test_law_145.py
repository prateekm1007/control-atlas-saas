"""
Unit tests for LAW-145 (Chirality).

Tests verify that the Cα chirality check correctly identifies
D-amino acid violations via the improper dihedral N-CA-C-CB.

CRITICAL: N, CA, C must NOT be collinear. If they are, dihedral_deg
returns 0.0 (degenerate case) and chirality cannot be determined.
All test geometries use non-collinear backbone atoms.

Verified geometry (from coord_math.dihedral_deg):
  N(-0.5,0.8,0) CA(0,0,0) C(1.2,0.7,0):
    CB(0,-1,-0.5) → dihedral = -149.94° (L-amino acid)
    CB(0,-1,+0.5) → dihedral = +149.94° (D-amino acid)
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


def _l_residue(res_seq):
    """Create a single residue with L-amino acid chirality.
    N-CA-C are non-collinear. CB below plane → negative dihedral.
    """
    return [
        _make_atom("N",  "N", [-0.5, 0.8, 0.0], res_seq=res_seq),
        _make_atom("CA", "C", [0.0, 0.0, 0.0],  res_seq=res_seq),
        _make_atom("C",  "C", [1.2, 0.7, 0.0],  res_seq=res_seq),
        _make_atom("O",  "O", [1.5, 1.5, 0.0],  res_seq=res_seq),
        _make_atom("CB", "C", [0.0, -1.0, -0.5], res_seq=res_seq),
    ]


def _d_residue(res_seq):
    """Create a single residue with D-amino acid chirality.
    Same backbone as _l_residue but CB mirrored → positive dihedral.
    """
    return [
        _make_atom("N",  "N", [-0.5, 0.8, 0.0], res_seq=res_seq),
        _make_atom("CA", "C", [0.0, 0.0, 0.0],  res_seq=res_seq),
        _make_atom("C",  "C", [1.2, 0.7, 0.0],  res_seq=res_seq),
        _make_atom("O",  "O", [1.5, 1.5, 0.0],  res_seq=res_seq),
        _make_atom("CB", "C", [0.0, -1.0, 0.5],  res_seq=res_seq),
    ]


class TestLaw145:
    def test_l_amino_acid_pass(self):
        """L-amino acid chirality → improper dihedral < 0 → PASS."""
        atoms = _l_residue(1)

        # Verify geometry produces negative improper dihedral
        improper = dihedral_deg(
            np.array([-0.5, 0.8, 0.0]),
            np.array([0.0, 0.0, 0.0]),
            np.array([1.2, 0.7, 0.0]),
            np.array([0.0, -1.0, -0.5])
        )
        assert improper < 0.0, f"Expected negative improper, got {improper}"

        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-145"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0
        assert res["sample"] == 1

    def test_d_amino_acid_veto(self):
        """D-amino acid chirality → improper dihedral > 0 → VETO."""
        atoms = _d_residue(1)

        # Verify geometry produces positive improper dihedral
        improper = dihedral_deg(
            np.array([-0.5, 0.8, 0.0]),
            np.array([0.0, 0.0, 0.0]),
            np.array([1.2, 0.7, 0.0]),
            np.array([0.0, -1.0, 0.5])
        )
        assert improper > 0.0, f"Expected positive improper, got {improper}"

        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-145"]
        assert res["status"] == "VETO"
        assert res["observed"] == 1
        assert res["sample"] == 1

    def test_glycine_excluded(self):
        """GLY has no CB → should be skipped, sample = 0."""
        atoms = [
            _make_atom("N",  "N", [0.0, 0.0, 0.0],  res_name="GLY", res_seq=1),
            _make_atom("CA", "C", [1.47, 0.0, 0.0],  res_name="GLY", res_seq=1),
            _make_atom("C",  "C", [2.5, 1.0, 0.0],   res_name="GLY", res_seq=1),
            _make_atom("O",  "O", [2.5, 1.0, 1.2],   res_name="GLY", res_seq=1),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-145"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0
        assert res["sample"] == 0

    def test_mixed_l_and_d(self):
        """Two L-residues + one D-residue → 1 violation → VETO."""
        atoms = _l_residue(1) + _d_residue(2) + _l_residue(3)
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-145"]
        assert res["status"] == "VETO"
        assert res["observed"] == 1
        assert res["sample"] == 3

    def test_missing_cb_skipped(self):
        """Residue without CB atom → skipped (not counted as violation)."""
        atoms = [
            _make_atom("N",  "N", [-0.5, 0.8, 0.0], res_seq=1),
            _make_atom("CA", "C", [0.0, 0.0, 0.0],  res_seq=1),
            _make_atom("C",  "C", [1.2, 0.7, 0.0],  res_seq=1),
            _make_atom("O",  "O", [1.5, 1.5, 0.0],  res_seq=1),
            # No CB atom
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(atoms))
        res = results["LAW-145"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0
        assert res["sample"] == 0

    def test_uses_station_sop_threshold(self):
        """Verify threshold comes from station_sop.py."""
        from tos.governance.station_sop import LAW_CANON
        assert LAW_CANON["LAW-145"]["threshold"] == 0.0
        assert LAW_CANON["LAW-145"]["operator"] == "="
        assert LAW_CANON["LAW-145"]["unit"] == "count"
        assert LAW_CANON["LAW-145"]["type"] == "Count"

    def test_simple_backbone_pass(self, simple_backbone):
        """simple_backbone fixture has L-chirality CB positions → should PASS."""
        results, _, _ = Tier1Measurements.run_full_audit(_make_struct(simple_backbone))
        res = results["LAW-145"]
        assert res["status"] == "PASS"
        assert res["observed"] == 0
        # 2 non-GLY residues with CB (res 1 ALA and res 3 ALA)
        assert res["sample"] == 2

    def test_improper_dihedral_math_verification(self):
        """Directly verify the improper dihedral sign convention.

        L-chirality: CB below N-CA-C plane → negative.
        D-chirality: CB above N-CA-C plane → positive.
        Mirror symmetry: |L| == |D|.
        """
        n  = np.array([-0.5, 0.8, 0.0])
        ca = np.array([0.0, 0.0, 0.0])
        c  = np.array([1.2, 0.7, 0.0])

        # L: CB below plane (z < 0)
        cb_l = np.array([0.0, -1.0, -0.5])
        l_val = dihedral_deg(n, ca, c, cb_l)
        assert l_val < 0.0, f"L-chirality should be negative, got {l_val}"

        # D: CB above plane (z > 0)
        cb_d = np.array([0.0, -1.0, 0.5])
        d_val = dihedral_deg(n, ca, c, cb_d)
        assert d_val > 0.0, f"D-chirality should be positive, got {d_val}"

        # Mirror symmetry
        assert abs(abs(l_val) - abs(d_val)) < 0.1
