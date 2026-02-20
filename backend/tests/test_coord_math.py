import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.coord_math import (
    distance, angle_deg, dihedral_deg,
    get_atom_pos, get_chains, get_sequential_residues,
    get_residue_atoms, CHI1_ATOMS, NO_CHI1_RESIDUES,
)

class TestDistance:
    def test_known_distance(self):
        a = np.array([0.0, 0.0, 0.0])
        b = np.array([3.0, 4.0, 0.0])
        assert abs(distance(a, b) - 5.0) < 1e-10

    def test_zero_distance(self):
        a = np.array([1.0, 2.0, 3.0])
        assert abs(distance(a, a)) < 1e-10

class TestAngle:
    def test_right_angle(self):
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 0.0, 0.0])
        c = np.array([0.0, 1.0, 0.0])
        assert abs(angle_deg(a, b, c) - 90.0) < 1e-6

class TestDihedral:
    def test_return_range(self):
        np.random.seed(42)
        for _ in range(100):
            pts = [np.random.randn(3) for _ in range(4)]
            dih = dihedral_deg(*pts)
            assert -180.0 <= dih <= 180.0

class TestGetAtomPos:
    def test_found(self, simple_backbone):
        pos = get_atom_pos(simple_backbone, "A", 1, "CA")
        assert pos is not None
        assert abs(pos[0] - 1.47) < 1e-6

class TestGetResidueAtoms:
    def test_returns_all_atoms(self, simple_backbone):
        atoms_dict = get_residue_atoms(simple_backbone, "A", 1)
        assert "N" in atoms_dict
        assert "CA" in atoms_dict
        assert "C" in atoms_dict
        assert "O" in atoms_dict
        assert "CB" in atoms_dict
        assert len(atoms_dict) == 5

class TestChi1Definitions:
    def test_all_standard_residues_covered(self):
        standard = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                     "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                     "THR", "TRP", "TYR", "VAL"}
        covered = set(CHI1_ATOMS.keys()) | NO_CHI1_RESIDUES
        assert standard == covered
