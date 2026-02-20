"""
Test fixtures for Toscanini physics engine tests.
Provides mock atoms and cached benchmark structures.

NOTE: b_iso is set to 90.0 by default because for predicted structures,
b_iso serves as the pLDDT confidence proxy. The coverage filter in
run_full_audit excludes residues with mean b_iso < 70 from the "core"
set. Using b_iso=20.0 would cause all residues to be excluded, making
LAW-100, LAW-120, LAW-125, and LAW-135 produce zero samples.

NOTE: CB positions use z=-1.0 to produce negative N-CA-C-CB improper
dihedrals, consistent with L-amino acid chirality (LAW-145).
Verified: dihedral_deg(N, CA, C, CB) < 0 for both residues.
"""
import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.ingestion.processor import Atom


@pytest.fixture
def mock_atom_factory():
    """Factory to create Atom objects at specified positions."""
    def _make(atom_name, element, pos, res_name="ALA", res_seq=1,
              chain_id="A", insertion_code="", b_iso=90.0):
        return Atom(
            atom_name=atom_name,
            element=element,
            pos=np.array(pos, dtype=float),
            res_name=res_name,
            res_seq=res_seq,
            chain_id=chain_id,
            insertion_code=insertion_code,
            b_iso=b_iso,
        )
    return _make


@pytest.fixture
def simple_backbone(mock_atom_factory):
    """
    A minimal 3-residue backbone with known geometry.

    Residue 1 (ALA): N at origin, CA at (1.47, 0, 0), C at (2.5, 1.0, 0)
    Residue 2 (GLY): N at (3.8, 1.0, 0), CA at (5.27, 1.0, 0), C at (6.3, 2.0, 0)
    Residue 3 (ALA): N at (7.6, 2.0, 0), CA at (9.07, 2.0, 0), C at (10.1, 3.0, 0)

    CA-CA distances: ~3.8-3.9 Angstroms (normal peptide spacing)
    b_iso=90.0 to pass predicted coverage filter (pLDDT >= 70)

    CB chirality: z=-1.0 gives L-amino acid chirality (negative N-CA-C-CB dihedral)
    Verified: res1 CB(1.47,-1,-1) → dihedral=-125.66°, res3 CB(9.07,1,-1) → dihedral=-125.66°
    """
    make = mock_atom_factory
    return [
        # Residue 1
        make("N",  "N", [0.0, 0.0, 0.0],    res_name="ALA", res_seq=1),
        make("CA", "C", [1.47, 0.0, 0.0],    res_name="ALA", res_seq=1),
        make("C",  "C", [2.5, 1.0, 0.0],     res_name="ALA", res_seq=1),
        make("O",  "O", [2.5, 1.0, 1.2],     res_name="ALA", res_seq=1),
        make("CB", "C", [1.47, -1.0, -1.0],  res_name="ALA", res_seq=1),
        # Residue 2 (GLY — no CB, no chirality)
        make("N",  "N", [3.8, 1.0, 0.0],     res_name="GLY", res_seq=2),
        make("CA", "C", [5.27, 1.0, 0.0],    res_name="GLY", res_seq=2),
        make("C",  "C", [6.3, 2.0, 0.0],     res_name="GLY", res_seq=2),
        make("O",  "O", [6.3, 2.0, 1.2],     res_name="GLY", res_seq=2),
        # Residue 3
        make("N",  "N", [7.6, 2.0, 0.0],     res_name="ALA", res_seq=3),
        make("CA", "C", [9.07, 2.0, 0.0],    res_name="ALA", res_seq=3),
        make("C",  "C", [10.1, 3.0, 0.0],    res_name="ALA", res_seq=3),
        make("O",  "O", [10.1, 3.0, 1.2],    res_name="ALA", res_seq=3),
        make("CB", "C", [9.07, 1.0, -1.0],   res_name="ALA", res_seq=3),
    ]


@pytest.fixture
def broken_backbone(mock_atom_factory):
    """
    A backbone with a chain break between residues 2 and 3.
    CA(2) to CA(3) distance is ~10 Angstroms (should trigger LAW-160).
    """
    make = mock_atom_factory
    return [
        make("N",  "N", [0.0, 0.0, 0.0],     res_name="ALA", res_seq=1),
        make("CA", "C", [1.47, 0.0, 0.0],     res_name="ALA", res_seq=1),
        make("C",  "C", [2.5, 1.0, 0.0],      res_name="ALA", res_seq=1),
        make("N",  "N", [3.8, 1.0, 0.0],      res_name="GLY", res_seq=2),
        make("CA", "C", [5.27, 1.0, 0.0],     res_name="GLY", res_seq=2),
        make("C",  "C", [6.3, 2.0, 0.0],      res_name="GLY", res_seq=2),
        # GAP: residue 3 is 10 Angstroms away
        make("N",  "N", [15.0, 2.0, 0.0],     res_name="ALA", res_seq=3),
        make("CA", "C", [16.47, 2.0, 0.0],    res_name="ALA", res_seq=3),
        make("C",  "C", [17.5, 3.0, 0.0],     res_name="ALA", res_seq=3),
    ]


@pytest.fixture
def multichain_atoms(mock_atom_factory):
    """Atoms from two chains: A and B."""
    make = mock_atom_factory
    return [
        make("CA", "C", [0.0, 0.0, 0.0], res_name="ALA", res_seq=1, chain_id="A"),
        make("CA", "C", [3.8, 0.0, 0.0], res_name="GLY", res_seq=2, chain_id="A"),
        make("CA", "C", [50.0, 0.0, 0.0], res_name="ALA", res_seq=1, chain_id="B"),
        make("CA", "C", [53.8, 0.0, 0.0], res_name="GLY", res_seq=2, chain_id="B"),
    ]


@pytest.fixture
def disulfide_pair(mock_atom_factory):
    """Two CYS residues with SG atoms at ideal disulfide distance (2.033 A)."""
    make = mock_atom_factory
    return [
        make("N",  "N", [0.0, 0.0, 0.0],    res_name="CYS", res_seq=5),
        make("CA", "C", [1.47, 0.0, 0.0],    res_name="CYS", res_seq=5),
        make("CB", "C", [1.47, 1.53, 0.0],   res_name="CYS", res_seq=5),
        make("SG", "S", [1.47, 3.33, 0.0],   res_name="CYS", res_seq=5),
        make("N",  "N", [10.0, 0.0, 0.0],    res_name="CYS", res_seq=20),
        make("CA", "C", [11.47, 0.0, 0.0],   res_name="CYS", res_seq=20),
        make("CB", "C", [11.47, 1.53, 0.0],  res_name="CYS", res_seq=20),
        make("SG", "S", [1.47, 5.363, 0.0],  res_name="CYS", res_seq=20),
    ]
