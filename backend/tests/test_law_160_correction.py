"""
WP-SCI-1 regression tests for LAW-160 sequence adjacency (directive §3).

The historical LAW-160 measured CA-CA between CONSECUTIVE ENTRIES of a
sorted residue list — a missing residue span (e.g. 45 -> 67) became a
physical chain-break veto. Definition v2 measures only immediately
sequential pairs (coord_math.residues_are_adjacent): insertion codes
honored, chains never joined, missing spans skipped, genuine coordinate
discontinuities between adjacent-numbered residues still veto.
"""
import copy
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Atom, Structure, ConfidenceProfile


def _atom(name, pos, res_seq, chain="A", icode="", res_name="ALA", b_iso=90.0):
    return Atom(atom_name=name, element="C", pos=np.array(pos, dtype=float),
                res_name=res_name, res_seq=res_seq, chain_id=chain,
                insertion_code=icode, b_iso=b_iso)


def _ca_residue(seq, x, chain="A", icode="", with_ca=True):
    """A residue with backbone atoms; CA placed at (x, 0, 0)."""
    atoms = [
        _atom("N", [x, -1.0, 0.0], seq, chain, icode),
        _atom("C", [x, 1.0, 0.0], seq, chain, icode),
        _atom("O", [x, 1.0, 1.2], seq, chain, icode),
    ]
    if with_ca:
        atoms.insert(1, _atom("CA", [x, 0.0, 0.0], seq, chain, icode))
    return atoms


def _struct(atoms):
    return Structure(atoms=atoms, audit_id="WPSCI1-160",
                     confidence=ConfidenceProfile([], False, "predicted"))


def _audit(atoms):
    results, _, _ = Tier1Measurements.run_full_audit(_struct(atoms))
    return results["LAW-160"]


class TestLaw160Adjacency:

    def test_ordinary_adjacent_residues_measured(self):
        atoms = _ca_residue(1, 0.0) + _ca_residue(2, 3.8) + _ca_residue(3, 7.6)
        res = _audit(atoms)
        assert res["status"] == "PASS"
        assert res["sample"] == 2
        assert res["observed"] == pytest.approx(3.8, abs=0.01)

    def test_missing_loop_is_not_a_chain_break(self):
        """Sequence 1-2-3 then 45-46-47: the 3->45 gap must be skipped,
        not measured as a ~160 A 'chain break'."""
        atoms = (_ca_residue(1, 0.0) + _ca_residue(2, 3.8) +
                 _ca_residue(3, 7.6) +
                 _ca_residue(45, 170.0) + _ca_residue(46, 173.8) +
                 _ca_residue(47, 177.6))
        res = _audit(atoms)
        assert res["status"] == "PASS"
        assert res["sample"] == 4
        assert res["observed"] == pytest.approx(3.8, abs=0.01)

    def test_huge_numbering_gap_with_far_ca_not_a_veto(self):
        """Even with the gap-flanking CA atoms 40 A apart, a numbering
        gap is absent density, not a physical break."""
        atoms = _ca_residue(45, 0.0) + _ca_residue(67, 40.0)
        res = _audit(atoms)
        assert res["status"] == "NOT_MEASURED"
        assert res["sample"] == 0
        assert res["observed"] is None

    def test_insertion_code_chain(self):
        """Residues 5, 5A, 6 form an adjacency chain of two pairs."""
        atoms = (_ca_residue(5, 0.0, icode="") +
                 _ca_residue(5, 3.8, icode="A") +
                 _ca_residue(6, 7.6, icode=""))
        res = _audit(atoms)
        assert res["status"] == "PASS"
        assert res["sample"] == 2

    def test_multiple_chains_never_joined(self):
        """Chains A and B each measured internally; the inter-chain
        distance is irrelevant."""
        atoms = (_ca_residue(1, 0.0, chain="A") + _ca_residue(2, 3.8, chain="A") +
                 _ca_residue(1, 60.0, chain="B") + _ca_residue(2, 63.8, chain="B"))
        res = _audit(atoms)
        assert res["status"] == "PASS"
        assert res["sample"] == 2  # one pair per chain

    def test_missing_ca_atoms_skipped(self):
        """Residue 2 lacks CA: pairs 1-2 and 2-3 are skipped; 3-4 measured."""
        atoms = (_ca_residue(1, 0.0, with_ca=False) +
                 _ca_residue(2, 3.8, with_ca=False) +
                 _ca_residue(3, 7.6, with_ca=True) +
                 _ca_residue(4, 11.4, with_ca=True))
        res = _audit(atoms)
        assert res["status"] == "PASS"
        assert res["sample"] == 1
        assert res["observed"] == pytest.approx(3.8, abs=0.01)

    def test_all_ca_missing_not_measured(self):
        atoms = (_ca_residue(1, 0.0, with_ca=False) +
                 _ca_residue(2, 3.8, with_ca=False))
        res = _audit(atoms)
        assert res["status"] == "NOT_MEASURED"
        assert res["observed"] is None
        assert res["sample"] == 0

    def test_genuine_coordinate_discontinuity_still_vetoes(self):
        """Adjacent-numbered residues with a 10 A CA-CA separation are a
        real chain break and must veto."""
        atoms = _ca_residue(1, 0.0) + _ca_residue(2, 3.8) + _ca_residue(3, 13.8)
        res = _audit(atoms)
        assert res["status"] == "VETO"
        assert res["observed"] == pytest.approx(10.0, abs=0.01)
        assert res["sample"] == 2

    def test_v23_breach_structures_were_numbering_gaps(self):
        """Regression anchor for the v24 finding: the v23 report listed
        2ACE/1PPE/1BRS/5A63/5T4O/6NB6/7JTL as LAW-160 breaches. Under
        correct adjacency, a pure numbering gap with normal spacing
        cannot breach — demonstrated here with the 45->67 shape."""
        atoms = (_ca_residue(44, 0.0) + _ca_residue(45, 3.8) +
                 _ca_residue(67, 170.0) + _ca_residue(68, 173.8))
        res = _audit(atoms)
        assert res["status"] == "PASS"
        assert res["sample"] == 2
