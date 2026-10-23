"""
Unit tests for LAW-125 (Ramachandran) degenerate torsion handling.
"""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Atom, ConfidenceProfile, Structure


def _make_atom(atom_name, element, pos, res_name="ALA", res_seq=1,
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


def _backbone_triple(res_seqs=(1, 2, 3)):
    atoms = []
    for seq in res_seqs:
        y = float(seq - 1) * 1.0
        atoms.extend([
            _make_atom("N", "N", [3.8 * (seq - 1), y, 0.0], res_seq=seq),
            _make_atom("CA", "C", [3.8 * (seq - 1) + 1.47, y, 0.0], res_seq=seq),
            _make_atom("C", "C", [3.8 * (seq - 1) + 2.5, y + 1.0, 0.0], res_seq=seq),
            _make_atom("O", "O", [3.8 * (seq - 1) + 2.5, y + 1.0, 1.2], res_seq=seq),
        ])
    return atoms


def _audit(atoms):
    struct = Structure(
        atoms=atoms,
        audit_id="LAW125",
        confidence=ConfidenceProfile([], False, "predicted"),
    )
    results, _, _ = Tier1Measurements.run_full_audit(struct)
    return results["LAW-125"]


class TestLaw125DegenerateTorsions:
    def test_single_degenerate_phi_is_not_a_valid_pair(self):
        atoms = _backbone_triple()
        atoms[6].pos = atoms[5].pos.copy()
        res = _audit(atoms)
        assert res["status"] == "NOT_MEASURED"
        assert res["observed"] is None
        assert res["sample"] == 0
        assert "phi/psi" in res["reason"]

    def test_single_degenerate_psi_is_not_a_valid_pair(self):
        atoms = _backbone_triple()
        atoms[9].pos = atoms[7].pos.copy()
        res = _audit(atoms)
        assert res["status"] == "NOT_MEASURED"
        assert res["observed"] is None
        assert res["sample"] == 0

    def test_collinear_torsion_is_not_a_zero_degree_observation(self):
        atoms = _backbone_triple()
        for a in atoms:
            a.pos = np.array([a.pos[0], 0.0, 0.0])
        res = _audit(atoms)
        assert res["status"] == "NOT_MEASURED"
        assert res["observed"] is None
        assert res["sample"] == 0
