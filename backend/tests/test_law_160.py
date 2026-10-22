"""
Unit tests for LAW-160 (Chain Integrity).
"""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Atom, ConfidenceProfile, Structure


def _ca_atom(res_seq, x, chain_id="A", icode="", res_name="ALA"):
    return Atom(
        atom_name="CA",
        element="C",
        pos=np.array([x, 0.0, 0.0], dtype=float),
        res_name=res_name,
        res_seq=res_seq,
        chain_id=chain_id,
        insertion_code=icode,
        b_iso=90.0,
    )


def _audit(atoms):
    struct = Structure(
        atoms=atoms,
        audit_id="LAW160",
        confidence=ConfidenceProfile([], False, "predicted"),
    )
    results, _, _ = Tier1Measurements.run_full_audit(struct)
    return results["LAW-160"]


def test_law_160_pass(simple_backbone):
    """simple_backbone has ~3.8-3.95Å CA-CA spacing — should PASS."""
    struct = Structure(atoms=simple_backbone, audit_id="TEST",
                       confidence=ConfidenceProfile([], False, "predicted"))
    results, coverage, _ = Tier1Measurements.run_full_audit(struct)

    res = results["LAW-160"]
    assert res["status"] == "PASS"
    # Fixture geometry gives CA-CA between 3.7 and 4.0
    assert 3.7 <= res["observed"] <= 4.0
    # Sample should be CA pairs, not total residues
    assert res["sample"] == 2


def test_law_160_fail(broken_backbone):
    """broken_backbone has a ~10Å gap — should VETO."""
    struct = Structure(atoms=broken_backbone, audit_id="TEST",
                       confidence=ConfidenceProfile([], False, "predicted"))
    results, coverage, _ = Tier1Measurements.run_full_audit(struct)

    res = results["LAW-160"]
    assert res["status"] == "VETO"
    assert res["observed"] > 9.0


def test_law_160_gap_separated_residues_are_not_evaluated():
    res = _audit([_ca_atom(1, 0.0), _ca_atom(3, 8.0)])
    assert res == {
        "observed": None,
        "status": "NOT_MEASURED",
        "sample": 0,
        "reason": "no eligible sequential CA pairs with both coordinates present",
    }


def test_law_160_insertion_codes_are_sequential():
    res = _audit([
        _ca_atom(1, 0.0),
        _ca_atom(1, 3.8, icode="A"),
        _ca_atom(2, 7.6),
    ])
    assert res["status"] == "PASS"
    assert res["sample"] == 2
    assert res["observed"] == pytest.approx(3.8, abs=1e-9)


def test_law_160_chain_boundaries_are_not_measured():
    res = _audit([
        _ca_atom(1, 0.0, chain_id="A"),
        _ca_atom(2, 3.8, chain_id="B"),
    ])
    assert res["status"] == "NOT_MEASURED"
    assert res["sample"] == 0


def test_law_160_missing_coordinate_is_not_a_measurement():
    res = _audit([
        _ca_atom(1, 0.0),
        Atom(
            atom_name="CA", element="C",
            pos=np.array([3.8, 0.0, 0.0], dtype=float),
            res_name="ALA", res_seq=2, chain_id="A",
            insertion_code="", b_iso=90.0,
        ),
    ])
    assert res["status"] == "PASS"
    assert res["sample"] == 1


def test_law_160_genuine_chain_break_vetoes():
    res = _audit([
        _ca_atom(1, 0.0),
        _ca_atom(2, 3.8),
        _ca_atom(3, 20.0),
    ])
    assert res["status"] == "VETO"
    assert res["observed"] == pytest.approx(16.2, abs=1e-9)
    assert res["sample"] == 2
