"""
Unit tests for LAW-160 (Chain Integrity).
"""
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Structure, ConfidenceProfile


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


# ═══════════════════════════════════════════════════════════════════════
# Round-24 (auditor direction 2): sequence-adjacency eligibility.
#
# get_sequential_residues() returns the residues that are PRESENT (sorted,
# de-duplicated). It does NOT guarantee that two successive entries are
# consecutive in the biological sequence, so measuring the C-alpha distance
# between successive entries can mistake a MODELLING GAP (missing intervening
# residues) for a physical chain break - a false veto. These tests pin the
# corrected eligibility, the retained true-break detection, the
# insertion-code semantics (UNCHANGED - the existing _is_sequential() semantics
# are preserved verbatim), and the sample count.
# ═══════════════════════════════════════════════════════════════════════


def _audit(atoms):
    struct = Structure(atoms=atoms, audit_id="TEST",
                       confidence=ConfidenceProfile([], False, "predicted"))
    results, coverage, _ = Tier1Measurements.run_full_audit(struct)
    return results


def test_law_160_missing_residues_are_not_a_chain_break(mock_atom_factory):
    """Residues 1, 2 and 4 with 3 absent: the 2->4 pair spans a modelling gap.

    Before the correction the C-alpha distance across the gap (here 10 A) was
    reported as a chain break (a false veto). The pair is now ineligible, so
    only the genuine 1->2 pair is measured.
    """
    make = mock_atom_factory
    atoms = [
        make("N",  "N", [0.0, 0.0, 0.0],    res_name="ALA", res_seq=1),
        make("CA", "C", [1.47, 0.0, 0.0],    res_name="ALA", res_seq=1),
        make("C",  "C", [2.5, 1.0, 0.0],     res_name="ALA", res_seq=1),
        make("N",  "N", [3.8, 1.0, 0.0],     res_name="GLY", res_seq=2),
        make("CA", "C", [5.27, 1.0, 0.0],    res_name="GLY", res_seq=2),
        make("C",  "C", [6.3, 2.0, 0.0],     res_name="GLY", res_seq=2),
        # residue 3 is ABSENT from the model; residue 4 sits 10 A away
        make("N",  "N", [15.0, 2.0, 0.0],    res_name="ALA", res_seq=4),
        make("CA", "C", [16.47, 2.0, 0.0],   res_name="ALA", res_seq=4),
        make("C",  "C", [17.5, 3.0, 0.0],    res_name="ALA", res_seq=4),
    ]
    res = _audit(atoms)["LAW-160"]
    # only the consecutive 1->2 pair is eligible and measured
    assert res["sample"] == 1
    assert res["status"] == "PASS"
    assert 3.7 <= res["observed"] <= 4.0


def test_law_160_true_chain_break_is_still_veto(broken_backbone):
    """A real, excessive C-alpha distance between CONSECUTIVE residues still vetoes."""
    res = _audit(broken_backbone)["LAW-160"]
    assert res["status"] == "VETO"
    assert res["observed"] > 9.0
    assert res["sample"] == 2


def test_law_160_insertion_code_adjacency_semantics(mock_atom_factory):
    """1, 2A, 2B, 3 - the insertion-code semantics are the PRESERVED ones.

    The existing _is_sequential() semantics are unchanged: a pair is adjacent
    when the residue number advances by one and neither residue carries an
    insertion code, when the residue number is unchanged and the insertion code
    advances, or when the residue number advances by one out of an
    insertion-coded residue. So 2A->2B and 2B->3 are adjacent, while 1->2A is
    not (the base residue 2 is itself insertion-coded, and the preserved
    semantics do not treat a missing base residue as adjacency).
    """
    make = mock_atom_factory
    atoms = [
        make("CA", "C", [1.47, 0.0, 0.0],  res_name="ALA", res_seq=1),
        make("CA", "C", [5.27, 1.0, 0.0],  res_name="GLY", res_seq=2, insertion_code="A"),
        make("CA", "C", [9.07, 2.0, 0.0],  res_name="GLY", res_seq=2, insertion_code="B"),
        make("CA", "C", [12.87, 3.0, 0.0], res_name="ALA", res_seq=3),
    ]
    res = _audit(atoms)["LAW-160"]
    assert res["sample"] == 2          # (2A,2B) and (2B,3)
    assert res["status"] == "PASS"


def test_law_160_insertion_code_transition_to_next_residue(mock_atom_factory):
    """1, 2A, 3: the transition OUT of an insertion-coded residue is adjacent."""
    make = mock_atom_factory
    atoms = [
        make("CA", "C", [1.47, 0.0, 0.0],  res_name="ALA", res_seq=1),
        make("CA", "C", [5.27, 1.0, 0.0],  res_name="GLY", res_seq=2, insertion_code="A"),
        make("CA", "C", [9.07, 2.0, 0.0],  res_name="ALA", res_seq=3),
    ]
    res = _audit(atoms)["LAW-160"]
    assert res["sample"] == 1          # (2A,3); 1->2A is not adjacent
    assert res["status"] == "PASS"


def test_law_160_gap_after_insertion_code_is_skipped(mock_atom_factory):
    """1, 2, 2A, 4: the 2A->4 pair spans a gap and is not eligible."""
    make = mock_atom_factory
    atoms = [
        make("CA", "C", [1.47, 0.0, 0.0],  res_name="ALA", res_seq=1),
        make("CA", "C", [5.27, 1.0, 0.0],  res_name="GLY", res_seq=2),
        make("CA", "C", [9.07, 2.0, 0.0],  res_name="GLY", res_seq=2, insertion_code="A"),
        make("CA", "C", [16.47, 3.0, 0.0], res_name="ALA", res_seq=4),
    ]
    res = _audit(atoms)["LAW-160"]
    assert res["sample"] == 2          # (1,2) and (2,2A); (2A,4) is a gap
    assert res["status"] == "PASS"
    assert 3.7 <= res["observed"] <= 4.0


def test_law_160_sample_counts_only_pairs_with_both_atoms(mock_atom_factory):
    """A residue without a C-alpha atom contributes no pair and no sample."""
    make = mock_atom_factory
    atoms = [
        make("N",  "N", [0.0, 0.0, 0.0],    res_name="ALA", res_seq=1),
        make("CA", "C", [1.47, 0.0, 0.0],    res_name="ALA", res_seq=1),
        make("C",  "C", [2.5, 1.0, 0.0],     res_name="ALA", res_seq=1),
        make("N",  "N", [3.8, 1.0, 0.0],     res_name="GLY", res_seq=2),
        make("C",  "C", [6.3, 2.0, 0.0],      res_name="GLY", res_seq=2),  # no CA
        make("N",  "N", [7.6, 2.0, 0.0],     res_name="ALA", res_seq=3),
        make("CA", "C", [9.07, 2.0, 0.0],    res_name="ALA", res_seq=3),
        make("C",  "C", [10.1, 3.0, 0.0],    res_name="ALA", res_seq=3),
        make("N",  "N", [11.6, 3.0, 0.0],    res_name="ALA", res_seq=4),
        make("CA", "C", [13.07, 3.0, 0.0],   res_name="ALA", res_seq=4),
        make("C",  "C", [14.1, 4.0, 0.0],    res_name="ALA", res_seq=4),
    ]
    res = _audit(atoms)["LAW-160"]
    # (1,2) and (2,3) have no C-alpha on residue 2; only (3,4) is measured
    assert res["sample"] == 1


def test_law_160_unmeasurable_is_not_reported_as_a_pass(mock_atom_factory):
    """No eligible pair => the law is NOT measured; no fabricated PASS.

    The engine emits no LAW-160 measurement, so the adjudicator applies its
    existing missing-measurement default (status FAIL, sample 0) instead of
    the engine claiming a perfect 0.0 A gap. A dedicated "not measured" state
    that does not read as a scientific failure is proposed in the audit
    evidence and awaits the owner's decision - this point is OPEN.
    """
    make = mock_atom_factory
    atoms = [
        make("N",  "N", [0.0, 0.0, 0.0],    res_name="ALA", res_seq=1),
        make("CA", "C", [1.47, 0.0, 0.0],    res_name="ALA", res_seq=1),
        make("C",  "C", [2.5, 1.0, 0.0],     res_name="ALA", res_seq=1),
    ]
    results = _audit(atoms)
    assert "LAW-160" not in results
