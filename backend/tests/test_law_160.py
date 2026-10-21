"""
Unit tests for LAW-160 (Chain Integrity).

The two original tests are kept unchanged. The regressions added for the
Round-26 LAW-160 correction cover:

* a sequence GAP (nonconsecutive residue numbers) must not be measured as a
  sequential pair, and its distance must not trigger a veto;
* only eligible, actually measured Cα pairs contribute to `sample`;
* insertion-code adjacency and the insertion-code transitions follow the
  established adjacency rules already used by the other Tier-1 laws;
* an absent metric is never reported as a valid normal result.
"""
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Structure, ConfidenceProfile


def _struct(atoms):
    return Structure(atoms=atoms, audit_id="TEST",
                     confidence=ConfidenceProfile([], False, "predicted"))


def _law160(atoms):
    results, _coverage, _ = Tier1Measurements.run_full_audit(_struct(atoms))
    return results["LAW-160"]


def _ca(make, res_seq, x, chain_id="A", insertion_code=""):
    return make("CA", "C", [x, 0.0, 0.0], res_name="ALA", res_seq=res_seq,
              chain_id=chain_id, insertion_code=insertion_code)


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


# --------------------------------------------------------------------------
# Round-26 LAW-160 correction regressions
# --------------------------------------------------------------------------

def test_sequence_gap_is_not_measured_as_a_sequential_pair(mock_atom_factory):
    """Residues 1 and 2 are consecutive; residue 67 is far away.

    The 2 -> 67 pair is NOT consecutive in the sequence (residues 3..66 are
    missing), so it is not a chain-integrity measurement. Its distance must
    not be reported and must not veto the structure.
    """
    make = mock_atom_factory
    atoms = [
        _ca(make, 1, 0.0),
        _ca(make, 2, 3.8),
        _ca(make, 67, 40.0),
    ]
    res = _law160(atoms)
    assert res["sample"] == 1, "only the consecutive 1-2 pair may be measured"
    assert res["observed"] <= 4.0, "the sequence gap must not be reported as a distance"
    assert res["status"] == "PASS"


def test_only_eligible_measured_pairs_count_in_sample(mock_atom_factory):
    """Two consecutive pairs plus one sequence gap: sample counts the two."""
    make = mock_atom_factory
    atoms = [
        _ca(make, 1, 0.0),
        _ca(make, 2, 3.8),
        _ca(make, 3, 7.6),
        _ca(make, 40, 60.0),   # nonconsecutive with 3
    ]
    res = _law160(atoms)
    assert res["sample"] == 2, "only the 1-2 and 2-3 pairs are eligible and measured"
    assert res["observed"] <= 4.0
    assert res["status"] == "PASS"


def test_insertion_code_adjacency_is_measured(mock_atom_factory):
    """10 -> 10A -> 11 are consecutive under the established adjacency rules.

    (same residue number, later insertion code) and (next residue number after
    an insertion-coded residue) are both adjacent, so both Cα pairs are
    measured.
    """
    make = mock_atom_factory
    atoms = [
        _ca(make, 10, 0.0),
        _ca(make, 10, 3.8, insertion_code="A"),
        _ca(make, 11, 7.6),
    ]
    res = _law160(atoms)
    assert res["sample"] == 2, "both insertion-code adjacency pairs must be measured"
    assert res["observed"] <= 4.0
    assert res["status"] == "PASS"


def test_insertion_code_transition_after_plain_residue_is_not_adjacent(mock_atom_factory):
    """11 -> 12A is NOT adjacent under the established rules.

    The rule set treats a residue number plus one as adjacent only when the
    FIRST residue carries the insertion code (or neither does). 11 (no
    insertion code) -> 12A is therefore not a consecutive pair and its
    distance must not be measured; the genuine 10-11 pair still is.
    """
    make = mock_atom_factory
    atoms = [
        _ca(make, 10, 0.0),
        _ca(make, 11, 3.8),
        _ca(make, 12, 40.0, insertion_code="A"),
    ]
    res = _law160(atoms)
    assert res["sample"] == 1, "the 11 -> 12A transition is not a sequential pair"
    assert res["observed"] <= 4.0
    assert res["status"] == "PASS"


def test_zero_sample_is_not_reported_as_pass(mock_atom_factory):
    """A structure with no eligible Cα pair has NO measurement.

    Reporting observed=0/sample=0/PASS would present an absent metric as a
    valid normal result. The adjudicator's existing absent-measurement
    representation is status "FAIL" with observed/sample 0.
    """
    make = mock_atom_factory
    atoms = [_ca(make, 5, 0.0)]
    res = _law160(atoms)
    assert res["sample"] == 0
    assert res["observed"] == 0.0
    assert res["status"] == "FAIL", "an absent measurement must not be reported as PASS"


def test_zero_sample_two_nonconsecutive_residues_is_not_pass(mock_atom_factory):
    """Two residues that are not consecutive yield no measurement either."""
    make = mock_atom_factory
    atoms = [
        _ca(make, 1, 0.0),
        _ca(make, 9, 25.0),
    ]
    res = _law160(atoms)
    assert res["sample"] == 0
    assert res["status"] == "FAIL"


def test_excessive_distance_between_consecutive_residues_still_vetoes(mock_atom_factory):
    """The correction must not weaken detection: consecutive residues that are
    genuinely too far apart still VETO."""
    make = mock_atom_factory
    atoms = [
        _ca(make, 1, 0.0),
        _ca(make, 2, 10.0),
    ]
    res = _law160(atoms)
    assert res["sample"] == 1
    assert res["observed"] > 9.0
    assert res["status"] == "VETO"
