"""P2: source-backed coverage contract — returned coverage, LAW-105, and the
adjudicator/API path must agree on the same measurement."""
import numpy as np
import pytest

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Atom, ConfidenceProfile, IngestionMetadata, Structure


def _mk(name, chain, seq, icode, x, label_seq, src):
    return Atom(
        atom_name=name, element="C", pos=np.array([x, 0.0, 0.0]),
        res_name="ALA", res_seq=seq, chain_id=chain, insertion_code=icode,
        b_iso=20.0, label_seq_id=label_seq, source_order=src,
    )


def _ala(chain, seq, icode, x, label_seq, src):
    return [
        _mk("N", chain, seq, icode, x, label_seq, src),
        _mk("CA", chain, seq, icode, x + 1.4, label_seq, src),
        _mk("C", chain, seq, icode, x + 2.0, label_seq, src),
        _mk("O", chain, seq, icode, x + 2.5, label_seq, src),
    ]


def _structure(atoms, source_sequences=None):
    meta = IngestionMetadata()
    if source_sequences:
        meta.source_sequences = source_sequences
    return Structure(atoms=atoms, audit_id="x",
                     confidence=ConfidenceProfile([], True, "xray", None),
                     metadata=meta)


def test_returned_coverage_agrees_with_law105_source_backed():
    atoms = _ala("A", 1, "", 0.0, 1, 0) + _ala("A", 3, "", 3.0, 3, 1)
    s = _structure(atoms, {"A": {"positions": 3, "basis": "label_seq_id"}})
    results, coverage, is_exp = Tier1Measurements.run_full_audit(s)
    law105 = results["LAW-105"]
    assert law105["source_backed"] is True
    assert law105["coverage_denominator"] == 3
    assert law105["coverage_numerator"] == 2
    # The returned coverage is the SAME measurement reported by LAW-105.
    assert coverage == pytest.approx(law105["observed"], abs=0.01)
    assert coverage == pytest.approx(66.6667, abs=0.01)
    assert is_exp is True


def test_missing_source_sequence_never_becomes_authoritative_pass():
    atoms = _ala("A", 1, "", 0.0, None, 0) + _ala("A", 2, "", 3.0, None, 1)
    s = _structure(atoms)   # no source sequence
    results, coverage, _ = Tier1Measurements.run_full_audit(s)
    law105 = results["LAW-105"]
    assert law105["source_backed"] is False
    assert "not source-backed" in law105.get("note", "")
    assert coverage == pytest.approx(law105["observed"], abs=0.01)


def test_coverage_cannot_exceed_100():
    atoms = _ala("A", 1, "", 0.0, 1, 0) + _ala("A", 2, "", 3.0, 2, 1)
    s = _structure(atoms, {"A": {"positions": 2, "basis": "label_seq_id"}})
    results, coverage, _ = Tier1Measurements.run_full_audit(s)
    assert coverage <= 100.0
    assert results["LAW-105"]["observed"] <= 100.0
