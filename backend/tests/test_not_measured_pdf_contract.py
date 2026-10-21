"""
End-to-end tests for the Round-24 NOT_MEASURED contract.

These are RENDERING-LEVEL assertions: the PDF is generated and its content streams are
decoded, so the tests fail if the report does not actually show the truth about an
unmeasured law. No extra dependency is required (the PDF text is read from the
compressed content streams directly), so the module is collectable in the CI
environment that installs only pytest and the backend requirements.
"""
import io
import os
import re
import sys
import zlib

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pydantic import ValidationError

from tos.forensic_artifacts.pdf_generator import generate_v21_dossier
from tos.schemas.response_v1 import (
    AdjudicationVerdict,
    Characterization,
    ConfidenceMeta,
    Governance,
    GovernanceFingerprint,
    PhysicsMeasurement,
    Provenance,
    StrategicMath,
    Tier1,
    Tier3,
    ToscaniniResponse,
)


# ────────────────────────────── helpers ──────────────────────────────

def _pdf_text(pdf_bytes: bytes) -> str:
    """The visible text of a generated PDF, read from its content streams.

    The PDF is generated with Flate compression, so the page content is
    decompressed and the text-showing operators are flattened. Parentheses (PDF
    string delimiters) are removed and whitespace is collapsed, so a phrase is
    found regardless of how the generator split it into text-showing operators.
    """
    chunks = []
    for match in re.finditer(rb"stream\r?\n(.*?)\r?\nendstream", pdf_bytes, re.S):
        raw = match.group(1)
        try:
            data = zlib.decompress(raw)
        except zlib.error:
            data = raw
        chunks.append(data.decode("latin-1", "replace"))
    text = "\n".join(chunks)
    text = re.sub(r"[()]", " ", text)
    text = re.sub(r"\\[0-9]{3}", " ", text)
    return re.sub(r"\s+", " ", text)


def _law(law_id, title, status, method, observed, threshold=4.5, operator="<=",
         measured=True, deviation=None, sample_size=1, units="A", scope="all"):
    return {
        "law_id": law_id, "title": title, "status": status, "method": method,
        "measured": measured, "observed": observed, "threshold": threshold,
        "operator": operator, "units": units, "deviation": deviation,
        "sample_size": sample_size, "scope": scope,
    }


def _payload(laws, verdict_binary="PASS", confidence=95.0, det_passed=None):
    det = [l for l in laws if l.get("method") == "deterministic"]
    return {
        "verdict": {
            "binary": verdict_binary, "deterministic_score": 100, "advisory_score": 100,
            "physical_score": 100, "confidence_score": confidence,
            "det_passed": len([l for l in det if l["status"] == "PASS"]),
            "det_total": len(det), "heur_passed": 1, "heur_total": 1, "coverage_pct": 90.0,
        },
        "governance": {
            "audit_id": "TEST-AUDIT", "station_version": "2.1", "timestamp_utc": "2026-10-21T00:00:00Z",
            "governance_fingerprint": {
                "canon_hash": "6a9cd4b4349b81de", "matrix_hash": "ed1731ace09a9a3bd24d3659ab6f3fd18",
                "matrix_schema_version": "1.0", "policy_ref": "PIL-CAL-03",
            },
        },
        "provenance": {"source": "test", "hash": "0" * 64, "byte_count": 0},
        "tier1": {"laws": laws},
        "tier3": {"probability": 50},
        "characterization": {"total_atoms": 100, "total_residues": 12, "source_type": "predicted",
                             "resolution": 2.0, "method": "Standard"},
        "confidence_meta": {"source_type": "predicted", "provenance_method": "plddt", "mean": 90.0,
                            "data_available": True},
        "strategic_math": {"s6": 80.0, "w_arch": 20.0, "m_s8": 15.0, "architecture": "MONOMER"},
    }


def _render(laws, **kw):
    pdf = generate_v21_dossier(_payload(laws, **kw))
    assert isinstance(pdf, (bytes, bytearray)), "the generator must return PDF bytes, got %r" % type(pdf)
    assert bytes(pdf)[:5] == b"%PDF-", "the generator did not return a PDF"
    return _pdf_text(bytes(pdf))


def _measurement(**kw):
    base = dict(law_id="LAW-160", title="Chain Integrity", status="PASS", method="deterministic",
                 measured=True, observed=3.9, threshold=4.5, operator="<=", units="A",
                 deviation="0.6", sample_size=2, scope="all")
    base.update(kw)
    return PhysicsMeasurement(**base)


# ─────────────────── A. response-model consistency validation ───────────────────

class TestMeasurementModelConsistency:
    def test_an_ordinary_measured_row_is_accepted(self):
        m = _measurement()
        assert m.measured is True and m.observed == 3.9

    def test_a_truthful_unmeasured_row_is_accepted(self):
        m = _measurement(status="NOT_MEASURED", measured=False, observed=None, deviation=None,
                         sample_size=0)
        assert m.status == "NOT_MEASURED" and m.observed is None and m.sample_size == 0

    def test_measured_false_with_a_number_is_rejected(self):
        with pytest.raises(ValidationError):
            _measurement(status="NOT_MEASURED", measured=False, observed=0.0, deviation=None,
                         sample_size=0)

    def test_measured_true_with_a_null_observation_is_rejected(self):
        with pytest.raises(ValidationError):
            _measurement(observed=None, deviation=None)

    def test_not_measured_with_a_non_zero_sample_is_rejected(self):
        with pytest.raises(ValidationError):
            _measurement(status="NOT_MEASURED", measured=False, observed=None, deviation=None,
                         sample_size=3)

    def test_not_measured_with_a_deviation_is_rejected(self):
        with pytest.raises(ValidationError):
            _measurement(status="NOT_MEASURED", measured=False, observed=None, deviation=1.0,
                         sample_size=0)

    def test_not_measured_that_claims_to_be_measured_is_rejected(self):
        with pytest.raises(ValidationError):
            _measurement(status="NOT_MEASURED", measured=True, observed=None, deviation=None,
                         sample_size=0)

    def test_an_unmeasured_row_with_another_status_is_rejected(self):
        with pytest.raises(ValidationError):
            _measurement(status="FAIL", measured=False, observed=None, deviation=None, sample_size=0)

    def test_the_observation_field_stays_required(self):
        with pytest.raises(ValidationError):
            PhysicsMeasurement(law_id="LAW-160", title="Chain Integrity", status="PASS",
                               method="deterministic", measured=True, threshold=4.5, operator="<=",
                               units="A", deviation=0.0, sample_size=1, scope="all")


# ─────────────────── full response serialization ───────────────────

def test_the_response_serializes_an_unmeasured_law():
    laws = [_measurement(),
            _measurement(law_id="LAW-150", title="Rotamer", status="NOT_MEASURED", measured=False,
                         observed=None, deviation=None, sample_size=0)]
    response = ToscaniniResponse(
        verdict=AdjudicationVerdict(binary="INDETERMINATE", deterministic_score=50,
                                    advisory_score=100, physical_score=100, confidence_score=90.0,
                                    det_passed=1, det_total=2, heur_passed=1, heur_total=1,
                                    coverage_pct=90.0),
        governance=Governance(
            audit_id="TEST", station_version="2.1", timestamp_utc="2026-10-21T00:00:00Z",
            governance_fingerprint=GovernanceFingerprint(
                canon_hash="6a9cd4b4349b81de", matrix_hash="ed1731ace09a9a3bd24d3659ab6f3fd18",
                matrix_schema_version="1.0", policy_ref="PIL-CAL-03")),
        provenance=Provenance(source="test", hash="0" * 64, byte_count=0),
        tier1=Tier1(laws=laws), tier3=Tier3(probability=50),
        characterization=Characterization(total_atoms=10, total_residues=2, source_type="predicted",
                                         resolution=2.0, method="Standard"),
        confidence_meta=ConfidenceMeta(source_type="predicted", provenance_method="plddt", mean=90.0,
                                       data_available=True),
        strategic_math=StrategicMath(s6=80.0, w_arch=20.0, m_s8=15.0, architecture="MONOMER"))
    dumped = response.model_dump()
    row = [l for l in dumped["tier1"]["laws"] if l["law_id"] == "LAW-150"][0]
    assert row["status"] == "NOT_MEASURED"
    assert row["measured"] is False
    assert row["observed"] is None
    assert row["deviation"] is None
    assert row["sample_size"] == 0


# ─────────────────── B. PDF reporting paths ───────────────────

def test_a_compliant_structure_with_an_unmeasured_law_shows_the_dedicated_section():
    laws = [_law("LAW-160", "Chain Integrity", "PASS", "deterministic", 3.9),
            _law("LAW-150", "Rotamer Outliers", "NOT_MEASURED", "deterministic", None,
                 measured=False, deviation=None, sample_size=0)]
    text = _render(laws)
    assert "NOT MEASURED" in text
    assert "LAW-150" in text
    assert "could not be evaluated" in text


def test_the_compliance_wording_does_not_imply_every_law_was_measured():
    laws = [_law("LAW-160", "Chain Integrity", "PASS", "deterministic", 3.9),
            _law("LAW-150", "Rotamer Outliers", "NOT_MEASURED", "deterministic", None,
                 measured=False, deviation=None, sample_size=0)]
    text = _render(laws)
    assert "All MEASURED deterministic coordinate measurements within thresholds" in text
    assert "All deterministic coordinate measurements within thresholds." not in text


def test_an_unmeasured_advisory_law_never_renders_none():
    laws = [_law("LAW-160", "Chain Integrity", "NOT_MEASURED", "advisory_experimental", None,
                 measured=False, deviation=None, sample_size=0),
            _law("LAW-120", "Bond Angle", "PASS", "deterministic", 2.1)]
    text = _render(laws)
    assert "not measured" in text
    assert "None" not in text


def test_a_true_deterministic_veto_stays_the_veto_cause_next_to_an_unmeasured_law():
    laws = [_law("LAW-120", "Bond Angle", "VETO", "deterministic", 30.0),
            _law("LAW-160", "Chain Integrity", "NOT_MEASURED", "deterministic", None,
                 measured=False, deviation=None, sample_size=0)]
    text = _render(laws, verdict_binary="VETO")
    assert "DETERMINISTIC VETO" in text
    assert "LAW-120" in text
    veto_line = [ln for ln in text.split("DETERMINISTIC VETO") if "LAW-120" in ln][0][:120]
    assert "LAW-160" not in veto_line, "an unmeasured law must not be named as a VETO cause"


def test_an_unmeasured_law_is_neither_a_violation_nor_a_refinement_target():
    laws = [_law("LAW-120", "Bond Angle", "PASS", "deterministic", 2.1),
            _law("LAW-160", "Chain Integrity", "NOT_MEASURED", "deterministic", None,
                 measured=False, deviation=None, sample_size=0)]
    text = _render(laws)
    assert "DETERMINISTIC VIOLATION" not in text
    assert "LAW-160" not in text.split("NOT MEASURED")[0]


def test_a_structure_without_unmeasured_laws_keeps_the_plain_compliance_wording():
    laws = [_law("LAW-160", "Chain Integrity", "PASS", "deterministic", 3.9)]
    text = _render(laws)
    assert "All deterministic coordinate measurements within thresholds." in text
    assert "NOT MEASURED" not in text
