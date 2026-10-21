"""
Round-25: the NOT_MEASURED contract, end to end.

Covers the PDF reporting semantics (the rendered content, not the helper),
the response-contract cross-field validation, the top-level response
serialisation, and the API/batch preservation of the nullable observation.
"""
import os
import sys

import pytest
from pydantic import ValidationError

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.schemas.response_v1 import (
    AdjudicationVerdict,
    PhysicsMeasurement,
    ToscaniniResponse,
)
from tos.governance.station_sop import LAW_CANON

import tos.forensic_artifacts.pdf_generator as pdfgen


# ─────────────────────────────────────────────────────────── PDF rendering

def _law(law_id, status, observed, method="deterministic", measured=True,
         deviation="0.0", sample=10):
    canon = LAW_CANON[law_id]
    return {
        "law_id": law_id, "title": canon["title"], "status": status,
        "method": method, "measured": measured, "observed": observed,
        "threshold": canon["threshold"], "operator": canon["operator"],
        "units": canon["unit"], "deviation": deviation, "sample_size": sample,
        "scope": canon["scope"], "principle": canon.get("principle", "N/A"),
    }


def _payload(laws, verdict="PASS", coverage=95.0):
    return {
        "verdict": {"binary": verdict, "deterministic_score": 90,
                     "advisory_score": 0, "physical_score": 90,
                     "confidence_score": 90.0, "det_passed": 8, "det_total": 12,
                     "heur_passed": 0, "heur_total": 0, "coverage_pct": coverage,
                     "suppression_reason": None, "primary_cause": None,
                     "recommended_action": None, "veto_law": None},
        "tier1": {"laws": laws},
        "tier3": {"probability": 80},
        "characterization": {"source_type": "predicted", "total_atoms": 100,
                             "chains": 1, "residues": 10, "coverage_pct": coverage},
        "provenance": {"source": "test", "hash": "0" * 64, "byte_count": 1},
        "confidence_meta": {},
        "governance": {"audit_id": "TEST-NOT-MEASURED", "governance_fingerprint": {}},
        "strategic_math": {"s6": 0.9, "w_arch": 1.0, "m_s8": 0.0,
                           "p_score": 90, "suppression": None},
    }


def _render(payload):
    """Render the dossier and return the text that was actually drawn."""
    drawn = []
    orig_cell = pdfgen.ToscaniniDossier.cell
    orig_multi = pdfgen.ToscaniniDossier.multi_cell

    def cell(self, *a, **kw):
        txt = kw.get("txt", a[2] if len(a) > 2 else None)
        if txt is not None:
            drawn.append(str(txt))
        return orig_cell(self, *a, **kw)

    def multi_cell(self, *a, **kw):
        txt = kw.get("txt", a[2] if len(a) > 2 else None)
        if txt is not None:
            drawn.append(str(txt))
        return orig_multi(self, *a, **kw)

    pdfgen.ToscaniniDossier.cell = cell
    pdfgen.ToscaniniDossier.multi_cell = multi_cell
    try:
        pdfgen.generate_v21_dossier(payload)
    finally:
        pdfgen.ToscaniniDossier.cell = orig_cell
        pdfgen.ToscaniniDossier.multi_cell = orig_multi
    return drawn


def test_pdf_renders_the_unmeasured_section_without_violations():
    """An otherwise compliant structure still SHOWS the NOT MEASURED section."""
    laws = [_law("LAW-160", "NOT_MEASURED", None, measured=False,
                 deviation=None, sample=0)]
    for lid in ("LAW-100", "LAW-120", "LAW-125", "LAW-130", "LAW-135",
                "LAW-145", "LAW-150", "LAW-155", "LAW-170", "LAW-182",
                "LAW-195", "LAW-200"):
        laws.append(_law(lid, "PASS", 0.1))
    drawn = _render(_payload(laws))
    text = "\n".join(drawn)
    assert "NOT MEASURED" in text, "the unmeasured section was not rendered"
    assert "LAW-160" in text
    assert "not measured" in text
    # no invented observation
    assert "None" not in text


def test_pdf_never_names_an_unmeasured_law_as_a_veto_cause():
    """A genuine VETO names only the real violation."""
    laws = [_law("LAW-160", "NOT_MEASURED", None, measured=False,
                 deviation=None, sample=0),
            _law("LAW-120", "VETO", 24.85)]
    drawn = _render(_payload(laws, verdict="VETO"))
    text = "\n".join(drawn)
    veto_lines = [t for t in drawn if "DETERMINISTIC VETO:" in t]
    assert veto_lines, "no veto banner was rendered"
    assert "LAW-120" in veto_lines[0]
    assert "LAW-160" not in veto_lines[0]
    # the unmeasured law is still reported, but not as a violation
    assert "NOT MEASURED" in text


def test_pdf_renders_an_advisory_unmeasured_law_as_not_measured():
    """LAW-160 is advisory_experimental for the experimental modalities."""
    laws = [_law("LAW-160", "NOT_MEASURED", None, method="advisory_experimental",
                 measured=False, deviation=None, sample=0),
            _law("LAW-100", "PASS", 0.1, method="advisory_experimental")]
    drawn = _render(_payload(laws))
    text = "\n".join(drawn)
    assert "None" not in text, "an unmeasured law was rendered as None"
    assert "not measured" in text


def test_pdf_qualifies_the_compliance_ledger_when_something_is_unmeasured():
    """The ledger must not imply that every deterministic law was measured."""
    laws = [_law("LAW-160", "NOT_MEASURED", None, measured=False,
                 deviation=None, sample=0)]
    for lid in ("LAW-100", "LAW-120", "LAW-125", "LAW-130", "LAW-135",
                "LAW-145", "LAW-150", "LAW-155", "LAW-170", "LAW-182",
                "LAW-195", "LAW-200"):
        laws.append(_law(lid, "PASS", 0.1))
    drawn = _render(_payload(laws))
    text = "\n".join(drawn)
    assert "All deterministic coordinate measurements within thresholds." not in text, \
        "the unqualified ledger wording is still used with an unmeasured law"
    assert "NOT MEASURED" in text


# ─────────────────────────────────────────── response-contract validation

def _row(**kw):
    base = dict(law_id="LAW-160", title="Chain Integrity", status="PASS",
                 method="deterministic", measured=True, observed=3.8, threshold=4.5,
                 operator="<=", units="A", deviation="0.0", sample_size=2,
                 scope="sequential C-alphas", principle="N/A")
    base.update(kw)
    return base


def test_the_model_rejects_contradictory_measurement_rows():
    """measured=false with a number, measured=true with null, and an
    inconsistent NOT_MEASURED row are all rejected."""
    with pytest.raises(ValidationError):
        PhysicsMeasurement(**_row(status="NOT_MEASURED", measured=False,
                                  observed=3.8, deviation=None, sample_size=0))
    with pytest.raises(ValidationError):
        PhysicsMeasurement(**_row(measured=True, observed=None))
    with pytest.raises(ValidationError):
        PhysicsMeasurement(**_row(status="NOT_MEASURED", measured=True,
                                  observed=None, deviation=None, sample_size=0))
    with pytest.raises(ValidationError):
        PhysicsMeasurement(**_row(status="NOT_MEASURED", measured=False,
                                  observed=None, deviation="0.0", sample_size=0))
    with pytest.raises(ValidationError):
        PhysicsMeasurement(**_row(status="NOT_MEASURED", measured=False,
                                  observed=None, deviation=None, sample_size=3))
    # the consistent unmeasured row is accepted
    ok = PhysicsMeasurement(**_row(status="NOT_MEASURED", measured=False,
                                   observed=None, deviation=None, sample_size=0))
    assert ok.observed is None and ok.measured is False


def _minimal(model, overrides=None):
    """A minimal VALID instance of a pydantic model (recursively)."""
    data = {}
    for name, f in model.model_fields.items():
        if not f.is_required():
            continue
        ann = f.annotation
        if isinstance(ann, type) and hasattr(ann, "model_fields"):
            data[name] = _minimal(ann)
        elif ann is float:
            data[name] = 0.0
        elif ann is int:
            data[name] = 0
        elif ann is bool:
            data[name] = False
        else:
            data[name] = "x"
    if overrides:
        data.update(overrides)
    return data


def test_the_top_level_response_accepts_and_serializes_an_unmeasured_law():
    """A complete response with an unmeasured law validates and round-trips."""
    laws = [_row(), _row(status="NOT_MEASURED", measured=False, observed=None,
                         deviation=None, sample_size=0)]
    payload = _minimal(
        ToscaniniResponse,
        {
            "verdict": _minimal(AdjudicationVerdict, {
                "binary": "INDETERMINATE", "deterministic_score": 80,
                "advisory_score": 0, "physical_score": 80,
                "confidence_score": 80.0, "det_passed": 8, "det_total": 12,
                "heur_passed": 0, "heur_total": 0, "coverage_pct": 95.0}),
            "tier1": {"laws": laws},
        },
    )
    resp = ToscaniniResponse(**payload)
    dumped = resp.model_dump()
    unmeasured = [l for l in dumped["tier1"]["laws"] if l["status"] == "NOT_MEASURED"]
    assert len(unmeasured) == 1
    assert unmeasured[0]["observed"] is None
    assert unmeasured[0]["measured"] is False
    assert unmeasured[0]["sample_size"] == 0
    # the observation is NOT converted to zero
    measured = [l for l in dumped["tier1"]["laws"] if l["status"] == "PASS"]
    assert measured[0]["observed"] == 3.8
