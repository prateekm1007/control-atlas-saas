"""LAW-182 advisory, end-to-end contract tests (P4).

LAW-182 is a heuristic (non-deterministic) law: it must never appear among
deterministic blocking laws, never flip the authoritative verdict to VETO, and
be labeled advisory in the response/PDF surfaces. Missing or failed LAW-182
measurements must not be misrepresented as favorable scientific evidence, and
other deterministic failures must continue to veto normally.
"""
import pytest

from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws
from tos.governance.station_sop import LAW_CANON
from tos.schemas.response_v1 import PhysicsMeasurement


def _passing():
    return {
        lid: {
            "observed": 0.0,
            "deviation": "0.0",
            "sample": 100,
            "status": "PASS",
        }
        for lid in LAW_CANON
    }


def _adjudicate(measurements, coverage=100.0):
    return adjudicate_laws(AdjudicationInput(
        raw_measurements=measurements,
        coverage=coverage,
        is_experimental=False,
        method_type="predicted",
        architecture="NONE",
    ))


def test_law_182_adverse_heuristic_never_vetoes():
    m = _passing()
    m["LAW-182"] = {"observed": 0.1, "deviation": "999", "sample": 10, "status": "VETO"}
    res = _adjudicate(m)
    assert res.verdict == "PASS"
    row = next(r for r in res.law_rows if r["law_id"] == "LAW-182")
    assert row["method"] == "heuristic"
    assert row["status"] == "FAIL (Advisory)"
    # It must not appear among deterministic blocking laws.
    assert "LAW-182" not in " ".join(res.failing_deterministic)


def test_law_182_missing_is_not_favorable_evidence():
    m = _passing()
    del m["LAW-182"]
    res = _adjudicate(m)
    assert res.verdict == "PASS"
    row = next(r for r in res.law_rows if r["law_id"] == "LAW-182")
    assert row["method"] == "heuristic"
    assert row["status"] == "NOT_MEASURED"
    # Only the other heuristics (LAW-155, LAW-200) still pass; LAW-182 missing
    # adds nothing to the advisory pass count.
    assert res.heur_passed == 2
    assert res.advisory_score == 66


def test_law_182_failed_is_not_favorable_evidence():
    m = _passing()
    m["LAW-182"] = {
        "observed": None, "deviation": None, "sample": 0,
        "status": "MEASUREMENT_FAILED", "reason": "no hydrophobic CA coordinates",
    }
    res = _adjudicate(m)
    assert res.verdict == "PASS"
    row = next(r for r in res.law_rows if r["law_id"] == "LAW-182")
    assert row["status"] == "MEASUREMENT_FAILED"
    assert res.heur_passed == 2


def test_deterministic_failure_still_vetoes_with_advisory_law182():
    m = _passing()
    m["LAW-182"] = {"observed": 0.1, "deviation": "999", "sample": 10, "status": "VETO"}
    m["LAW-130"] = {"observed": 99.0, "deviation": "999", "sample": 10, "status": "VETO"}
    res = _adjudicate(m)
    assert res.verdict == "VETO"
    assert any("LAW-130" in item for item in res.failing_deterministic)


def test_law_182_advisory_row_is_non_deterministic_for_pdf_table():
    """The PDF advisory table selects rows where method != 'deterministic'; an
    adverse LAW-182 row must land there and be labeled Heuristic/Advisory."""
    m = _passing()
    m["LAW-182"] = {"observed": 0.1, "deviation": "999", "sample": 10, "status": "VETO"}
    row = next(r for r in _adjudicate(m).law_rows if r["law_id"] == "LAW-182")
    assert row["method"] != "deterministic"


def test_response_schema_accepts_fail_advisory():
    p = PhysicsMeasurement(
        law_id="LAW-182",
        title="Hydrophobic Burial",
        status="FAIL (Advisory)",
        method="heuristic",
        observed=0.1,
        threshold=0.3,
        operator=">=",
        units="ratio",
        deviation="999",
        sample_size=10,
        scope="chain",
        principle="P",
    )
    assert p.status == "FAIL (Advisory)"