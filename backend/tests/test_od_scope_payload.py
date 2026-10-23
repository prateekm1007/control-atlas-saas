"""End-to-end payload/PDF checks for honest missing and retired states."""
import base64
import io
import os
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")

import main as main_mod
from tos.governance.station_sop import LAW_CANON

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"
KEY = {"X-API-Key": "test-key-xyz"}


def _measurements():
    rows = {
        lid: {
            "observed": 0.0,
            "deviation": "0.0",
            "sample": 100,
            "status": "PASS",
        }
        for lid in LAW_CANON
    }
    rows["LAW-130"] = {
        "observed": None,
        "deviation": None,
        "sample": 0,
        "status": "NOT_MEASURED",
        "reason": "synthetic missing-clash measurement",
    }
    return rows


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("TOSCANINI_DATA_DIR", str(tmp_path / "od-payload"))
    return TestClient(main_mod.app)


def test_missing_law_is_honest_in_response_and_pdf(client, monkeypatch):
    monkeypatch.setattr(
        main_mod.Tier1Measurements,
        "run_full_audit",
        staticmethod(lambda structure, user_intent="NONE": (_measurements(), 100.0, False)),
    )

    response = client.post(
        "/ingest",
        data={"mode": "experimental", "candidate_id": "OD_SCOPE_TEST"},
        files={"file": ("real_pdb_1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
        headers=KEY,
    )
    assert response.status_code == 200, response.text
    payload = response.json()

    laws = {law["law_id"]: law for law in payload["tier1"]["laws"]}
    assert laws["LAW-130"]["status"] == "NOT_MEASURED"
    assert laws["LAW-130"]["observed"] is None
    assert payload["verdict"]["binary"] == "INDETERMINATE"
    assert "LAW-130" in payload["verdict"]["reason"]
    assert laws["LAW-155"]["status"] == "RETIRED"
    assert laws["LAW-155"]["method"] == "retired"
    assert payload["confidence_meta"]["schema_version"] == "confidence-meta/1.0.0"
    assert payload["confidence_meta"]["source_type"] == "B-factor"
    assert payload["confidence_meta"]["data_available"] is False
    assert payload["confidence_meta"]["mean"] is None
    assert payload["verdict"]["confidence_score"] is None
    assert payload["confidence_meta"]["b_factor_mean"] is not None
    assert payload["confidence_meta"]["b_factor_unit"] == "B-factor"
    assert payload["confidence_meta"]["b_factor_observation_count"] > 0

    pdf = base64.b64decode(payload["pdf_b64"])
    assert pdf.startswith(b"%PDF")
    assert payload.get("pdf_error") is None

    pypdf = pytest.importorskip("pypdf")
    text = "\n".join(page.extract_text() or "" for page in pypdf.PdfReader(io.BytesIO(pdf)).pages)
    assert "NOT_MEASURED" in text
    assert "RETIRED" in text
    assert "MEASURED LAW-105 RESIDUE COVERAGE" in text
    assert "EPISTEMIC GAP" in text
    assert "DETERMINISTIC VIOLATION(S)" not in text
    assert "All deterministic coordinate measurements within thresholds" not in text


def test_malformed_not_applicable_attempt_agrees_across_payload_and_pdf(client, monkeypatch):
    rows = {
        lid: {"observed": 0.0, "deviation": "0.0", "sample": 100, "status": "PASS"}
        for lid in LAW_CANON
    }
    rows["LAW-195"] = {"observed": None, "deviation": None, "sample": 0,
                       "status": "NOT_APPLICABLE"}
    monkeypatch.setattr(
        main_mod.Tier1Measurements,
        "run_full_audit",
        staticmethod(lambda structure, user_intent="NONE": (rows, 100.0, False)),
    )

    response = client.post(
        "/ingest",
        data={"mode": "experimental", "candidate_id": "OD_SCOPE_MALFORMED_NA"},
        files={"file": ("real_pdb_1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
        headers=KEY,
    )
    assert response.status_code == 200, response.text
    payload = response.json()
    law = next(l for l in payload["tier1"]["laws"] if l["law_id"] == "LAW-195")
    assert law["status"] == "MEASUREMENT_FAILED"
    assert law.get("reason") is not None
    assert "cannot prove NOT_APPLICABLE" in law["reason"]
    assert payload["verdict"]["binary"] == "INDETERMINATE"
    assert "LAW-195" in payload["verdict"]["reason"]
    assert payload["verdict"]["deterministic_score"] < 100

    pdf = base64.b64decode(payload["pdf_b64"])
    pypdf = pytest.importorskip("pypdf")
    text = "\n".join(page.extract_text() or "" for page in pypdf.PdfReader(io.BytesIO(pdf)).pages)
    assert "MEASUREMENT_FAILED" in text
    assert "LAW-195" in text
