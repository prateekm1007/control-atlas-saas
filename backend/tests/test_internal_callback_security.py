"""
WP-3 — /refinement/internal-callback security tests.

Covers the brief's required cases:
  - no signature          -> 401
  - bad signature         -> 401
  - stale timestamp       -> 401
  - malformed hex         -> 422
  - valid signed request  -> 200 and exactly one comparison record

The job store is monkeypatched (no Redis in tests). The original audit is
seeded via the real audit store. The refined structure is the REAL 1CRN
structure (backend/tests/fixtures/real_pdb_1crn.pdb), not a fabricated one.
"""
import hashlib
import hmac as hmac_mod
import json
import os
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

import main as main_mod

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"
SECRET = "test-internal-callback-secret"
ORIGINAL_AUDIT_ID = "WPCB0001"

client = TestClient(main_mod.app, raise_server_exceptions=False)


def _signed_post(body: bytes):
    """POST with a valid HMAC signature over the exact body bytes."""
    ts = str(int(time.time()))
    return client.post("/refinement/internal-callback", content=body,
                       headers={"Content-Type": "application/json",
                                "X-Timestamp": ts,
                                "X-Signature": _sign(body, ts)})


def _sign(body: bytes, timestamp: str, secret: str = SECRET) -> str:
    return hmac_mod.new(
        secret.encode(), timestamp.encode() + b"." + body, hashlib.sha256
    ).hexdigest()


def _body(job_id: str = "JOB00001", hex_override: str = None) -> bytes:
    hexval = hex_override if hex_override is not None else FIXTURE.read_bytes().hex()
    return json.dumps({
        "original_audit_id": ORIGINAL_AUDIT_ID,
        "refined_pdb_hex": hexval,
        "job_id": job_id,
        "protocol": "openmm",
        "user_email": "worker@example.com",
    }).encode()


@pytest.fixture(scope="module", autouse=True)
def _env_and_stores(tmp_path_factory):
    os.environ["INTERNAL_CALLBACK_SECRET"] = SECRET
    data_dir = tmp_path_factory.mktemp("wp3-data")
    os.environ["TOSCANINI_DATA_DIR"] = str(data_dir)
    # Seed the original audit through the real store.
    from tos.storage.audit_store import store_audit_result
    store_audit_result(ORIGINAL_AUDIT_ID, {
        "verdict": {"binary": "VETO", "deterministic_score": 10, "coverage_pct": 55.0},
        "governance": {"audit_id": ORIGINAL_AUDIT_ID},
    })
    yield


@pytest.fixture(autouse=True)
def _job_store(monkeypatch):
    """Fake the Redis job store: JOB00001 exists with a key hash."""
    def fake_get_job_record(job_id):
        if job_id == "JOB00001":
            return {"job_id": job_id, "created_by_key_hash": "abc123", "state": "queued"}
        return None
    monkeypatch.setattr(main_mod, "_get_job_record", fake_get_job_record)


class TestSignatureEnforcement:

    def test_no_signature_rejected_401(self):
        r = client.post("/refinement/internal-callback", content=_body(),
                        headers={"Content-Type": "application/json"})
        assert r.status_code == 401

    def test_bad_signature_rejected_401(self):
        body = _body()
        ts = str(int(time.time()))
        r = client.post("/refinement/internal-callback", content=body,
                        headers={"Content-Type": "application/json",
                                 "X-Timestamp": ts,
                                 "X-Signature": "0" * 64})
        assert r.status_code == 401

    def test_wrong_secret_rejected_401(self):
        body = _body()
        ts = str(int(time.time()))
        r = client.post("/refinement/internal-callback", content=body,
                        headers={"Content-Type": "application/json",
                                 "X-Timestamp": ts,
                                 "X-Signature": _sign(body, ts, secret="attacker-secret")})
        assert r.status_code == 401

    def test_stale_timestamp_rejected_401(self):
        body = _body()
        ts = str(int(time.time()) - 301)  # older than the 300 s window
        r = client.post("/refinement/internal-callback", content=body,
                        headers={"Content-Type": "application/json",
                                 "X-Timestamp": ts,
                                 "X-Signature": _sign(body, ts)})
        assert r.status_code == 401

    def test_tampered_body_rejected_401(self):
        """Signature computed over a different body must not validate."""
        body = _body()
        ts = str(int(time.time()))
        sig = _sign(body, ts)
        tampered = _body(hex_override=(FIXTURE.read_bytes() + b"\nREMARK tampered\n").hex())
        r = client.post("/refinement/internal-callback", content=tampered,
                        headers={"Content-Type": "application/json",
                                 "X-Timestamp": ts, "X-Signature": sig})
        assert r.status_code == 401

    def test_unconfigured_secret_returns_503(self, monkeypatch):
        monkeypatch.setenv("INTERNAL_CALLBACK_SECRET", "")
        r = client.post("/refinement/internal-callback", content=_body(),
                        headers={"Content-Type": "application/json"})
        assert r.status_code == 503


class TestBodyValidation:

    def test_malformed_hex_rejected_422(self):
        body = _body(hex_override="NOT-VALID-HEX")
        r = _signed_post(body)
        assert r.status_code == 422

    def test_hex_without_atom_records_rejected_422(self):
        body = _body(hex_override=b"REMARK nothing here\n".hex())
        r = _signed_post(body)
        assert r.status_code == 422

    def test_missing_field_rejected_422(self):
        body = json.dumps({"original_audit_id": ORIGINAL_AUDIT_ID,
                           "refined_pdb_hex": FIXTURE.read_bytes().hex()}).encode()
        r = _signed_post(body)  # job_id missing
        assert r.status_code == 422

    def test_unknown_job_rejected_422(self):
        r = _signed_post(_body(job_id="JOB_DOES_NOT_EXIST"))
        assert r.status_code == 422

    def test_unknown_original_audit_rejected_422(self, monkeypatch, tmp_path):
        # Point the audit store at an empty dir -> original not found.
        monkeypatch.setenv("TOSCANINI_DATA_DIR", str(tmp_path / "empty"))
        r = _signed_post(_body())
        assert r.status_code == 422


class TestValidCallback:

    def test_valid_request_creates_exactly_one_comparison(self, tmp_path, monkeypatch):
        """End-to-end with the real 1CRN structure: 200 + exactly one
        comparison record, stamped with the job's key hash."""
        from tos.storage.audit_store import store_audit_result
        from tos.storage import comparisons as cmp_mod

        # Fresh, isolated data dir for exact record counting.
        data_dir = tmp_path / "wp3-valid"
        monkeypatch.setenv("TOSCANINI_DATA_DIR", str(data_dir))
        store_audit_result(ORIGINAL_AUDIT_ID, {
            "verdict": {"binary": "VETO", "deterministic_score": 10, "coverage_pct": 55.0},
            "governance": {"audit_id": ORIGINAL_AUDIT_ID},
        })

        r = _signed_post(_body())
        assert r.status_code == 200, r.text
        payload = r.json()
        assert payload["status"] == "success"
        assert payload["refined_audit_id"]

        comp_dir = data_dir / "comparisons"
        records = list(comp_dir.glob("*.json"))
        assert len(records) == 1, f"expected exactly 1 comparison, found {len(records)}"
        record = json.loads(records[0].read_text())
        assert record["created_by_key_hash"] == "abc123"
        assert record["refinement_method"] == "managed_gpu"
