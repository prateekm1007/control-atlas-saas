"""
WP-6 — NKG input validation + public field whitelist.

record_audit rejects: empty source bytes, unknown verdicts, and the
SHA-256-of-empty-input hash. Duplicates (same artifact hash + outcome)
are suppressed idempotently. The public /nkg endpoint exposes only the
customer field whitelist — no artifact hashes, no governance IDs.
"""
import json
import os

import pytest
from fastapi.testclient import TestClient

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")

import main as main_mod
import tos.security.auth as _auth
import tos.nkg.manager as nkg_mod
from tos.nkg.manager import NKGManager, EMPTY_SOURCE_SHA256

if not _auth.ALLOWED_KEYS:
    _auth.ALLOWED_KEYS = {"test-key-xyz"}

client = TestClient(main_mod.app, raise_server_exceptions=False)
KEY = {"X-API-Key": "test-key-xyz"}

VALID_HASH = "a" * 64


def _payload(verdict="VETO", source_hash=VALID_HASH, byte_count=1024, **extra):
    return {
        "verdict": {"binary": verdict, "physical_score": 88.0, "confidence_score": 77.0},
        "governance": {"program_id": "PROG-1", "state": "audited"},
        "provenance": {"hash": source_hash, "byte_count": byte_count},
        "tier3": {"probability": 0.42},
        "architecture": {"authoritative_category": "S8"},
        **extra,
    }


@pytest.fixture()
def nkg(tmp_path, monkeypatch):
    """Fresh NKG manager on a temp ledger; reset the module singleton."""
    monkeypatch.setenv("TOSCANINI_NKG_PATH", str(tmp_path / "nkg" / "piu_moat.jsonl"))
    monkeypatch.setattr(nkg_mod, "_nkg", None)
    manager = nkg_mod.get_nkg()
    yield manager
    monkeypatch.setattr(nkg_mod, "_nkg", None)


def _line_count(manager) -> int:
    if not manager.path.exists():
        return 0
    return sum(1 for l in open(manager.path) if l.strip())


class TestRecordAuditValidation:

    def test_missing_verdict_rejected_not_written(self, nkg):
        payload = _payload()
        del payload["verdict"]
        with pytest.raises(ValueError):
            nkg.record_audit(payload)
        assert _line_count(nkg) == 0, "empty payload must not be recorded"

    def test_unknown_verdict_rejected_not_written(self, nkg):
        with pytest.raises(ValueError):
            nkg.record_audit(_payload(verdict="BANANA"))
        assert _line_count(nkg) == 0

    def test_empty_source_hash_rejected_not_written(self, nkg):
        with pytest.raises(ValueError):
            nkg.record_audit(_payload(source_hash=EMPTY_SOURCE_SHA256))
        assert _line_count(nkg) == 0

    def test_missing_source_hash_rejected_not_written(self, nkg):
        with pytest.raises(ValueError):
            nkg.record_audit(_payload(source_hash=""))
        assert _line_count(nkg) == 0

    def test_zero_byte_source_rejected_not_written(self, nkg):
        with pytest.raises(ValueError):
            nkg.record_audit(_payload(byte_count=0))
        assert _line_count(nkg) == 0

    def test_valid_input_written_once(self, nkg):
        assert nkg.record_audit(_payload()) is True
        assert _line_count(nkg) == 1
        record = json.loads(nkg.path.read_text().splitlines()[0])
        assert record["outcome_category"] == "VETO"
        assert record["artifact_hash"] == VALID_HASH[:nkg_mod.AuditConstants.HASH_LEN]

    def test_ineterminate_verdict_recorded(self, nkg):
        assert nkg.record_audit(_payload(verdict="INDETERMINATE")) is True
        assert _line_count(nkg) == 1

    def test_duplicate_suppressed_idempotent(self, nkg):
        assert nkg.record_audit(_payload()) is True
        assert nkg.record_audit(_payload()) is False  # same hash + outcome
        assert _line_count(nkg) == 1, "duplicates must not be written twice"

    def test_same_structure_different_verdict_appends(self, nkg):
        assert nkg.record_audit(_payload(verdict="VETO")) is True
        assert nkg.record_audit(_payload(verdict="PASS")) is True
        assert _line_count(nkg) == 2


class TestPublicNkgEndpoint:

    def test_only_whitelisted_fields_exposed(self, tmp_path, monkeypatch):
        monkeypatch.setenv("TOSCANINI_NKG_PATH", str(tmp_path / "nkg" / "pub.jsonl"))
        monkeypatch.setattr(nkg_mod, "_nkg", None)
        manager = nkg_mod.get_nkg()
        manager.record_audit(_payload(verdict="VETO"))

        r = client.get("/nkg", headers=KEY)
        assert r.status_code == 200
        body = r.json()
        assert body["total_vetoes"] == 1
        rec = body["vetoes"][0]
        assert set(rec) <= set(main_mod.PUBLIC_NKG_FIELDS)
        for forbidden in ("artifact_hash", "program_id", "program_state", "s8_tax_applied"):
            assert forbidden not in rec, f"{forbidden} must not be exposed publicly"
        assert rec["outcome_category"] == "VETO"
        monkeypatch.setattr(nkg_mod, "_nkg", None)
