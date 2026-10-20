"""
WP-5 — unique, unguessable IDs.

- Audit ID: full SHA-256 of the canonical coordinate hash (64 hex chars),
  deterministic (identical structures share an ID — owner decision).
- Job ID: full uuid4().hex (128 bits).
- Legacy 8-character IDs must still resolve against the audit store.
"""
import hashlib
import os
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")

import main as main_mod
import tos.security.auth as _auth
from tos.storage.audit_store import store_audit_result, get_audit_result

if not _auth.ALLOWED_KEYS:
    _auth.ALLOWED_KEYS = {"test-key-xyz"}

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"
client = TestClient(main_mod.app, raise_server_exceptions=False)
KEY = {"X-API-Key": "test-key-xyz"}


@pytest.fixture(autouse=True)
def _fresh_data(tmp_path, monkeypatch):
    monkeypatch.setenv("TOSCANINI_DATA_DIR", str(tmp_path / "wp5"))


class TestAuditIds:

    def _ingest_1crn(self):
        return client.post(
            "/ingest",
            data={"mode": "experimental", "candidate_id": "1CRN_TEST"},
            files={"file": ("real_pdb_1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
            headers=KEY,
        )

    def test_audit_id_is_full_sha256(self):
        r = self._ingest_1crn()
        assert r.status_code == 200, r.text
        gov = r.json()["governance"]
        audit_id = gov["audit_id"]
        assert len(audit_id) == 64
        int(audit_id, 16)  # valid hex
        assert audit_id == audit_id.lower()

    def test_audit_id_derives_deterministically_from_coord_hash(self):
        """OD-07 (Round 14 correction): deterministic derivation from the
        coordinate hash COMBINED with the canon and matrix hashes — the
        bare sha256(coord_hash) form is superseded (it kept the same ID
        across law-semantics changes)."""
        r = self._ingest_1crn()
        body = r.json()
        coord_hash = body["provenance"]["hash"]
        fp = body["governance"]["governance_fingerprint"]
        from tos.governance.audit_identity import compute_audit_id
        assert body["governance"]["audit_id"] == compute_audit_id(
            coord_hash, fp["canon_hash"], fp["matrix_hash"])
        assert body["governance"]["audit_id"] != \
            hashlib.sha256(coord_hash.encode()).hexdigest()

    def test_identical_structures_share_one_audit_id(self):
        """Owner decision: deterministic dedup kept (brief §4.4)."""
        a = self._ingest_1crn().json()["governance"]["audit_id"]
        b = self._ingest_1crn().json()["governance"]["audit_id"]
        assert a == b
        # The audit store holds exactly one record under that ID.
        assert get_audit_result(a) is not None

    def test_new_id_stored_and_resolvable(self):
        r = self._ingest_1crn()
        audit_id = r.json()["governance"]["audit_id"]
        assert get_audit_result(audit_id) is not None


class TestLegacyIds:

    def test_legacy_8char_id_still_resolves(self):
        """Existing 8-char records resolve unchanged (migration path)."""
        legacy = "AB12CD34"
        store_audit_result(legacy, {"verdict": {"binary": "VETO"}, "legacy": True})
        rec = get_audit_result(legacy)
        assert rec is not None and rec["legacy"] is True

    def test_legacy_and_new_ids_coexist(self):
        store_audit_result("AB12CD34", {"verdict": {"binary": "VETO"}, "legacy": True})
        new_id = "e" * 64
        store_audit_result(new_id, {"verdict": {"binary": "PASS"}, "legacy": False})
        assert get_audit_result("AB12CD34")["legacy"] is True
        assert get_audit_result(new_id)["legacy"] is False


class TestJobIds:

    def test_job_id_is_full_uuid_hex(self):
        """B2 submit (execution engine unavailable in the brain container)
        falls back to the local generator — must be full uuid4().hex."""
        r = client.post(
            "/refinement/submit",
            data={"audit_id": "SOMEAUDIT", "protocol": "openmm", "user_email": "u@x.com"},
            files={"file": ("real_pdb_1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
            headers=KEY,
        )
        assert r.status_code == 200, r.text
        job_id = r.json()["job_id"]
        assert len(job_id) == 32
        int(job_id, 16)  # valid hex
        assert job_id == job_id.lower()
