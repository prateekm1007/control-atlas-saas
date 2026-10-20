"""
OD-07 — Audit identity and legacy lookup (auditor, delegated by owner,
2026-10-20; ruling REAFFIRMED and corrected in Round 14).

Rules under test:
- The audit ID is the SHA-256 of ONE canonical serialization that
  includes the structure's coordinate hash, the canon hash, and the
  matrix hash. Same coordinates under a different canon or matrix
  produce a DIFFERENT ID (the ruling test, added Round 14). Identical
  input under the same semantics shares one ID.
- The separate audit_identity object of the first implementation is
  REMOVED (it left canon/matrix beside the ID instead of inside it);
  the binding values live in governance_fingerprint.
- A legacy lookup resolves 8-character references read-only: both
  pre-WP-5 filename records and the deterministic
  governance.legacy_audit_id of newer records. No writes on GET.
- Key scope applies, fail-closed: /audit lookups only return records
  created by the calling key; unattributed (legacy/anonymous) records
  are not visible to any key.
- No uuid4()[:8] identifier exists at the ingestion layer.
"""
import hashlib
import io
import os
import zipfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz,other-key-abc")

import main as main_mod
import tos.security.auth as _auth
from tos.storage.audit_store import (
    store_audit_result, get_audit_result,
    resolve_audit_ref, compute_legacy_audit_id,
)
from tos.governance.audit_identity import compute_audit_id
from tos.saas.key_store import hash_key
from tos.governance.station_sop import LAW_CANON_HASH
from tos.governance.modality_matrix import compute_matrix_hash
from tos.ingestion.processor import IngestionProcessor

# OD-07 tests need a second static key for the cross-key (deny) case.
# Adding a key can only widen static acceptance; the primary key used by
# every other test module remains present.
_auth.ALLOWED_KEYS = set(_auth.ALLOWED_KEYS) | {"other-key-abc"}

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"
KEY_A = {"X-API-Key": "test-key-xyz"}
KEY_B = {"X-API-Key": "other-key-abc"}

client = TestClient(main_mod.app, raise_server_exceptions=False)


@pytest.fixture(autouse=True)
def _fresh_data(tmp_path, monkeypatch):
    monkeypatch.setenv("TOSCANINI_DATA_DIR", str(tmp_path / "od07"))


@pytest.fixture(autouse=True)
def _second_key_present():
    """test_api_key_guard.py reloads tos.security.auth mid-session, which
    re-executes the module and resets ALLOWED_KEYS to the env-derived set.
    Re-add the second static key before every test so the cross-key cases
    exercise the LOOKUP scope (not a 401 auth failure)."""
    _auth.ALLOWED_KEYS = set(_auth.ALLOWED_KEYS) | {"other-key-abc"}
    yield


def _ingest_1crn():
    return client.post(
        "/ingest",
        data={"mode": "experimental", "candidate_id": "1CRN_TEST"},
        files={"file": ("real_pdb_1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
        headers=KEY_A,
    )


class TestAuditIdentity:

    def test_audit_id_binds_coord_canon_and_matrix(self):
        """OD-07 ruling: the ID is the SHA-256 of the canonical
        serialization binding the coordinate hash, the canon hash, and
        the matrix hash (the binding values carried by the fingerprint)."""
        r = _ingest_1crn()
        assert r.status_code == 200, r.text
        body = r.json()
        gov = body["governance"]
        fp = gov["governance_fingerprint"]
        assert gov["audit_id"] == compute_audit_id(
            body["provenance"]["hash"], fp["canon_hash"], fp["matrix_hash"])
        assert len(gov["audit_id"]) == 64
        int(gov["audit_id"], 16)
        assert gov["audit_id"] == gov["audit_id"].lower()
        assert fp["canon_hash"] == LAW_CANON_HASH
        assert fp["matrix_hash"] == compute_matrix_hash()
        # The superseded separate identity object is gone.
        assert "audit_identity" not in gov

    def test_identical_input_shares_identity(self):
        a = _ingest_1crn().json()["governance"]["audit_id"]
        b = _ingest_1crn().json()["governance"]["audit_id"]
        assert a == b

    def test_audit_id_is_not_the_bare_coord_hash(self):
        """Round 14 correction: the bare sha256(coord_hash) form is the
        DEVIATION (same ID under different canon/matrix). It is
        superseded and must not reappear."""
        body = _ingest_1crn().json()
        coord_hash = body["provenance"]["hash"]
        assert body["governance"]["audit_id"] != \
            hashlib.sha256(coord_hash.encode()).hexdigest()

    def test_same_coords_different_canon_different_id(self):
        """THE ruling test (Round 14, direction 2): same coordinates,
        different canon hash, DIFFERENT ID."""
        digest = "a" * 64
        matrix = "b" * 64
        id1 = compute_audit_id(digest, "c" * 16, matrix)
        id2 = compute_audit_id(digest, "d" * 16, matrix)
        assert id1 != id2
        assert compute_audit_id(digest, "c" * 16, matrix) == id1  # deterministic

    def test_same_coords_different_matrix_different_id(self):
        """Ruling, matrix component: same coordinates, different matrix
        hash, DIFFERENT ID."""
        digest = "a" * 64
        canon = "c" * 16
        assert compute_audit_id(digest, canon, "b" * 64) != \
            compute_audit_id(digest, canon, "e" * 64)

    def test_same_input_different_canon_changes_response_id(self, monkeypatch):
        """Integration form of the ruling test: the SAME uploaded file
        audited under a different canon hash (simulated law-semantics
        change) yields a different audit_id in the response."""
        first = _ingest_1crn().json()["governance"]["audit_id"]
        real_canon = main_mod.LAW_CANON_HASH
        monkeypatch.setattr(main_mod, "LAW_CANON_HASH", "f" * 16,
                            raising=True)
        try:
            second = _ingest_1crn().json()["governance"]["audit_id"]
        finally:
            monkeypatch.setattr(main_mod, "LAW_CANON_HASH", real_canon,
                                raising=True)
        assert first != second

    def test_legacy_ref_is_deterministic_pre_wp5_scheme(self):
        body = _ingest_1crn().json()
        coord_hash = body["provenance"]["hash"]
        legacy = body["governance"]["legacy_audit_id"]
        assert len(legacy) == 8
        assert legacy == legacy.upper()
        int(legacy, 16)
        # Exactly the pre-WP-5 main.py formula: md5(coord_hash)[:8].upper()
        assert legacy == hashlib.md5(coord_hash.encode()).hexdigest()[:8].upper()
        assert legacy == compute_legacy_audit_id(coord_hash)


class TestNoIngestionIdentifier:

    def test_ingestion_layer_carries_no_identifier(self):
        """OD-07: uuid4()[:8] removed from the ingestion path — the
        parser-level audit_id is empty and NOT random."""
        content = FIXTURE.read_bytes()
        s1 = IngestionProcessor.run(content, "x.pdb", "Audit", "experimental")
        s2 = IngestionProcessor.run(content, "x.pdb", "Audit", "experimental")
        assert s1.audit_id == ""
        assert s2.audit_id == ""
        assert s1.audit_id == s2.audit_id  # nothing random


class TestKeyScopedLookup:

    def test_record_attributed_to_calling_key(self):
        r = _ingest_1crn()
        audit_id = r.json()["governance"]["audit_id"]
        rec = get_audit_result(audit_id)
        assert rec is not None
        assert rec.get("created_by_key_hash") == hash_key("test-key-xyz")

    def test_lookup_by_canonical_id_for_owner(self):
        audit_id = _ingest_1crn().json()["governance"]["audit_id"]
        r = client.get(f"/audit/{audit_id}", headers=KEY_A)
        assert r.status_code == 200 and r.json()["status"] == "success"
        assert r.json()["audit"]["governance"]["audit_id"] == audit_id

    def test_lookup_denied_for_other_key(self):
        """Key scope: another configured key cannot read key A's record."""
        audit_id = _ingest_1crn().json()["governance"]["audit_id"]
        r = client.get(f"/audit/{audit_id}", headers=KEY_B)
        assert r.status_code == 200  # endpoint reachable...
        assert r.json()["status"] == "not_found"  # ...record not visible

    def test_anonymous_records_visible_to_no_key(self):
        """Fail-closed: a record without key attribution (e.g. created by
        an internal flow) is not visible to any keyed lookup."""
        some_id = "f" * 64
        store_audit_result(some_id, {"governance": {"audit_id": some_id}})
        for key in (KEY_A, KEY_B):
            r = client.get(f"/audit/{some_id}", headers=key)
            assert r.json()["status"] == "not_found"
        # Internal (no key filter) callers still resolve it.
        assert resolve_audit_ref(some_id)["governance"]["audit_id"] == some_id


class TestLegacyLookup:

    def test_legacy_8char_resolves_for_owner_key(self):
        body = _ingest_1crn().json()
        legacy = body["governance"]["legacy_audit_id"]
        canonical = body["governance"]["audit_id"]
        r = client.get(f"/audit/{legacy}", headers=KEY_A)
        assert r.status_code == 200 and r.json()["status"] == "success"
        assert r.json()["audit"]["governance"]["audit_id"] == canonical

    def test_legacy_8char_denied_for_other_key(self):
        legacy = _ingest_1crn().json()["governance"]["legacy_audit_id"]
        r = client.get(f"/audit/{legacy}", headers=KEY_B)
        assert r.json()["status"] == "not_found"

    def test_pre_wp5_filename_record_fail_closed_via_api(self):
        """Pre-WP-5 records (8-char filename, no attribution) stay
        resolvable internally but are not exposed to any key."""
        legacy = "AB12CD34"
        store_audit_result(legacy, {"verdict": {"binary": "VETO"},
                                    "legacy": True})
        assert get_audit_result(legacy)["legacy"] is True  # internal path
        for key in (KEY_A, KEY_B):
            r = client.get(f"/audit/{legacy}", headers=key)
            assert r.json()["status"] == "not_found"

    def test_invalid_ref_not_found(self):
        for bad in ("xyz", "", "g" * 20):
            r = client.get(f"/audit/{bad}" if bad else "/audit/%20", headers=KEY_A)
            assert r.json()["status"] == "not_found"

    def test_legacy_lookup_writes_nothing(self):
        """Round 14 reaffirmation: the legacy lookup is READ-ONLY — a GET
        must not create, rewrite, or rename any stored record (stored-
        record rewrites belong to the explicit migration script
        scripts/migrate_audit_ids_od07.py, never to a read path)."""
        body = _ingest_1crn().json()
        legacy = body["governance"]["legacy_audit_id"]
        from tos.storage.audit_store import _get_audit_dir
        before = sorted(p.name for p in _get_audit_dir().iterdir())
        r = client.get(f"/audit/{legacy}", headers=KEY_A)
        assert r.status_code == 200 and r.json()["status"] == "success"
        after = sorted(p.name for p in _get_audit_dir().iterdir())
        assert before == after


class TestRefinementPdbLegacyResolution:

    def test_pdb_served_via_legacy_ref(self):
        body = _ingest_1crn()
        legacy = body.json()["governance"]["legacy_audit_id"]
        r = client.get(f"/refinement/pdb/{legacy}", headers=KEY_A)
        assert r.status_code == 200, r.text
        assert r.content == FIXTURE.read_bytes()

    def test_pdb_served_via_canonical_id(self):
        audit_id = _ingest_1crn().json()["governance"]["audit_id"]
        r = client.get(f"/refinement/pdb/{audit_id}", headers=KEY_A)
        assert r.status_code == 200
        assert r.content == FIXTURE.read_bytes()

    def test_pdb_unknown_ref_404(self):
        r = client.get("/refinement/pdb/ZZZZZZZZ", headers=KEY_A)
        assert r.status_code == 404


class TestBatchKeyAttribution:

    def test_batch_records_carry_key_hash(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as zf:
            zf.writestr("1crn.pdb", FIXTURE.read_bytes())
        r = client.post(
            "/v1/batch",
            data={"mode": "experimental"},
            files={"file": ("batch.zip", buf.getvalue(), "application/zip")},
            headers=KEY_A,
        )
        assert r.status_code == 200, r.text
        results = r.json()["results"]
        assert results and results[0]["success"], r.text
        audit_id = results[0]["response"]["governance"]["audit_id"]
        rec = get_audit_result(audit_id)
        assert rec is not None
        assert rec.get("created_by_key_hash") == hash_key("test-key-xyz")
