"""P7 regression tests: canonical OD-07 audit-ID binding across the API path.

The audit_id in the API response must equal compute_audit_id(coord_hash,
canon_hash, matrix_hash) — a full SHA-256 of a canonical serialization that
binds all three identity components. Identical inputs yield identical IDs;
changing any one component changes the ID. This pins the runtime API path,
not just the pure function.
"""
import hashlib
import os
from pathlib import Path
import pytest

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")
import main as main_mod
import tos.security.auth as _auth
from fastapi.testclient import TestClient

if not _auth.ALLOWED_KEYS:
    _auth.ALLOWED_KEYS = {"test-key-xyz"}

KEY = {"X-API-Key": "test-key-xyz"}
FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"
client = TestClient(main_mod.app, raise_server_exceptions=False)


def _ingest(candidate_id):
    return client.post(
        "/ingest",
        data={"mode": "experimental", "candidate_id": candidate_id},
        files={"file": ("1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
        headers=KEY,
    )


class TestOD07ApiBinding:

    def test_audit_id_binds_coord_canon_matrix_in_response(self):
        r = _ingest("P7_A")
        assert r.status_code == 200, r.text
        body = r.json()
        gov = body["governance"]
        fp = gov["governance_fingerprint"]
        from tos.engine.identity import compute_audit_id
        expected = compute_audit_id(
            body["provenance"]["hash"], fp["canon_hash"], fp["matrix_hash"])
        assert gov["audit_id"] == expected, "audit_id must bind coord+canon+matrix"
        assert len(gov["audit_id"]) == 64
        assert all(c in "0123456789abcdef" for c in gov["audit_id"])
        # must NOT be the bare sha256(coord_hash) form
        bare = hashlib.sha256(body["provenance"]["hash"].encode()).hexdigest()
        assert gov["audit_id"] != bare

    def test_identical_structures_share_audit_id(self):
        b1 = _ingest("P7_B1").json()
        b2 = _ingest("P7_B2").json()
        # identical coordinate bytes -> identical coordinate hash -> same audit id
        assert b1["provenance"]["hash"] == b2["provenance"]["hash"]
        assert b1["governance"]["audit_id"] == b2["governance"]["audit_id"]
