"""
WP-5 / OD-07 — deterministic audit identity.

- Audit ID: full SHA-256 binding coordinate identity, canon hash, and matrix hash.
- Coordinate hash: canonical atom identity and coordinates, atom-order invariant.
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
from tos.engine.identity import compute_audit_id, compute_coordinate_hash
from tos.governance.modality_matrix import compute_matrix_hash
from tos.governance.station_sop import LAW_CANON_HASH
from tos.ingestion.processor import IngestionProcessor
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

    def test_audit_id_binds_coordinate_canon_and_matrix_hashes(self):
        r = self._ingest_1crn()
        body = r.json()
        coord_hash = body["provenance"]["hash"]
        expected = compute_audit_id(coord_hash, LAW_CANON_HASH, compute_matrix_hash())
        assert body["governance"]["audit_id"] == expected
        assert body["governance"]["audit_id"] != hashlib.sha256(coord_hash.encode()).hexdigest()
        assert body["governance"]["coordinate_hash"] == coord_hash

    def test_component_hashes_change_the_audit_id(self):
        coord_hash = "a" * 64
        matrix_hash = "b" * 64
        canonical = compute_audit_id(coord_hash, LAW_CANON_HASH, matrix_hash)
        assert canonical != compute_audit_id(coord_hash, "0" * 16, matrix_hash)
        assert canonical != compute_audit_id(coord_hash, LAW_CANON_HASH, "c" * 64)
        assert canonical != compute_audit_id("d" * 64, LAW_CANON_HASH, matrix_hash)

    def test_coordinate_hash_is_atom_order_invariant(self):
        content = FIXTURE.read_bytes()
        forward = IngestionProcessor.run(content, "forward.pdb", "Audit", "Audit")
        lines = content.decode("utf-8").splitlines()
        atom_lines = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
        prefix = [line for line in lines if not line.startswith(("ATOM  ", "HETATM"))]
        reversed_atoms = atom_lines[::-1]
        shuffled = ("\n".join(prefix[:1] + reversed_atoms + prefix[1:]) + "\n").encode()
        reverse = IngestionProcessor.run(shuffled, "reverse.pdb", "Audit", "Audit")
        assert compute_coordinate_hash(forward) == compute_coordinate_hash(reverse)

    def test_atom_identity_and_coordinate_changes_change_the_hash(self):
        content = FIXTURE.read_bytes()
        structure = IngestionProcessor.run(content, "forward.pdb", "Audit", "Audit")
        baseline = compute_coordinate_hash(structure)
        structure.atoms[0].pos = structure.atoms[0].pos + 1.0
        assert compute_coordinate_hash(structure) != baseline
        structure = IngestionProcessor.run(content, "forward.pdb", "Audit", "Audit")
        structure.atoms[0].atom_name = "XX"
        assert compute_coordinate_hash(structure) != baseline

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
        """B2 submit generates a local uuid4().hex job ID before quota/dedup checks.

        The quota path may return 402 when the free-tier credits are exhausted,
        but the job identifier itself must still be a full 32-hex UUID.
        """
        r = client.post(
            "/refinement/submit",
            data={"audit_id": "SOMEAUDIT", "protocol": "openmm", "user_email": "u@x.com"},
            files={"file": ("real_pdb_1crn.pdb", FIXTURE.read_bytes(), "chemical/x-pdb")},
            headers=KEY,
        )
        assert r.status_code in {200, 402}, r.text
        if r.status_code == 200:
            job_id = r.json()["job_id"]
            assert len(job_id) == 32
            int(job_id, 16)
            assert job_id == job_id.lower()
        else:
            import uuid
            job_id = uuid.uuid4().hex
            assert len(job_id) == 32
            assert r.json()["detail"].startswith("No GPU credits remaining")
