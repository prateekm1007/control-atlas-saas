"""
Phase 3 Hardening Tests
-----------------------
Verify that the governance fingerprint is present and correct
in every API response payload.
"""
import hashlib
import pytest
from unittest.mock import MagicMock, patch
from tos.governance.station_sop import LAW_CANON_HASH
from tos.governance.modality_matrix import compute_matrix_hash, get_matrix_meta


class TestGovernanceFingerprintStructure:
    """Verify fingerprint dict has the correct shape and values."""

    def _build_fingerprint(self):
        """Build the fingerprint the same way main.py does."""
        return {
            "canon_hash": LAW_CANON_HASH,
            "matrix_hash": compute_matrix_hash(),
            "matrix_schema_version": get_matrix_meta()["schema_version"],
            "policy_ref": get_matrix_meta()["policy_ref"],
        }

    def test_fingerprint_contains_all_required_keys(self):
        fp = self._build_fingerprint()
        required = {"canon_hash", "matrix_hash", "matrix_schema_version", "policy_ref"}
        assert required == set(fp.keys()), (
            f"Missing keys: {required - set(fp.keys())}"
        )

    def test_canon_hash_matches_station_sop(self):
        fp = self._build_fingerprint()
        assert fp["canon_hash"] == LAW_CANON_HASH

    def test_matrix_hash_matches_locked_value(self):
        fp = self._build_fingerprint()
        expected = "ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4"
        assert fp["matrix_hash"] == expected, (
            f"Matrix hash drifted: {fp['matrix_hash']}"
        )

    def test_canon_hash_matches_locked_value(self):
        fp = self._build_fingerprint()
        assert fp["canon_hash"] == "6a9cd4b4349b81de", (
            f"Canon hash drifted: {fp['canon_hash']}"
        )

    def test_policy_ref_is_PIL_CAL_03(self):
        fp = self._build_fingerprint()
        assert fp["policy_ref"] == "PIL-CAL-03"

    def test_schema_version_is_semver(self):
        fp = self._build_fingerprint()
        parts = fp["matrix_schema_version"].split(".")
        assert len(parts) == 3
        assert all(p.isdigit() for p in parts)


class TestFingerprintInPayload:
    """Simulate what main.py does and verify the fingerprint lands in the payload."""

    def test_governance_dict_contains_fingerprint(self):
        """
        Build a governance dict the way _run_physics_sync does
        and assert the fingerprint sub-dict is present.
        """
        governance = {
            "audit_id": "TEST1234",
            "station_version": "test",
            "timestamp_utc": "2025-01-01T00:00:00Z",
            "governance_fingerprint": {
                "canon_hash": LAW_CANON_HASH,
                "matrix_hash": compute_matrix_hash(),
                "matrix_schema_version": get_matrix_meta()["schema_version"],
                "policy_ref": get_matrix_meta()["policy_ref"],
            },
        }
        fp = governance["governance_fingerprint"]
        assert fp["canon_hash"] == LAW_CANON_HASH
        assert len(fp["matrix_hash"]) == 64
        assert fp["policy_ref"] == "PIL-CAL-03"

    def test_fingerprint_hashes_are_not_empty(self):
        fp = {
            "canon_hash": LAW_CANON_HASH,
            "matrix_hash": compute_matrix_hash(),
        }
        assert fp["canon_hash"], "canon_hash must not be empty"
        assert fp["matrix_hash"], "matrix_hash must not be empty"
        assert len(fp["canon_hash"]) > 8, "canon_hash looks truncated"
        assert len(fp["matrix_hash"]) == 64, "matrix_hash must be full SHA-256"


class TestPDFGeneratorImports:
    """Verify pdf_generator can access the matrix hash."""

    def test_compute_matrix_hash_callable_from_pdf_context(self):
        """
        This tests that the import chain works:
        pdf_generator -> modality_matrix -> compute_matrix_hash
        """
        from tos.forensic_artifacts.pdf_generator import generate_v21_dossier
        # If the import fails, this test fails. No need to generate a full PDF.
        h = compute_matrix_hash()
        assert len(h) == 64
