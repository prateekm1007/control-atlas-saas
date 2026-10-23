"""
Tests for usage telemetry logging.
Tests the logger functions directly — no HTTP, no file I/O.
"""
import hashlib
import pytest
from unittest.mock import MagicMock, patch, call
from tos.telemetry.usage_logger import (
    build_usage_record,
    log_usage,
    _hash_api_key,
)


def _sample_payload():
    """Minimal valid payload matching ToscaniniResponse structure."""
    return {
        "verdict": {
            "binary": "PASS",
            "deterministic_score": 83,
            "coverage_pct": 100.0,
        },
        "provenance": {
            "hash": "7cba39b01b5e27def61e3499bffda8d595095bd2afa7c20d1c21c69122e264ec",
            "byte_count": 49491,
        },
        "governance": {
            "audit_id": "TEST1234",
            "governance_fingerprint": {
                "canon_hash": "6a9cd4b4349b81de",
                "matrix_hash": "ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4",
            },
        },
    }


class TestBuildUsageRecord:

    def test_contains_all_required_fields(self):
        record = build_usage_record(_sample_payload(), "/ingest")
        required = {
            "timestamp", "structure_hash", "verdict", "det_score",
            "coverage_pct", "api_key_hash", "canon_hash",
            "matrix_hash", "endpoint", "byte_count",
        }
        assert required == set(record.keys()), (
            f"Missing: {required - set(record.keys())}"
        )

    def test_verdict_captured(self):
        record = build_usage_record(_sample_payload(), "/ingest")
        assert record["verdict"] == "PASS"

    def test_det_score_captured(self):
        record = build_usage_record(_sample_payload(), "/ingest")
        assert record["det_score"] == 83

    def test_canon_hash_matches_payload(self):
        record = build_usage_record(_sample_payload(), "/ingest")
        assert record["canon_hash"] == "6a9cd4b4349b81de"

    def test_matrix_hash_truncated(self):
        record = build_usage_record(_sample_payload(), "/ingest")
        assert len(record["matrix_hash"]) == 16

    def test_structure_hash_truncated(self):
        record = build_usage_record(_sample_payload(), "/ingest")
        assert len(record["structure_hash"]) == 16

    def test_endpoint_recorded(self):
        record = build_usage_record(_sample_payload(), "/v1/batch")
        assert record["endpoint"] == "/v1/batch"

    def test_no_api_key_returns_none(self):
        record = build_usage_record(_sample_payload(), "/ingest", api_key=None)
        assert record["api_key_hash"] is None

    def test_api_key_hashed_not_raw(self):
        record = build_usage_record(
            _sample_payload(), "/ingest", api_key="secret-key-123"
        )
        assert record["api_key_hash"] is not None
        assert record["api_key_hash"] != "secret-key-123"
        assert len(record["api_key_hash"]) == 16

    def test_byte_count_captured(self):
        record = build_usage_record(_sample_payload(), "/ingest")
        assert record["byte_count"] == 49491


class TestHashApiKey:

    def test_none_returns_none(self):
        assert _hash_api_key(None) is None

    def test_empty_returns_none(self):
        assert _hash_api_key("") is None

    def test_returns_truncated_sha256(self):
        result = _hash_api_key("test-key")
        assert len(result) == 16
        # Verify it's actually a SHA-256 prefix
        expected = hashlib.sha256("test-key".encode()).hexdigest()[:16]
        assert result == expected

    def test_different_keys_different_hashes(self):
        h1 = _hash_api_key("key-alpha")
        h2 = _hash_api_key("key-beta")
        assert h1 != h2


class TestLogUsage:

    def test_calls_nkg_log_usage(self):
        mock_nkg = MagicMock()
        log_usage(mock_nkg, _sample_payload(), "/ingest")
        mock_nkg.log_usage.assert_called_once()

        record = mock_nkg.log_usage.call_args[0][0]
        assert record["verdict"] == "PASS"
        assert record["endpoint"] == "/ingest"

    def test_failure_does_not_raise(self):
        mock_nkg = MagicMock()
        mock_nkg.log_usage.side_effect = IOError("disk full")
        # Should not raise
        log_usage(mock_nkg, _sample_payload(), "/ingest")

    def test_batch_logging_per_structure(self):
        """Simulate what main.py does for batch: log each structure."""
        mock_nkg = MagicMock()
        payloads = [_sample_payload() for _ in range(3)]

        for p in payloads:
            log_usage(mock_nkg, p, "/v1/batch")

        assert mock_nkg.log_usage.call_count == 3
        for c in mock_nkg.log_usage.call_args_list:
            record = c[0][0]
            assert record["endpoint"] == "/v1/batch"

    def test_api_key_passed_through_hashed(self):
        mock_nkg = MagicMock()
        log_usage(mock_nkg, _sample_payload(), "/ingest", api_key="my-secret")

        record = mock_nkg.log_usage.call_args[0][0]
        assert record["api_key_hash"] is not None
        assert record["api_key_hash"] != "my-secret"
