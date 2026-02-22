"""
Usage Telemetry Logger
-----------------------
Builds and writes usage records to the NKG usage log.

Design constraints:
  - Never blocks or delays the API response
  - Never stores raw API keys (SHA-256 hash only)
  - Never modifies adjudication output
  - Failure is swallowed with a warning
"""

from __future__ import annotations

import hashlib
import logging
from datetime import datetime, timezone
from typing import Optional

logger = logging.getLogger("toscanini.telemetry")


def _hash_api_key(raw_key: Optional[str]) -> Optional[str]:
    """Hash the API key for logging. Never store raw keys."""
    if not raw_key:
        return None
    return hashlib.sha256(raw_key.encode()).hexdigest()[:16]


def build_usage_record(
    payload: dict,
    endpoint: str,
    api_key: Optional[str] = None,
) -> dict:
    """
    Build a usage telemetry record from an adjudication payload.

    Args:
        payload:  The full ToscaniniResponse dict
        endpoint: "/ingest" or "/v1/batch"
        api_key:  Raw API key (will be hashed before storage)

    Returns:
        A dict ready for JSONL serialization.
    """
    verdict = payload.get("verdict", {})
    provenance = payload.get("provenance", {})
    governance = payload.get("governance", {})
    fingerprint = governance.get("governance_fingerprint", {})

    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "structure_hash": provenance.get("hash", "UNKNOWN")[:16],
        "verdict": verdict.get("binary", "ERROR"),
        "det_score": verdict.get("deterministic_score", 0),
        "coverage_pct": verdict.get("coverage_pct", 0),
        "api_key_hash": _hash_api_key(api_key),
        "canon_hash": fingerprint.get("canon_hash", "UNKNOWN"),
        "matrix_hash": fingerprint.get("matrix_hash", "UNKNOWN")[:16],
        "endpoint": endpoint,
        "byte_count": provenance.get("byte_count", 0),
    }


def log_usage(
    nkg_manager,
    payload: dict,
    endpoint: str,
    api_key: Optional[str] = None,
) -> None:
    """
    Build and write a usage record. Swallows all exceptions.

    Args:
        nkg_manager: The NKG manager instance (has log_usage method)
        payload:     Full adjudication payload
        endpoint:    Endpoint string
        api_key:     Raw API key (optional, hashed before storage)
    """
    try:
        record = build_usage_record(payload, endpoint, api_key)
        nkg_manager.log_usage(record)
    except Exception as e:
        logger.warning(f"Usage logging failed (non-fatal): {e}")
