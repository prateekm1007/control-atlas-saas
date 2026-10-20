"""
Toscanini Phase B2 — Audit Result Store
Stores audit results for retrieval during comparison.
Phase 1: JSON file-based. Phase 2: Redis/PostgreSQL.

Path is evaluated lazily via _get_audit_dir() to support
TOSCANINI_DATA_DIR env var in both Docker and host test environments.

OD-07 (auditor, delegated by owner, 2026-10-20): the audit ID is the
full 64-hex SHA-256 of the canonical input, combined with the canon
hash and matrix hash. Identical input shares an ID. This module also
provides the LEGACY LOOKUP for 8-character IDs: pre-WP-5 records were
stored under their old 8-char md5-based IDs as filenames, and records
created after WP-5 carry their deterministic legacy ref
(``governance.legacy_audit_id``) so old references keep resolving.
Key scope applies (fail-closed): keyed lookups only return records
that carry a matching ``created_by_key_hash`` — the same convention as
``tos.storage.comparisons.list_comparisons_by_key`` (legacy records
without key attribution are not visible to any key).
"""
import hashlib
import json
import os
from pathlib import Path
from typing import Optional, Dict

AUDIT_TTL_DAYS = 30  # Keep audits for 30 days

# OD-07 legacy reference length (the pre-WP-5 md5[:8] scheme).
LEGACY_AUDIT_ID_LEN = 8
# OD-07 canonical identity length (full SHA-256 hex digest).
CANONICAL_AUDIT_ID_LEN = 64


def _get_audit_dir() -> Path:
    """Lazy evaluation — read TOSCANINI_DATA_DIR at call time, not import time."""
    d = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "audits"
    d.mkdir(parents=True, exist_ok=True)
    return d


def store_audit_result(audit_id: str, result: dict) -> None:
    """Store full audit result for later retrieval."""
    audit_file = _get_audit_dir() / f"{audit_id}.json"
    with open(audit_file, "w") as f:
        json.dump(result, f, indent=2, sort_keys=True, default=str)


def get_audit_result(audit_id: str) -> Optional[Dict]:
    """Retrieve stored audit result by ID."""
    audit_file = _get_audit_dir() / f"{audit_id}.json"
    if audit_file.exists():
        with open(audit_file, "r") as f:
            return json.load(f)
    return None


def compute_legacy_audit_id(coord_hash: str) -> str:
    """OD-07: deterministic legacy 8-char reference for a canonical input.

    Exactly the pre-WP-5 ``main.py`` formula —
    ``md5(coord_hash)[:8].upper()`` — so a client holding an old-style
    8-char reference for a given input can always resolve the record
    that input now maps to. Informational only: it is never the
    canonical identity (the 64-hex SHA-256 is) and never the storage
    filename for records created after WP-5.
    """
    return hashlib.md5(coord_hash.encode()).hexdigest()[:8].upper()


def _record_key_hash(record: Dict) -> str:
    """The key attribution stored on a record ('' when unattributed)."""
    return record.get("created_by_key_hash", "") or ""


def _record_legacy_id(record: Dict) -> str:
    """The stored legacy 8-char reference of a record ('' when absent)."""
    gov = record.get("governance", {}) or {}
    return gov.get("legacy_audit_id", "") or ""


def _key_scoped(record: Dict, key_hash: Optional[str]) -> bool:
    """OD-07 key scope, fail-closed.

    ``key_hash=None``  → internal/trusted caller: no key filtering.
    ``key_hash=str``   → the record MUST carry a matching
                         ``created_by_key_hash``; unattributed (legacy
                         or anonymous) records are invisible, mirroring
                         the WP-2 comparisons convention.
    """
    if key_hash is None:
        return True
    if not key_hash:
        return False
    return bool(_record_key_hash(record)) and _record_key_hash(record) == key_hash


def resolve_audit_ref(ref: str,
                      key_hash: Optional[str] = None) -> Optional[Dict]:
    """OD-07 lookup: resolve a canonical 64-hex ID or a legacy 8-char ID.

    Resolution order:
      1. 64-hex ref      → canonical filename ``<ref>.json``.
      2. 8-char legacy ref →
         a. pre-WP-5 records stored under their legacy filename
            ``<REF>.json`` (key-scoped: these carry no key attribution,
            so they are only returned to internal callers), then
         b. post-WP-5 records whose ``governance.legacy_audit_id``
            matches the ref (key-scoped against ``created_by_key_hash``).
    Anything else (wrong length/format) resolves to None.

    Returns the stored record dict, or None when not found / not
    visible to the given key (fail-closed).
    """
    if not ref:
        return None
    ref = ref.strip()

    if len(ref) == CANONICAL_AUDIT_ID_LEN:
        record = get_audit_result(ref)
        if record is not None and _key_scoped(record, key_hash):
            return record
        return None

    if len(ref) == LEGACY_AUDIT_ID_LEN:
        # (a) pre-WP-5 filename scheme (legacy IDs were uppercase).
        record = get_audit_result(ref.upper())
        if record is not None and _key_scoped(record, key_hash):
            return record
        # (b) scan for a post-WP-5 record carrying this legacy ref.
        for audit_file in sorted(_get_audit_dir().glob("*.json")):
            try:
                with open(audit_file, "r") as f:
                    candidate = json.load(f)
            except (OSError, json.JSONDecodeError):
                continue  # corrupted record — never fatal for other refs
            if _record_legacy_id(candidate).upper() != ref.upper():
                continue
            if _key_scoped(candidate, key_hash):
                return candidate
        return None

    # Unknown reference format — not found (fail-closed).
    return None
