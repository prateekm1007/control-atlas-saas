"""
Audit identity — OD-07 (auditor ruling, Round 13; reaffirmed Round 14).

The audit ID is the SHA-256 of ONE canonical serialization that includes
the structure's coordinate hash, the law canon hash, and the modality
matrix hash. Same coordinates under a different canon or matrix produce
a DIFFERENT ID — an audit of the same structure under different law
semantics is a different audit. Identical canonical input under the same
law semantics shares one ID (the WP-5 dedup behavior is kept, not
superseded).

Superseded forms (must not reappear):
- pre-WP-5 8-character md5 refs (retained only as the informational
  ``governance.legacy_audit_id`` for read-only legacy resolution);
- the bare ``sha256(coord_hash)`` form (same ID under different
  canon/matrix — the deviation identified on this branch's first
  implementation, corrected here);
- the withdrawn Round-12 UUID formulation.

This module is the single derivation point so the ingestion path, the
tests, and the migration script can never disagree.
"""
import hashlib
from typing import Optional

_HEX64 = set("0123456789abcdefABCDEF")


def compute_audit_id(structure_digest: str, canon_hash: str, matrix_hash: str) -> str:
    """Full 64-hex SHA-256 over the canonical serialization binding the
    coordinate digest, the canon hash, and the matrix hash.

    The serialization is fixed-order and domain-tagged; every component
    is length-delimited by '|' and no component may contain '|' (hex
    digests never do), so distinct inputs cannot collide by reordering.
    """
    preimage = "|".join([
        "od07-audit-id-v1",
        str(structure_digest),
        str(canon_hash),
        str(matrix_hash),
    ])
    return hashlib.sha256(preimage.encode("utf-8")).hexdigest()


def is_full_audit_id(value: str) -> bool:
    """True when ``value`` is a full 64-char hexadecimal audit ID."""
    if not isinstance(value, str) or len(value) != 64:
        return False
    return all(c in _HEX64 for c in value)


def audit_id_from_record(record: dict) -> Optional[str]:
    """Derive the OD-07 audit ID from a stored audit record's own
    provenance digest and governance fingerprint. Returns None when the
    record predates those fields (never guessed)."""
    try:
        structure_digest = record["provenance"]["hash"]
        fingerprint = record["governance"]["governance_fingerprint"]
        canon_hash = fingerprint["canon_hash"]
        matrix_hash = fingerprint["matrix_hash"]
    except (KeyError, TypeError):
        return None
    if not structure_digest or not canon_hash or not matrix_hash:
        return None
    return compute_audit_id(structure_digest, canon_hash, matrix_hash)
