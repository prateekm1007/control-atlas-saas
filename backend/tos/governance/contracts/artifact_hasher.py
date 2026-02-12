import hashlib
import json
from typing import Dict, Any
from tos.utils.type_guards import force_bytes, force_str

def sanitize_for_json(obj: Any) -> Any:
    """Institutional Guard: Recursively ensures JSON serializability and byte-safety."""
    if isinstance(obj, dict):
        return {force_str(k): sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize_for_json(i) for i in obj]
    if isinstance(obj, (bytes, bytearray)):
        return obj.hex()
    return obj

ARTIFACT_HASH_ALLOWLIST = {
    "verdict", "provenance", "tier1", "tier3", 
    "architecture", "governance", "definitions", "notary_record"
}

def compute_artifact_hash(payload: Dict[str, Any]) -> str:
    """PILLAR 19: Notary Seal. Guaranteed Byte-Safe and Name-Aligned."""
    filtered = {
        k: sanitize_for_json(payload[k])
        for k in ARTIFACT_HASH_ALLOWLIST if k in payload
    }

    # json.dumps returns a string.
    canonical_json = json.dumps(
        filtered,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True
    )

    # Use the canonical force_bytes before hashing
    return hashlib.sha256(force_bytes(canonical_json)).hexdigest()
