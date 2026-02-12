from typing import Dict, Any
def validate_payload(p: Dict[str, Any]):
    for f in ["verdict", "governance", "provenance", "tier1", "tier3", "architecture", "semantic_version"]:
        if f not in p: raise ValueError(f"Missing {f}")
def freeze_payload_shape(p: Dict[str, Any]):
    p.setdefault("witness_reports", {}); p.setdefault("definitions", {})
    return p
