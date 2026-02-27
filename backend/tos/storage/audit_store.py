"""
Toscanini Phase B2 — Audit Result Store
Stores audit results for retrieval during comparison.
Phase 1: JSON file-based. Phase 2: Redis/PostgreSQL.
"""
import json
from pathlib import Path
from typing import Optional, Dict

import os
AUDIT_DIR = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "audits"
AUDIT_DIR.mkdir(parents=True, exist_ok=True)
AUDIT_TTL_DAYS = 30  # Keep audits for 30 days


def store_audit_result(audit_id: str, result: dict) -> None:
    """Store full audit result for later retrieval."""
    audit_file = AUDIT_DIR / f"{audit_id}.json"
    with open(audit_file, "w") as f:
        json.dump(result, f, indent=2, sort_keys=True, default=str)


def get_audit_result(audit_id: str) -> Optional[Dict]:
    """Retrieve stored audit result by ID."""
    audit_file = AUDIT_DIR / f"{audit_id}.json"
    if audit_file.exists():
        with open(audit_file, "r") as f:
            return json.load(f)
    return None
