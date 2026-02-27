"""
Toscanini Phase B1 — Comparison Metadata Storage
Simple JSON-based storage for audit comparisons (Phase 1).
Phase 2: migrate to Redis or PostgreSQL.

Path is evaluated lazily via _get_storage_dir() to support
TOSCANINI_DATA_DIR env var in both Docker and host test environments.
"""
import json
import os
from pathlib import Path
from typing import Optional, Dict, List
from datetime import datetime, timezone


def _get_storage_dir() -> Path:
    """Lazy evaluation — read TOSCANINI_DATA_DIR at call time, not import time."""
    d = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "comparisons"
    d.mkdir(parents=True, exist_ok=True)
    return d


def store_comparison(original_id: str, refined_id: str, metadata: Dict) -> None:
    """Store comparison metadata linking original and refined audits."""
    comparison_file = _get_storage_dir() / f"{original_id}_to_{refined_id}.json"
    data = {
        "original_audit_id": original_id,
        "refined_audit_id":  refined_id,
        "created_at":        datetime.now(timezone.utc).isoformat(),
        **metadata
    }
    with open(comparison_file, 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)


def get_comparison(original_id: str, refined_id: str) -> Optional[Dict]:
    """Retrieve comparison metadata."""
    comparison_file = _get_storage_dir() / f"{original_id}_to_{refined_id}.json"
    if comparison_file.exists():
        with open(comparison_file, 'r') as f:
            return json.load(f)
    return None


def list_comparisons_by_original(original_id: str) -> List[Dict]:
    """Find all refined audits derived from a given original audit."""
    comparisons = []
    for comp_file in _get_storage_dir().glob(f"{original_id}_to_*.json"):
        with open(comp_file, 'r') as f:
            comparisons.append(json.load(f))
    return sorted(comparisons, key=lambda x: x.get("created_at", ""), reverse=True)


def list_user_comparisons(user_email: str) -> List[Dict]:
    """List all comparisons for a user (if email was provided)."""
    comparisons = []
    for comp_file in _get_storage_dir().glob("*.json"):
        with open(comp_file, 'r') as f:
            data = json.load(f)
            if data.get("user_email") == user_email:
                comparisons.append(data)
    return sorted(comparisons, key=lambda x: x.get("created_at", ""), reverse=True)
