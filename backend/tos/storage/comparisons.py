"""
Toscanini Phase B1 — Comparison Metadata Storage
Simple JSON-based storage for audit comparisons (Phase 1).
Phase 2: migrate to Redis or PostgreSQL.
"""
import json
from pathlib import Path
from typing import Optional, Dict, List
from datetime import datetime, timezone

import os
STORAGE_DIR = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "comparisons"
STORAGE_DIR.mkdir(parents=True, exist_ok=True)

def store_comparison(original_id: str, refined_id: str, metadata: Dict) -> None:
    """
    Store comparison metadata linking original and refined audits.
    
    Args:
        original_id: Original audit ID (baseline)
        refined_id: Refined structure audit ID
        metadata: Additional data (user_email, uploaded_at, etc.)
    """
    comparison_file = STORAGE_DIR / f"{original_id}_to_{refined_id}.json"
    
    data = {
        "original_audit_id": original_id,
        "refined_audit_id": refined_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        **metadata
    }
    
    with open(comparison_file, 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)

def get_comparison(original_id: str, refined_id: str) -> Optional[Dict]:
    """
    Retrieve comparison metadata.
    
    Returns:
        Comparison dict or None if not found
    """
    comparison_file = STORAGE_DIR / f"{original_id}_to_{refined_id}.json"
    
    if comparison_file.exists():
        with open(comparison_file, 'r') as f:
            return json.load(f)
    
    return None

def list_comparisons_by_original(original_id: str) -> List[Dict]:
    """
    Find all refined audits derived from a given original audit.
    
    Useful for showing "refinement history" on original audit page.
    """
    comparisons = []
    
    for comp_file in STORAGE_DIR.glob(f"{original_id}_to_*.json"):
        with open(comp_file, 'r') as f:
            comparisons.append(json.load(f))
    
    return sorted(comparisons, key=lambda x: x.get("created_at", ""), reverse=True)

def list_user_comparisons(user_email: str) -> List[Dict]:
    """
    List all comparisons for a user (if email was provided).
    
    Beta feature for user dashboard.
    """
    comparisons = []
    
    for comp_file in STORAGE_DIR.glob("*.json"):
        with open(comp_file, 'r') as f:
            data = json.load(f)
            if data.get("user_email") == user_email:
                comparisons.append(data)
    
    return sorted(comparisons, key=lambda x: x.get("created_at", ""), reverse=True)
