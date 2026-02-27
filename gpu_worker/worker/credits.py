"""
Toscanini Phase B2 — Credit System
Tracks GPU compute credits per user.

Phase 1: File-based (no auth required)
Phase 2: Redis + user accounts

Credit model:
- Beta: 10 free managed runs per email
- Anonymous: 3 free managed runs per IP
- Cost tracking: GPU minutes per job
"""
import json
import os
import time
from pathlib import Path
from typing import Optional, Dict
from datetime import datetime, timezone

CREDITS_DIR   = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "credits"
CREDITS_DIR.mkdir(parents=True, exist_ok=True)

# Beta allocation
BETA_CREDITS_EMAIL = 10
BETA_CREDITS_ANON  = 3

# Cost per protocol (GPU minutes)
PROTOCOL_COST = {
    "openmm":  5,   # ~5 GPU minutes
    "rosetta": 3,   # ~3 GPU minutes
    "both":    10,  # ~10 GPU minutes
    "loop":    15,  # ~15 GPU minutes
}


def _credit_file(identifier: str) -> Path:
    """Get credit file path for a user identifier."""
    safe_id = identifier.replace("@", "_at_").replace(".", "_").replace("/", "_")
    return CREDITS_DIR / f"{safe_id}.json"


def get_credits(identifier: str) -> Dict:
    """
    Get credit balance for a user.

    Args:
        identifier: Email address or IP address

    Returns:
        Credit record dict
    """
    cf = _credit_file(identifier)

    if cf.exists():
        with open(cf, "r") as f:
            return json.load(f)

    # New user — initialize with beta allocation
    is_email    = "@" in identifier
    max_credits = BETA_CREDITS_EMAIL if is_email else BETA_CREDITS_ANON

    record = {
        "identifier":    identifier,
        "is_email":      is_email,
        "credits_total": max_credits,
        "credits_used":  0,
        "credits_remaining": max_credits,
        "gpu_minutes_used":  0,
        "jobs_submitted":    0,
        "jobs_completed":    0,
        "jobs_failed":       0,
        "created_at":    datetime.now(timezone.utc).isoformat(),
        "updated_at":    datetime.now(timezone.utc).isoformat(),
        "tier":          "beta"
    }

    _save_credits(identifier, record)
    return record


def _save_credits(identifier: str, record: Dict) -> None:
    """Persist credit record."""
    record["updated_at"] = datetime.now(timezone.utc).isoformat()
    cf = _credit_file(identifier)
    with open(cf, "w") as f:
        json.dump(record, f, indent=2, sort_keys=True)


def check_credits(identifier: str, protocol: str) -> Dict:
    """
    Check if user has sufficient credits for a protocol.

    Returns:
        dict with allowed, credits_remaining, cost, reason
    """
    record  = get_credits(identifier)
    cost    = PROTOCOL_COST.get(protocol, 5)
    allowed = record["credits_remaining"] > 0

    return {
        "allowed":           allowed,
        "credits_remaining": record["credits_remaining"],
        "credits_total":     record["credits_total"],
        "cost":              cost,
        "reason": (
            "Sufficient credits" if allowed
            else f"No credits remaining. Beta allocation: {record['credits_total']} runs."
        )
    }


def deduct_credits(identifier: str, protocol: str,
                   job_id: str, gpu_minutes: float = 0) -> Dict:
    """
    Deduct credits after job submission.

    Args:
        identifier:  User email or IP
        protocol:    Job protocol
        job_id:      Job ID for audit trail
        gpu_minutes: Actual GPU minutes used (0 = estimated)

    Returns:
        Updated credit record
    """
    record = get_credits(identifier)
    cost   = PROTOCOL_COST.get(protocol, 5)

    record["credits_used"]      += 1
    record["credits_remaining"]  = max(0, record["credits_remaining"] - 1)
    record["gpu_minutes_used"]  += gpu_minutes if gpu_minutes > 0 else cost
    record["jobs_submitted"]    += 1

    # Add job to history
    if "job_history" not in record:
        record["job_history"] = []

    record["job_history"].append({
        "job_id":     job_id,
        "protocol":   protocol,
        "cost":       cost,
        "submitted_at": datetime.now(timezone.utc).isoformat()
    })

    # Keep last 50 jobs only
    record["job_history"] = record["job_history"][-50:]

    _save_credits(identifier, record)
    return record


def record_job_completion(identifier: str, job_id: str,
                          success: bool, gpu_minutes: float = 0) -> None:
    """Update credit record after job completes."""
    if not identifier:
        return

    record = get_credits(identifier)

    if success:
        record["jobs_completed"] += 1
    else:
        record["jobs_failed"] += 1

    if gpu_minutes > 0:
        # Update with actual GPU time
        record["gpu_minutes_used"] = round(
            record.get("gpu_minutes_used", 0) + gpu_minutes, 2
        )

    _save_credits(identifier, record)


def get_usage_stats(identifier: str) -> Dict:
    """Get formatted usage statistics for a user."""
    record = get_credits(identifier)
    return {
        "identifier":        record["identifier"],
        "tier":              record["tier"],
        "credits_remaining": record["credits_remaining"],
        "credits_total":     record["credits_total"],
        "credits_used":      record["credits_used"],
        "gpu_minutes_used":  record.get("gpu_minutes_used", 0),
        "jobs_submitted":    record.get("jobs_submitted", 0),
        "jobs_completed":    record.get("jobs_completed", 0),
        "jobs_failed":       record.get("jobs_failed", 0),
        "member_since":      record.get("created_at", "unknown")[:10]
    }
