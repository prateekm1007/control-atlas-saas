"""
Toscanini Phase B2 — Cost Tracker
Tracks actual GPU compute costs per job.

Estimates based on AWS g4dn.xlarge pricing (~$0.52/hr on-demand).
"""
import json
import os
import time
from pathlib import Path
from datetime import datetime, timezone
from typing import Optional, Dict

COST_DIR     = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "costs"
COST_DIR.mkdir(parents=True, exist_ok=True)

# AWS g4dn.xlarge: $0.526/hr = $0.00877/min
GPU_COST_PER_MINUTE = 0.00877

# Protocol expected durations (minutes)
EXPECTED_DURATION = {
    "openmm":  5.0,
    "rosetta": 3.0,
    "both":    10.0,
    "loop":    15.0,
}


def estimate_cost(protocol: str) -> Dict:
    """
    Estimate job cost before execution.

    Returns:
        dict with estimated_minutes, estimated_usd
    """
    minutes = EXPECTED_DURATION.get(protocol, 5.0)
    usd     = round(minutes * GPU_COST_PER_MINUTE, 4)

    return {
        "protocol":          protocol,
        "estimated_minutes": minutes,
        "estimated_usd":     usd,
        "gpu_instance":      "g4dn.xlarge",
        "note":              "Estimate only. Actual cost depends on structure complexity."
    }


def record_job_cost(job_id: str, protocol: str,
                    start_time: float, end_time: float,
                    success: bool) -> Dict:
    """
    Record actual job cost after completion.

    Args:
        job_id:     Job identifier
        protocol:   Protocol used
        start_time: Unix timestamp of job start
        end_time:   Unix timestamp of job end
        success:    Whether job succeeded

    Returns:
        Cost record dict
    """
    duration_seconds = end_time - start_time
    duration_minutes = round(duration_seconds / 60, 2)
    actual_usd       = round(duration_minutes * GPU_COST_PER_MINUTE, 4)
    estimated        = estimate_cost(protocol)

    record = {
        "job_id":            job_id,
        "protocol":          protocol,
        "success":           success,
        "duration_seconds":  round(duration_seconds, 1),
        "duration_minutes":  duration_minutes,
        "actual_usd":        actual_usd,
        "estimated_usd":     estimated["estimated_usd"],
        "cost_accuracy_pct": round(
            actual_usd / estimated["estimated_usd"] * 100, 1
        ) if estimated["estimated_usd"] > 0 else 0,
        "recorded_at":       datetime.now(timezone.utc).isoformat()
    }

    # Save cost record
    cost_file = COST_DIR / f"{job_id}.json"
    with open(cost_file, "w") as f:
        json.dump(record, f, indent=2, sort_keys=True)

    return record


def get_total_costs(days: int = 30) -> Dict:
    """
    Aggregate total costs across all jobs.

    Args:
        days: Number of days to look back

    Returns:
        Aggregated cost summary
    """
    total_usd      = 0.0
    total_minutes  = 0.0
    total_jobs     = 0
    successful     = 0
    failed_jobs    = 0
    by_protocol    = {}

    cutoff = time.time() - (days * 86400)

    for cost_file in COST_DIR.glob("*.json"):
        try:
            with open(cost_file, "r") as f:
                rec = json.load(f)

            # Check if within time window
            rec_time = datetime.fromisoformat(
                rec.get("recorded_at", "2020-01-01T00:00:00+00:00")
            ).timestamp()

            if rec_time < cutoff:
                continue

            total_usd     += rec.get("actual_usd", 0)
            total_minutes += rec.get("duration_minutes", 0)
            total_jobs    += 1

            if rec.get("success"):
                successful += 1
            else:
                failed_jobs += 1

            proto = rec.get("protocol", "unknown")
            if proto not in by_protocol:
                by_protocol[proto] = {"jobs": 0, "usd": 0.0, "minutes": 0.0}
            by_protocol[proto]["jobs"]    += 1
            by_protocol[proto]["usd"]     += rec.get("actual_usd", 0)
            by_protocol[proto]["minutes"] += rec.get("duration_minutes", 0)

        except Exception:
            continue

    return {
        "period_days":        days,
        "total_jobs":         total_jobs,
        "successful_jobs":    successful,
        "failed_jobs":        failed_jobs,
        "success_rate_pct":   round(successful / total_jobs * 100, 1) if total_jobs > 0 else 0,
        "total_usd":          round(total_usd, 4),
        "total_gpu_minutes":  round(total_minutes, 1),
        "avg_cost_per_job":   round(total_usd / total_jobs, 4) if total_jobs > 0 else 0,
        "by_protocol":        by_protocol
    }
