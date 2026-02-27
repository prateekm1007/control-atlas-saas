"""
Toscanini Phase B2 — Job State Manager
Tracks refinement job states in Redis.

States: queued -> running -> success | failed | timeout
"""
import json
import time
import os
from datetime import datetime, timezone
from typing import Optional, Dict, List

REDIS_URL     = os.environ.get("REDIS_URL", "redis://redis:6379/0")
JOB_TTL       = 86400 * 7  # 7 days

STATE_QUEUED  = "queued"
STATE_RUNNING = "running"
STATE_SUCCESS = "success"
STATE_FAILED  = "failed"
STATE_TIMEOUT = "timeout"


def _r():
    import redis
    return redis.from_url(REDIS_URL, decode_responses=True)


def create_job(job_id: str, original_audit_id: str,
               protocol: str, user_email: Optional[str] = None) -> Dict:
    job = {
        "job_id":            job_id,
        "original_audit_id": original_audit_id,
        "protocol":          protocol,
        "user_email":        user_email,
        "state":             STATE_QUEUED,
        "created_at":        datetime.now(timezone.utc).isoformat(),
        "updated_at":        datetime.now(timezone.utc).isoformat(),
        "started_at":        None,
        "completed_at":      None,
        "refined_audit_id":  None,
        "comparison_url":    None,
        "error":             None,
        "logs":              None,
        "refinement_method": "managed_gpu"
    }
    try:
        r = _r()
        r.setex(f"job:{job_id}", JOB_TTL, json.dumps(job))
        r.lpush(f"audit_jobs:{original_audit_id}", job_id)
        r.expire(f"audit_jobs:{original_audit_id}", JOB_TTL)
    except Exception:
        pass  # Degrade gracefully if Redis unavailable
    return job


def update_job(job_id: str, **kwargs) -> Optional[Dict]:
    try:
        r = _r()
        raw = r.get(f"job:{job_id}")
        if not raw:
            return None
        job = json.loads(raw)
        job.update(kwargs)
        job["updated_at"] = datetime.now(timezone.utc).isoformat()
        if kwargs.get("state") == STATE_RUNNING and not job.get("started_at"):
            job["started_at"] = datetime.now(timezone.utc).isoformat()
        if kwargs.get("state") in (STATE_SUCCESS, STATE_FAILED, STATE_TIMEOUT):
            job["completed_at"] = datetime.now(timezone.utc).isoformat()
        r.setex(f"job:{job_id}", JOB_TTL, json.dumps(job))
        return job
    except Exception:
        return None


def get_job(job_id: str) -> Optional[Dict]:
    try:
        r = _r()
        raw = r.get(f"job:{job_id}")
        return json.loads(raw) if raw else None
    except Exception:
        return None


def get_jobs_for_audit(original_audit_id: str) -> List[Dict]:
    try:
        r = _r()
        job_ids = r.lrange(f"audit_jobs:{original_audit_id}", 0, -1)
        jobs = [get_job(jid) for jid in job_ids if get_job(jid)]
        return sorted(jobs, key=lambda x: x.get("created_at", ""), reverse=True)
    except Exception:
        return []


def get_job_duration(job: Dict) -> Optional[float]:
    if not job.get("started_at"):
        return None
    start = datetime.fromisoformat(job["started_at"])
    end   = datetime.fromisoformat(job["completed_at"]) if job.get("completed_at") else datetime.now(timezone.utc)
    return round((end - start).total_seconds(), 1)
