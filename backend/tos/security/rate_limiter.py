"""
Toscanini Phase B1 — Rate Limiting
Simple file-based rate limiter for API endpoints.
Phase 2: migrate to Redis with sliding window.
"""
import json
import time
from pathlib import Path
from typing import Dict, Tuple
from collections import defaultdict

RATE_LIMIT_DIR = Path("/app/data/rate_limits")
RATE_LIMIT_DIR.mkdir(parents=True, exist_ok=True)

# Rate limit rules per endpoint
RATE_LIMITS = {
    "/ingest": {"requests": 20, "window_seconds": 3600},           # 20 uploads/hour
    "/refinement/callback": {"requests": 5, "window_seconds": 3600}, # 5 callbacks/hour per IP
    "/refinement/submit": {"requests": 10, "window_seconds": 3600},  # 10 submits/hour
}


def _get_rate_file(ip: str, endpoint: str) -> Path:
    """Get rate limit tracking file for an IP + endpoint combo."""
    safe_ip = ip.replace(".", "_").replace(":", "_")
    safe_ep = endpoint.replace("/", "_")
    return RATE_LIMIT_DIR / f"{safe_ip}{safe_ep}.json"


def check_rate_limit(ip: str, endpoint: str) -> Tuple[bool, Dict]:
    """
    Check if request is within rate limit.
    
    Returns:
        (allowed: bool, info: dict with limit details)
    """
    limit_config = RATE_LIMITS.get(endpoint, {"requests": 100, "window_seconds": 3600})
    max_requests = limit_config["requests"]
    window = limit_config["window_seconds"]
    
    rate_file = _get_rate_file(ip, endpoint)
    now = time.time()
    
    # Load existing rate data
    if rate_file.exists():
        with open(rate_file, 'r') as f:
            data = json.load(f)
    else:
        data = {"requests": [], "endpoint": endpoint, "ip": ip}
    
    # Remove expired requests outside window
    data["requests"] = [t for t in data["requests"] if now - t < window]
    
    # Check limit
    current_count = len(data["requests"])
    allowed = current_count < max_requests
    
    if allowed:
        # Record this request
        data["requests"].append(now)
        with open(rate_file, 'w') as f:
            json.dump(data, f)
    
    return allowed, {
        "allowed": allowed,
        "current_count": current_count + (1 if allowed else 0),
        "max_requests": max_requests,
        "window_seconds": window,
        "reset_in_seconds": int(window - (now - min(data["requests"]))) if data["requests"] else window
    }
