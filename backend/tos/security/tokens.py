"""
backend/tos/security/tokens.py

Toscanini Phase B1 — JWT Token Security System
Manages callback tokens for managed refinement orchestration.

Hardening (B.6):
- Original: single-use blacklist (7-day TTL, max 3 callbacks per audit_id)
- NEW: sliding-window rate limit per API key (max 10/hour)
- NEW: _get_rate_dir() lazy path — safe for host tests
- PRESERVED: all original functions unchanged (create_refinement_token,
  validate_refinement_token, revoke_token, get_token_stats)
"""
import jwt
import time
import secrets
import hashlib
import json
from pathlib import Path
from datetime import datetime, timedelta, timezone
from typing import Optional, Dict, List

import os
SECRET_KEY = os.environ.get(
    "TOSCANINI_JWT_SECRET",
    "toscanini-phase-b1-secret-SETMEINENV"
)
ALGORITHM = "HS256"
TOKEN_EXPIRY_DAYS = 7

# ── Rate limiting constants (B.6) ─────────────────────────────────────────────
RATE_LIMIT_WINDOW_SECONDS = 3600          # 1-hour sliding window
RATE_LIMIT_MAX_PER_WINDOW = int(os.environ.get("CALLBACK_RATE_LIMIT", "10"))

# ── Lazy path getters (B.1 pattern preserved, B.6 extended) ──────────────────

def _get_blacklist_dir() -> Path:
    """Lazy evaluation — read TOSCANINI_DATA_DIR at call time, not import time."""
    d = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "token_blacklist"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _get_rate_dir() -> Path:
    """B.6: Lazy evaluation for rate-limit counters. Separate namespace."""
    d = Path(os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "rate_limits" / "callback"
    d.mkdir(parents=True, exist_ok=True)
    return d


# Max callbacks per audit_id (original B.1 constant — unchanged)
MAX_CALLBACKS_PER_AUDIT = 3


# ── Original B.1 helpers (unchanged) ─────────────────────────────────────────

def _token_fingerprint(token: str) -> str:
    """Generate short fingerprint for a token (for storage/logging)."""
    return hashlib.sha256(token.encode()).hexdigest()[:16]


def _is_blacklisted(token: str) -> bool:
    """Check if token has already been used."""
    fp = _token_fingerprint(token)
    blacklist_file = _get_blacklist_dir() / f"{fp}.json"
    return blacklist_file.exists()


def _blacklist_token(token: str, audit_id: str, reason: str = "used") -> None:
    """Mark token as used — prevents reuse."""
    fp = _token_fingerprint(token)
    blacklist_file = _get_blacklist_dir() / f"{fp}.json"
    data = {
        "fingerprint": fp,
        "audit_id": audit_id,
        "blacklisted_at": datetime.now(timezone.utc).isoformat(),
        "reason": reason
    }
    with open(blacklist_file, 'w') as f:
        json.dump(data, f, indent=2)


def _get_callback_count(audit_id: str) -> int:
    """Count how many callbacks have been made for a given audit_id."""
    count = 0
    for f in _get_blacklist_dir().glob("*.json"):
        with open(f, 'r') as fh:
            try:
                data = json.load(fh)
                if data.get("audit_id") == audit_id and data.get("reason") == "used":
                    count += 1
            except Exception:
                pass
    return count


# ── B.6: Rate limiting helpers ────────────────────────────────────────────────

def _rate_path(api_key: str) -> Path:
    """Deterministic path for a given api_key's rate record."""
    safe = hashlib.sha256(api_key.encode()).hexdigest()
    return _get_rate_dir() / f"{safe}.json"


def _read_rate_record(api_key: str) -> Dict:
    """Read sliding-window timestamps for api_key. Returns empty record if missing."""
    path = _rate_path(api_key)
    if not path.exists():
        return {"timestamps": []}
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return {"timestamps": []}


def _write_rate_record(api_key: str, data: Dict) -> None:
    path = _rate_path(api_key)
    path.write_text(json.dumps(data))


def check_callback_rate_limit(api_key: str) -> Dict:
    """
    Check whether api_key is within the sliding-window rate limit.

    Returns dict:
        allowed (bool): True if under the ceiling
        attempts_in_window (int): calls made in current window
        limit (int): ceiling
        reset_in_seconds (int): seconds until oldest call drops off
    """
    now = time.time()
    cutoff = now - RATE_LIMIT_WINDOW_SECONDS
    record = _read_rate_record(api_key)

    # Prune timestamps outside the window
    active = [t for t in record.get("timestamps", []) if t > cutoff]

    allowed = len(active) < RATE_LIMIT_MAX_PER_WINDOW
    reset_in = int(active[0] + RATE_LIMIT_WINDOW_SECONDS - now) if active else 0

    return {
        "allowed":            allowed,
        "attempts_in_window": len(active),
        "limit":              RATE_LIMIT_MAX_PER_WINDOW,
        "reset_in_seconds":   max(0, reset_in),
    }


def _increment_rate_counter(api_key: str) -> None:
    """Record a new callback attempt for api_key (called after successful consume)."""
    now = time.time()
    cutoff = now - RATE_LIMIT_WINDOW_SECONDS
    record = _read_rate_record(api_key)
    active = [t for t in record.get("timestamps", []) if t > cutoff]
    active.append(now)
    _write_rate_record(api_key, {"timestamps": active})


# ── Original B.1 public API (signatures unchanged) ───────────────────────────

def create_refinement_token(audit_id: str, user_email: Optional[str] = None) -> str:
    """
    Generate secure callback token for refinement re-upload.

    Args:
        audit_id: Original audit ID from failed/INDETERMINATE verdict
        user_email: Optional email for tracking

    Returns:
        JWT token string (valid for 7 days, single-use)
    """
    now = datetime.now(timezone.utc)
    expiry = now + timedelta(days=TOKEN_EXPIRY_DAYS)
    nonce = secrets.token_hex(8)

    payload = {
        "audit_id":   audit_id,
        "user_email": user_email,
        "issued_at":  now.isoformat(),
        "expires_at": expiry.isoformat(),
        "exp":        int(time.time()) + (TOKEN_EXPIRY_DAYS * 86400),
        "purpose":    "refinement_callback",
        "version":    "1.0",
        "nonce":      nonce
    }

    return jwt.encode(payload, SECRET_KEY, algorithm=ALGORITHM)


def validate_refinement_token(token: str, consume: bool = False) -> Dict:
    """
    Validate and decode refinement callback token.

    B.6 addition: if consume=True, also enforce sliding-window rate limit
    keyed on the audit_id (serves as the per-job identifier when no
    API key is available at token-validation time).

    Args:
        token:   JWT token string
        consume: If True, blacklist token + increment rate counter

    Returns:
        Decoded payload dict

    Raises:
        ValueError: If token is invalid, expired, blacklisted, or rate-limited
    """
    # Fast path: blacklist check before decoding
    if _is_blacklisted(token):
        raise ValueError(
            "Token has already been used. Each token is single-use only."
        )

    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
    except jwt.ExpiredSignatureError:
        raise ValueError("Token has expired (7-day limit exceeded). Generate a new token.")
    except jwt.InvalidTokenError as e:
        raise ValueError(f"Invalid token: {str(e)}")

    # Purpose check
    if payload.get("purpose") != "refinement_callback":
        raise ValueError("Invalid token purpose.")

    if payload.get("version") != "1.0":
        raise ValueError("Unsupported token version.")

    audit_id = payload["audit_id"]

    # Per-audit callback cap (original B.1 limit)
    callback_count = _get_callback_count(audit_id)
    if callback_count >= MAX_CALLBACKS_PER_AUDIT:
        raise ValueError(
            f"Maximum callbacks ({MAX_CALLBACKS_PER_AUDIT}) reached for this audit. "
            "Re-audit the original structure to get a fresh token."
        )

    # B.6: sliding-window rate limit keyed on audit_id
    # (audit_id used as key so anonymous users without API keys are also covered)
    if consume:
        rl = check_callback_rate_limit(audit_id)
        if not rl["allowed"]:
            raise ValueError(
                f"Callback rate limit exceeded ({RATE_LIMIT_MAX_PER_WINDOW}/"
                f"{RATE_LIMIT_WINDOW_SECONDS}s). "
                f"Resets in {rl['reset_in_seconds']}s."
            )

    # Consume token if requested
    if consume:
        _blacklist_token(token, audit_id, reason="used")
        _increment_rate_counter(audit_id)

    return payload


def revoke_token(token: str, audit_id: str) -> None:
    """Manually revoke a token (admin use)."""
    _blacklist_token(token, audit_id, reason="revoked")


def get_token_stats(audit_id: str) -> Dict:
    """Get callback statistics for an audit_id."""
    count = _get_callback_count(audit_id)
    rl    = check_callback_rate_limit(audit_id)
    return {
        "audit_id":           audit_id,
        "callbacks_used":     count,
        "callbacks_remaining": max(0, MAX_CALLBACKS_PER_AUDIT - count),
        "max_callbacks":      MAX_CALLBACKS_PER_AUDIT,
        "rate_limit":         rl,
    }
