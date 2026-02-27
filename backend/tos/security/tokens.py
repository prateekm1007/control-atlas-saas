"""
Toscanini Phase B1 — JWT Token Security System
Manages callback tokens for managed refinement orchestration.

Hardening (Week 4):
- Single-use enforcement via file-based blacklist
- Rate limiting per audit_id (max 3 callback attempts)
- Token fingerprinting for audit trail
"""
import jwt
import time
import secrets
import hashlib
import json
from pathlib import Path
from datetime import datetime, timedelta
from typing import Optional, Dict

# Secret key (production: load from environment variable)
import os
SECRET_KEY = os.environ.get(
    "TOSCANINI_JWT_SECRET",
    "toscanini-phase-b1-secret-" + secrets.token_hex(16)
)
ALGORITHM = "HS256"
TOKEN_EXPIRY_DAYS = 7

# Blacklist storage (Phase 1: file-based. Phase 2: Redis)
import os as _os
def _get_blacklist_dir() -> Path:
    """Lazy evaluation — read TOSCANINI_DATA_DIR at call time, not import time."""
    d = Path(_os.environ.get("TOSCANINI_DATA_DIR", "/app/data")) / "token_blacklist"
    d.mkdir(parents=True, exist_ok=True)
    return d

# Rate limit: max callbacks per original audit_id
MAX_CALLBACKS_PER_AUDIT = 3


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
        "blacklisted_at": datetime.utcnow().isoformat(),
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


def create_refinement_token(audit_id: str, user_email: Optional[str] = None) -> str:
    """
    Generate secure callback token for refinement re-upload.
    
    Args:
        audit_id: Original audit ID from failed/INDETERMINATE verdict
        user_email: Optional email for tracking
    
    Returns:
        JWT token string (valid for 7 days, single-use)
    """
    now = datetime.utcnow()
    expiry = now + timedelta(days=TOKEN_EXPIRY_DAYS)
    nonce = secrets.token_hex(8)  # Prevents token reuse even with same payload

    payload = {
        "audit_id": audit_id,
        "user_email": user_email,
        "issued_at": now.isoformat(),
        "expires_at": expiry.isoformat(),
        "exp": int(time.time()) + (TOKEN_EXPIRY_DAYS * 86400),
        "purpose": "refinement_callback",
        "version": "1.0",
        "nonce": nonce  # Unique per token
    }

    return jwt.encode(payload, SECRET_KEY, algorithm=ALGORITHM)


def validate_refinement_token(token: str, consume: bool = False) -> Dict:
    """
    Validate and decode refinement callback token.
    
    Args:
        token: JWT token string
        consume: If True, blacklist token after validation (single-use enforcement)
    
    Returns:
        Decoded payload dict
    
    Raises:
        ValueError: If token is invalid, expired, blacklisted, or rate-limited
    """
    # Check blacklist first (before decoding — fast path)
    if _is_blacklisted(token):
        raise ValueError("Token has already been used. Each token is single-use only.")

    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
    except jwt.ExpiredSignatureError:
        raise ValueError("Token has expired (7-day limit exceeded). Generate a new token.")
    except jwt.InvalidTokenError as e:
        raise ValueError(f"Invalid token: {str(e)}")

    # Verify purpose
    if payload.get("purpose") != "refinement_callback":
        raise ValueError("Invalid token purpose.")

    # Verify version
    if payload.get("version") != "1.0":
        raise ValueError("Unsupported token version.")

    # Rate limit check
    audit_id = payload["audit_id"]
    callback_count = _get_callback_count(audit_id)
    if callback_count >= MAX_CALLBACKS_PER_AUDIT:
        raise ValueError(
            f"Maximum callbacks ({MAX_CALLBACKS_PER_AUDIT}) reached for this audit. "
            f"Re-audit the original structure to get a fresh token."
        )

    # Consume token if requested
    if consume:
        _blacklist_token(token, audit_id, reason="used")

    return payload


def revoke_token(token: str, audit_id: str) -> None:
    """Manually revoke a token (admin use)."""
    _blacklist_token(token, audit_id, reason="revoked")


def get_token_stats(audit_id: str) -> Dict:
    """Get callback statistics for an audit_id."""
    count = _get_callback_count(audit_id)
    return {
        "audit_id": audit_id,
        "callbacks_used": count,
        "callbacks_remaining": max(0, MAX_CALLBACKS_PER_AUDIT - count),
        "max_callbacks": MAX_CALLBACKS_PER_AUDIT
    }
