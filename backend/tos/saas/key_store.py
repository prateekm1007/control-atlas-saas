"""
API Key Store — SQLite-backed key management with quota enforcement.
Designed for easy migration to PostgreSQL later.
"""
import sqlite3
import hashlib
import secrets
from datetime import date, datetime, timezone
from typing import Optional
from dataclasses import dataclass
from pathlib import Path

DB_PATH = Path(__file__).parent / "keys.db"

@dataclass
class APIKey:
    key_hash: str
    tier: str
    active: bool
    created_at: str
    last_used: Optional[str]

TIER_LIMITS = {
    "free":       {"daily": 10,     "batch_max": 25,  "gpu_runs": 3},
    "pro":        {"daily": 500,    "batch_max": 500,  "gpu_runs": 50},
    "enterprise": {"daily": 999999, "batch_max": 999999, "gpu_runs": 999999},
}

# ── Week B.3: GPU credit allocation by tier ───────────────────────────────────
GPU_TIER_CREDITS = {
    "free":       3,
    "pro":        50,
    "enterprise": 999999,
}

def _get_conn():
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    """Initialize database schema."""
    conn = _get_conn()
    conn.execute("""
        CREATE TABLE IF NOT EXISTS api_keys (
            key_hash TEXT PRIMARY KEY,
            tier TEXT NOT NULL,
            active INTEGER NOT NULL DEFAULT 1,
            created_at TEXT NOT NULL,
            last_used TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS quota_usage (
            key_hash TEXT NOT NULL,
            date TEXT NOT NULL,
            count INTEGER NOT NULL DEFAULT 0,
            PRIMARY KEY (key_hash, date),
            FOREIGN KEY (key_hash) REFERENCES api_keys(key_hash)
        )
    """)
    conn.commit()
    conn.close()

def hash_key(key: str) -> str:
    """SHA256 hash of API key."""
    return hashlib.sha256(key.encode()).hexdigest()

def create_key(tier: str = "free") -> str:
    """Generate new API key and store hash."""
    if tier not in TIER_LIMITS:
        raise ValueError(f"Invalid tier: {tier}")
    
    key = f"tos_{secrets.token_hex(32)}"
    key_hash_val = hash_key(key)
    
    conn = _get_conn()
    conn.execute(
        "INSERT INTO api_keys (key_hash, tier, created_at) VALUES (?, ?, ?)",
        (key_hash_val, tier, datetime.now(timezone.utc).isoformat())
    )
    conn.commit()
    conn.close()
    
    return key  # Return plaintext only once

def get_key(key_hash_val: str) -> Optional[APIKey]:
    """Retrieve key metadata."""
    conn = _get_conn()
    row = conn.execute(
        "SELECT * FROM api_keys WHERE key_hash = ?",
        (key_hash_val,)
    ).fetchone()
    conn.close()
    
    if not row:
        return None
    
    return APIKey(
        key_hash=row["key_hash"],
        tier=row["tier"],
        active=bool(row["active"]),
        created_at=row["created_at"],
        last_used=row["last_used"]
    )

def check_quota(key_hash_val: str) -> tuple[bool, int, int]:
    """
    Check if key is under daily quota.
    Returns: (allowed, used, limit)
    """
    key = get_key(key_hash_val)
    if not key or not key.active:
        return (False, 0, 0)
    
    tier = key.tier
    limit = TIER_LIMITS[tier]["daily"]
    today = date.today().isoformat()
    
    conn = _get_conn()
    row = conn.execute(
        "SELECT count FROM quota_usage WHERE key_hash = ? AND date = ?",
        (key_hash_val, today)
    ).fetchone()
    conn.close()
    
    used = row["count"] if row else 0
    allowed = used < limit
    
    return (allowed, used, limit)

def increment_usage(key_hash_val: str) -> None:
    """Atomically increment daily usage counter."""
    today = date.today().isoformat()
    
    conn = _get_conn()
    conn.execute(
        """
        INSERT INTO quota_usage (key_hash, date, count)
        VALUES (?, ?, 1)
        ON CONFLICT(key_hash, date)
        DO UPDATE SET count = count + 1
        """,
        (key_hash_val, today)
    )
    conn.execute(
        "UPDATE api_keys SET last_used = ? WHERE key_hash = ?",
        (datetime.now(timezone.utc).isoformat(), key_hash_val)
    )
    conn.commit()
    conn.close()

def revoke_key(key_hash_val: str) -> bool:
    """Revoke (deactivate) an API key."""
    conn = _get_conn()
    cursor = conn.execute(
        "UPDATE api_keys SET active = 0 WHERE key_hash = ?",
        (key_hash_val,)
    )
    conn.commit()
    conn.close()
    return cursor.rowcount > 0


def get_gpu_allocation_for_key(key_hash_val: str) -> int:
    """
    Return the GPU run allocation for a given API key.

    Used by /refinement/submit to set the credit ceiling
    based on the key tier rather than anonymous email/IP.

    Returns:
        int: number of GPU runs allowed (999999 = unlimited)
        Returns 3 (anon default) if key not found or inactive.
    """
    key = get_key(key_hash_val)
    if not key or not key.active:
        return GPU_TIER_CREDITS["free"]
    return GPU_TIER_CREDITS.get(key.tier, GPU_TIER_CREDITS["free"])


def get_tier_for_key(key_hash_val: str) -> str:
    """Return tier string for a given key hash. Returns 'free' if not found."""
    key = get_key(key_hash_val)
    if not key or not key.active:
        return "free"
    return key.tier

# Initialize DB on import
init_db()
