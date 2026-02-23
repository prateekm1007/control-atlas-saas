"""
TOSCANINI Station Authentication & Input Validation.
Supports single key (TOSCANINI_API_KEY) or multiple keys (TOSCANINI_API_KEYS).
"""
import os
import logging
from fastapi import Header, HTTPException, UploadFile
from tos.saas.key_store import hash_key, get_key, check_quota, increment_usage

logger = logging.getLogger("toscanini.security")

# Support both single and multi-key configurations
_single_key = os.getenv("TOSCANINI_API_KEY", "")
_multi_keys = os.getenv("TOSCANINI_API_KEYS", "")

ALLOWED_KEYS: set[str] = set()
if _multi_keys:
    ALLOWED_KEYS = {k.strip() for k in _multi_keys.split(",") if k.strip()}
elif _single_key:
    ALLOWED_KEYS = {_single_key}

OPEN_MODE = len(ALLOWED_KEYS) == 0

if OPEN_MODE:
    logger.warning("No API keys configured — running in OPEN MODE (dev only)")
else:
    logger.info(f"API key guard active: {len(ALLOWED_KEYS)} key(s) configured")

UPLOAD_LIMIT_BYTES = 25 * 1024 * 1024  # 25 MB
ALLOWED_EXTENSIONS = {".pdb", ".cif", ".mmcif"}
BATCH_ALLOWED_EXTENSIONS = {".zip"}


def verify_api_key(x_api_key: str = Header(None, alias="X-API-Key")):
    """FastAPI dependency for API key verification with quota enforcement."""
    # Legacy OPEN_MODE (dev only)
    if OPEN_MODE:
        return True
    
    # Check against static keys first (backward compat)
    if x_api_key and x_api_key in ALLOWED_KEYS:
        return True
    
    # DB-backed multi-tenant key validation
    if not x_api_key:
        logger.warning("Missing API key")
        raise HTTPException(status_code=401, detail="Missing API key")
    
    key_hash_val = hash_key(x_api_key)
    key = get_key(key_hash_val)
    
    if not key:
        logger.warning(f"Invalid API key: {key_hash_val[:16]}...")
        raise HTTPException(status_code=401, detail="Invalid API key")
    
    if not key.active:
        logger.warning(f"Revoked API key: {key_hash_val[:16]}...")
        raise HTTPException(status_code=401, detail="API key has been revoked")
    
    # Quota check
    allowed, used, limit = check_quota(key_hash_val)
    if not allowed:
        logger.warning(f"Quota exceeded for key: {used}/{limit}")
        raise HTTPException(
            status_code=429,
            detail=f"Daily quota exceeded: {used}/{limit} requests used"
        )
    
    # Increment usage atomically
    increment_usage(key_hash_val)
    
    logger.info(f"Authorized request: tier={key.tier}, usage={used+1}/{limit}")
    return True


def validate_upload(file: UploadFile) -> None:
    """Validate file extension before processing."""
    if file is None:
        return
    name = (file.filename or "").lower()
    ext = "." + name.rsplit(".", 1)[-1] if "." in name else ""
    if ext not in ALLOWED_EXTENSIONS:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported file type: {ext}. Allowed: {ALLOWED_EXTENSIONS}",
        )


def validate_batch_upload(file: UploadFile) -> None:
    """Validate batch upload is a ZIP file."""
    if file is None:
        raise HTTPException(status_code=400, detail="ZIP file required")
    name = (file.filename or "").lower()
    ext = "." + name.rsplit(".", 1)[-1] if "." in name else ""
    if ext not in BATCH_ALLOWED_EXTENSIONS:
        raise HTTPException(
            status_code=400,
            detail=f"Batch endpoint requires .zip file, got: {ext}",
        )


async def enforce_size_limit(content: bytes) -> bytes:
    """Reject oversized uploads."""
    if len(content) > UPLOAD_LIMIT_BYTES:
        raise HTTPException(
            status_code=413,
            detail=f"File exceeds {UPLOAD_LIMIT_BYTES // (1024*1024)}MB limit ({len(content)} bytes received)",
        )
    return content
