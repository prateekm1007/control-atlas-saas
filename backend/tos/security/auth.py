"""
TOSCANINI Station Authentication & Input Validation.
Supports single key (TOSCANINI_API_KEY) or multiple keys (TOSCANINI_API_KEYS).
"""
import os
import logging
from fastapi import Header, HTTPException, UploadFile

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
    """FastAPI dependency for API key verification."""
    if OPEN_MODE:
        return True
    if not x_api_key or x_api_key not in ALLOWED_KEYS:
        logger.warning("Unauthorized access attempt")
        raise HTTPException(status_code=401, detail="Invalid or missing API key")
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
