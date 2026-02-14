"""
TOSCANINI Station Authentication & Input Validation.
Separate API key from third-party credentials.
"""
import os
import logging
from fastapi import Header, HTTPException, UploadFile

logger = logging.getLogger("toscanini.security")

TOSCANINI_API_KEY = os.getenv("TOSCANINI_API_KEY", "")
UPLOAD_LIMIT_BYTES = 25 * 1024 * 1024  # 25 MB
ALLOWED_EXTENSIONS = {".pdb", ".cif", ".mmcif"}


def verify_api_key(x_api_key: str = Header(None, alias="X-API-Key")):
    """FastAPI dependency for mandatory API key verification."""
    if not TOSCANINI_API_KEY:
        logger.warning("TOSCANINI_API_KEY not set — running in OPEN MODE (dev only)")
        return True
    if x_api_key != TOSCANINI_API_KEY:
        logger.warning(f"Unauthorized access attempt")
        raise HTTPException(status_code=401, detail="Invalid or missing API key")
    return True


def validate_upload(file: UploadFile) -> None:
    """Validate file extension and size before processing."""
    if file is None:
        return
    name = (file.filename or "").lower()
    ext = "." + name.rsplit(".", 1)[-1] if "." in name else ""
    if ext not in ALLOWED_EXTENSIONS:
        raise HTTPException(status_code=400, detail=f"Unsupported file type: {ext}. Allowed: {ALLOWED_EXTENSIONS}")


async def enforce_size_limit(content: bytes) -> bytes:
    """Reject oversized uploads."""
    if len(content) > UPLOAD_LIMIT_BYTES:
        raise HTTPException(
            status_code=413,
            detail=f"File exceeds {UPLOAD_LIMIT_BYTES // (1024*1024)}MB limit ({len(content)} bytes received)"
        )
    return content
