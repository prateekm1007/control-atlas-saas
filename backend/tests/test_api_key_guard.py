"""
Tests for API key authentication guard.
Tests the auth module directly — no HTTP server needed.
"""
import os
import pytest
from unittest.mock import patch
from fastapi import HTTPException
from tos.security.auth import (
    verify_api_key,
    validate_upload,
    validate_batch_upload,
    OPEN_MODE,
)


class TestVerifyApiKey:

    def test_open_mode_allows_no_key(self):
        """When no keys are configured, any request passes."""
        with patch("tos.security.auth.OPEN_MODE", True):
            with patch("tos.security.auth.ALLOWED_KEYS", set()):
                result = verify_api_key(x_api_key=None)
                assert result is True

    def test_valid_key_passes(self):
        with patch("tos.security.auth.OPEN_MODE", False):
            with patch("tos.security.auth.ALLOWED_KEYS", {"test-key-123"}):
                result = verify_api_key(x_api_key="test-key-123")
                assert result is True

    def test_invalid_key_raises_401(self):
        with patch("tos.security.auth.OPEN_MODE", False):
            with patch("tos.security.auth.ALLOWED_KEYS", {"test-key-123"}):
                with pytest.raises(HTTPException) as exc_info:
                    verify_api_key(x_api_key="wrong-key")
                assert exc_info.value.status_code == 401

    def test_missing_key_raises_401(self):
        with patch("tos.security.auth.OPEN_MODE", False):
            with patch("tos.security.auth.ALLOWED_KEYS", {"test-key-123"}):
                with pytest.raises(HTTPException) as exc_info:
                    verify_api_key(x_api_key=None)
                assert exc_info.value.status_code == 401

    def test_empty_string_key_raises_401(self):
        with patch("tos.security.auth.OPEN_MODE", False):
            with patch("tos.security.auth.ALLOWED_KEYS", {"test-key-123"}):
                with pytest.raises(HTTPException) as exc_info:
                    verify_api_key(x_api_key="")
                assert exc_info.value.status_code == 401

    def test_multiple_keys_supported(self):
        with patch("tos.security.auth.OPEN_MODE", False):
            with patch("tos.security.auth.ALLOWED_KEYS", {"key-a", "key-b", "key-c"}):
                assert verify_api_key(x_api_key="key-a") is True
                assert verify_api_key(x_api_key="key-b") is True
                assert verify_api_key(x_api_key="key-c") is True
                with pytest.raises(HTTPException):
                    verify_api_key(x_api_key="key-d")


class TestValidateUpload:

    def test_pdb_file_accepted(self):
        from unittest.mock import MagicMock
        file = MagicMock()
        file.filename = "structure.pdb"
        validate_upload(file)  # Should not raise

    def test_cif_file_accepted(self):
        from unittest.mock import MagicMock
        file = MagicMock()
        file.filename = "structure.cif"
        validate_upload(file)

    def test_txt_file_rejected(self):
        from unittest.mock import MagicMock
        file = MagicMock()
        file.filename = "readme.txt"
        with pytest.raises(HTTPException) as exc_info:
            validate_upload(file)
        assert exc_info.value.status_code == 400

    def test_none_file_passes(self):
        validate_upload(None)  # Should not raise


class TestValidateBatchUpload:

    def test_zip_file_accepted(self):
        from unittest.mock import MagicMock
        file = MagicMock()
        file.filename = "structures.zip"
        validate_batch_upload(file)  # Should not raise

    def test_non_zip_rejected(self):
        from unittest.mock import MagicMock
        file = MagicMock()
        file.filename = "structures.tar.gz"
        with pytest.raises(HTTPException) as exc_info:
            validate_batch_upload(file)
        assert exc_info.value.status_code == 400

    def test_none_file_rejected(self):
        with pytest.raises(HTTPException) as exc_info:
            validate_batch_upload(None)
        assert exc_info.value.status_code == 400
