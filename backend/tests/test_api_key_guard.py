"""
Tests for API key authentication guard.
Tests the auth module directly — no HTTP server needed.

WP-2 additions:
- Open mode is never implicit (requires TOSCANINI_ALLOW_OPEN_MODE=1).
- Fail-closed behavior when no keys are configured (503 on protected routes).
- Parametrized coverage over every public API route: 401 without a key on
  protected routes, non-auth expected status with a valid key.
"""
import os

# ── WP-2: auth.py reads env at import time — configure BEFORE importing ─────
# Keys configured so the app under test runs in enforced mode.
os.environ["TOSCANINI_API_KEYS"] = "test-key-xyz"
os.environ.setdefault("TOSCANINI_DATA_DIR", os.path.abspath(
    os.path.join(os.path.dirname(__file__), ".tmp-guard-data")))
os.environ.pop("TOSCANINI_ALLOW_OPEN_MODE", None)  # open mode must be off by default

import tempfile
import pytest
from unittest.mock import patch, MagicMock
from fastapi import HTTPException, Depends, FastAPI
from fastapi.testclient import TestClient
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


# ═══════════════════════════════════════════════════════════════════════════
# WP-2: open mode is never implicit + fail-closed unconfigured deployment
# ═══════════════════════════════════════════════════════════════════════════

class TestOpenModeNeverImplicit:

    def test_open_mode_off_by_default(self):
        """No TOSCANINI_ALLOW_OPEN_MODE in env -> OPEN_MODE must be False
        even though this test session has no real keys of its own."""
        assert os.environ.get("TOSCANINI_ALLOW_OPEN_MODE", "") != "1"
        import tos.security.auth as auth_mod
        assert auth_mod.OPEN_MODE is False

    def test_open_mode_requires_explicit_flag(self):
        """OPEN_MODE == (no keys) AND (TOSCANINI_ALLOW_OPEN_MODE == '1').
        Verified by re-importing the module source with each env combo."""
        import importlib, sys
        import tos.security.auth as auth_mod

        src_backup = {k: os.environ.get(k) for k in
                      ("TOSCANINI_API_KEYS", "TOSCANINI_API_KEY", "TOSCANINI_ALLOW_OPEN_MODE")}
        try:
            # no keys, no flag -> closed
            for k in src_backup:
                os.environ.pop(k, None)
            importlib.reload(auth_mod)
            assert auth_mod.OPEN_MODE is False

            # no keys + flag -> explicit open mode
            os.environ["TOSCANINI_ALLOW_OPEN_MODE"] = "1"
            importlib.reload(auth_mod)
            assert auth_mod.OPEN_MODE is True

            # keys + flag -> flag alone must NOT open a configured deployment
            os.environ["TOSCANINI_API_KEYS"] = "k1"
            importlib.reload(auth_mod)
            assert auth_mod.OPEN_MODE is False
        finally:
            for k, v in src_backup.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
            importlib.reload(auth_mod)

    def test_unconfigured_deployment_returns_503_on_missing_key(self):
        """No keys configured, no open flag: a request WITHOUT a key hits the
        fail-closed 503 path; a request WITH any key falls through to DB
        validation and gets 401."""
        import tos.security.auth as auth_mod
        probe = FastAPI()

        @probe.get("/probe")
        def probe_route(_auth=Depends(auth_mod.verify_api_key)):
            return {"ok": True}

        with patch.object(auth_mod, "ALLOWED_KEYS", set()), \
             patch.object(auth_mod, "OPEN_MODE", False):
            client = TestClient(probe, raise_server_exceptions=False)
            r_no_key = client.get("/probe")
            assert r_no_key.status_code == 503
            r_bad_key = client.get("/probe", headers={"X-API-Key": "not-a-real-key"})
            assert r_bad_key.status_code == 401


# ═══════════════════════════════════════════════════════════════════════════
# WP-2: full route coverage — every route is public by explicit decision
# (/health, /laws, /definitions) or requires a key (401/503 without one).
# ═══════════════════════════════════════════════════════════════════════════

VALID_KEY = "test-key-xyz"  # matches TOSCANINI_API_KEYS set at module import

# (method, path, no_key_status, with_key_statuses)
#   no_key_status     — exact expected status without X-API-Key
#   with_key_statuses — acceptable statuses WITH a valid key (data-dependent
#                       outcomes like 404/422 allowed; auth failures excluded)
ROUTE_MATRIX = [
    # ── public by explicit decision (reference data / monitoring only) ──
    ("GET",  "/health",                       200, {200}),
    ("GET",  "/laws",                         200, {200}),
    ("GET",  "/definitions",                  200, {200}),
    # ── protected: core audit ──
    ("POST", "/ingest",                       401, {422}),          # form validation only
    ("POST", "/v1/batch",                     401, {422}),          # zip required
    ("GET",  "/audit/{audit_id}",             401, {200}),
    ("POST", "/search",                       401, {422}),          # query form required
    # ── protected: NKG / governance ──
    ("GET",  "/nkg",                          401, {200}),
    # ── protected: comparisons ──
    ("GET",  "/comparison/{o}/{r}",           401, {200}),     # legacy tuple-return renders as 200
    ("GET",  "/comparisons/by-original/{id}", 401, {200}),
    ("POST", "/compare",                      401, {422}),          # form required
    # ── protected: refinement (B1 token flow) ──
    ("GET",  "/refinement/pdb/{audit_id}",    401, {404}),
    ("POST", "/refinement/token",             401, {422}),          # form required
    ("POST", "/refinement/callback",          422, {422}),          # token-gated, not key-gated
    # ── protected: refinement (B2 GPU flow) ──
    ("POST", "/refinement/submit",            401, {422}),          # form required
    ("GET",  "/refinement/status/{job_id}",   401, {200}),
    ("GET",  "/refinement/jobs/{audit_id}",   401, {200}),
    # ── protected: history & usage (key-scoped) ──
    ("GET",  "/history/{user_email:path}",    401, {200}),   # template requested literally
    ("GET",  "/usage",                        401, {200, 401}),     # static keys lack a DB record
    # ── removed in WP-2 (was unreachable; exposed IPs). The path now falls
    #    through to the /history/{user_email:path} catch-all and is key-scoped.
    ("GET",  "/history/ip/1.2.3.4",           401, {200}),
    # ── WP-3: HMAC-secured. Without a signature: 401 when the server secret is
    #    configured, 503 when it is not (both fail-closed; never 2xx).
    ("POST", "/refinement/internal-callback", {401, 503}, {401, 503}),
]

# Rows pending a stacked work package — reported xfail until then.
XFAIL_PENDING = {}


def _route_entries():
    """Yield (method, path) for every matrix row."""
    for entry in ROUTE_MATRIX:
        yield entry[0], entry[1]


@pytest.fixture(scope="module")
def api_client():
    """TestClient over the real app, keys configured (module env)."""
    import main as main_mod
    return TestClient(main_mod.app, raise_server_exceptions=False)


class TestRouteCoverage:

    @pytest.mark.parametrize("method,path,no_key,with_key", ROUTE_MATRIX,
                             ids=[f"{m} {p}".replace("/", "_") for m, p in _route_entries()])
    def test_route_auth(self, api_client, method, path, no_key, with_key):
        """Every route: exact status without a key, non-auth status with one."""
        if (method, path) in XFAIL_PENDING:
            pytest.xfail(XFAIL_PENDING[(method, path)])
        no_key = {no_key} if isinstance(no_key, int) else no_key
        fn = getattr(api_client, method.lower())
        r_no_key = fn(path)
        assert r_no_key.status_code in no_key, \
            f"{method} {path}: no-key status {r_no_key.status_code} not in {no_key}"

        r_key = fn(path, headers={"X-API-Key": VALID_KEY})
        assert r_key.status_code in with_key, \
            f"{method} {path}: keyed status {r_key.status_code} not in {with_key}"

    def test_all_app_routes_covered(self, api_client):
        """Guards against route drift: every registered API route (excluding
        FastAPI internals) must appear in ROUTE_MATRIX."""
        covered = {(m, p.split("{")[0]) for m, p in _route_entries()}
        app = api_client.app
        registered = set()
        for route in app.routes:
            path = getattr(route, "path", "")
            if path in ("/openapi.json", "/docs", "/docs/oauth2-redirect", "/redoc"):
                continue
            for m in getattr(route, "methods", set()) - {"HEAD", "OPTIONS"}:
                registered.add((m, path))
        for m, path in registered:
            assert any(m == cm and path.startswith(cp) for cm, cp in covered), \
                f"Route {m} {path} is missing from ROUTE_MATRIX"

    def test_invalid_key_rejected_on_protected_route(self, api_client):
        r = api_client.get("/nkg", headers={"X-API-Key": "totally-wrong-key"})
        assert r.status_code == 401
