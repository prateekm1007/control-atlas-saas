"""
WP-13 (SEC-015) regression tests: JWT signing secret must fail closed.

Contract under test (backend/tos/security/tokens.py):
- TOSCANINI_JWT_SECRET is REQUIRED. The historical hardcoded development
  fallback ("toscanini-phase-b1-secret-SETMEINENV") is REMOVED.
- Missing, empty, whitespace-only, or short (<32 chars) configuration raises
  TokenSecretConfigurationError from BOTH create_refinement_token and
  validate_refinement_token.
- A compliant secret supports the normal round-trip.
- Tokens signed under a different secret are rejected (ValueError from the
  InvalidTokenError wrapper), never accepted.
- Token format/expiry/purpose/single-use semantics are unchanged (see also
  test_b1_tokens.py and test_b6_regression.py, which run under the explicit
  test-only secret injected by conftest).

No test in this file depends on a live production secret.
"""
import pytest
import jwt as pyjwt

from tos.security.tokens import (
    TokenSecretConfigurationError,
    _get_jwt_secret,
    create_refinement_token,
    validate_refinement_token,
)

TEST_SECRET = "test-only-wp13-secret-0123456789abcdef0123456789abcdef-NOT-FOR-PRODUCTION"
ALT_SECRET = "alternative-wp13-secret-0123456789abcdef0123456789abcdef-DIFFERENT"


# ── fail-closed: missing configuration ────────────────────────────────────────

def test_missing_secret_rejects_signing(monkeypatch):
    monkeypatch.delenv("TOSCANINI_JWT_SECRET", raising=False)
    with pytest.raises(TokenSecretConfigurationError):
        create_refinement_token("AUDIT-WP13-1")


def test_missing_secret_rejects_validation(monkeypatch):
    monkeypatch.setenv("TOSCANINI_JWT_SECRET", TEST_SECRET)
    token = create_refinement_token("AUDIT-WP13-2")
    monkeypatch.delenv("TOSCANINI_JWT_SECRET", raising=False)
    with pytest.raises(TokenSecretConfigurationError):
        validate_refinement_token(token)


def test_missing_secret_via_helper(monkeypatch):
    monkeypatch.delenv("TOSCANINI_JWT_SECRET", raising=False)
    with pytest.raises(TokenSecretConfigurationError, match="not configured"):
        _get_jwt_secret()


# ── fail-closed: invalid / inadequate configuration ──────────────────────────

@pytest.mark.parametrize("bad_secret", ["", "   ", "short", "a" * 31])
def test_inadequate_secret_rejected(monkeypatch, bad_secret):
    monkeypatch.setenv("TOSCANINI_JWT_SECRET", bad_secret)
    with pytest.raises(TokenSecretConfigurationError):
        create_refinement_token("AUDIT-WP13-3")
    with pytest.raises(TokenSecretConfigurationError):
        _get_jwt_secret()


def test_minimum_length_boundary_accepts(monkeypatch):
    monkeypatch.setenv("TOSCANINI_JWT_SECRET", "x" * 32)
    assert len(_get_jwt_secret()) == 32


# ── valid configuration supports the normal round-trip ───────────────────────

def test_valid_secret_round_trip(monkeypatch):
    monkeypatch.setenv("TOSCANINI_JWT_SECRET", TEST_SECRET)
    token = create_refinement_token("AUDIT-WP13-4", "ops@example.com",
                                    extra_claims={"created_by_key_hash": "ab" * 32})
    payload = validate_refinement_token(token)
    assert payload["audit_id"] == "AUDIT-WP13-4"
    assert payload["purpose"] == "refinement_callback"
    assert payload["version"] == "1.0"
    assert payload["created_by_key_hash"] == "ab" * 32


# ── foreign-secret tokens are rejected ───────────────────────────────────────

def test_foreign_secret_token_rejected(monkeypatch):
    monkeypatch.setenv("TOSCANINI_JWT_SECRET", ALT_SECRET)
    forged = pyjwt.encode(
        {"audit_id": "AUDIT-WP13-5", "purpose": "refinement_callback",
         "version": "1.0", "exp": 9999999999},
        ALT_SECRET, algorithm="HS256",
    )
    monkeypatch.setenv("TOSCANINI_JWT_SECRET", TEST_SECRET)
    with pytest.raises(ValueError):
        validate_refinement_token(forged)


def test_publicly_known_dev_fallback_literal_removed():
    """The removed fallback literal must no longer exist in the module source.

    (Setting it EXPLICITLY as an env value is an operator choice; the defect
    was the implicit default, proven removed by the missing-secret tests and
    by this source-level assertion.)
    """
    import inspect
    import tos.security.tokens as tokens_module
    assert "SETMEINENV" not in inspect.getsource(tokens_module)


# ── preserved semantics smoke (full coverage lives in B.1/B.6 tests) ─────────

def test_single_use_semantics_preserved(monkeypatch, tmp_path):
    monkeypatch.setenv("TOSCANINI_JWT_SECRET", TEST_SECRET)
    monkeypatch.setenv("TOSCANINI_DATA_DIR", str(tmp_path))
    token = create_refinement_token("AUDIT-WP13-7")
    assert validate_refinement_token(token, consume=True)["audit_id"] == "AUDIT-WP13-7"
    with pytest.raises(ValueError, match="single-use"):
        validate_refinement_token(token)
