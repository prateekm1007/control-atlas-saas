"""
Toscanini Phase B1 — Token System Test Suite.

WP-7: rewritten from a script-with-sys.exit into proper pytest tests. The
original module-level harness (helper `test()` collected as a test +
module-level sys.exit(1)) broke collection.
"""
import secrets

import pytest

from tos.security.tokens import (
    create_refinement_token,
    validate_refinement_token,
    revoke_token,
    get_token_stats,
)


def test_token_generation():
    token = create_refinement_token("AUDIT001", "test@test.com")
    assert token is not None and len(token) > 50


def test_token_validation():
    token = create_refinement_token("AUDIT002", "test@test.com")
    payload = validate_refinement_token(token)
    assert payload["audit_id"] == "AUDIT002"
    assert payload["user_email"] == "test@test.com"
    assert payload["purpose"] == "refinement_callback"


def test_token_without_email():
    token = create_refinement_token("AUDIT003")
    payload = validate_refinement_token(token)
    assert payload["audit_id"] == "AUDIT003"
    assert payload["user_email"] is None


def test_invalid_token_rejected():
    with pytest.raises(ValueError):
        validate_refinement_token("invalid.token.string")


def test_tampered_token_rejected():
    token = create_refinement_token("AUDIT005")
    tampered = token[:-5] + "XXXXX"
    with pytest.raises(ValueError):
        validate_refinement_token(tampered)


def test_single_use_enforcement():
    token = create_refinement_token("AUDIT006_UNIQUE_" + secrets.token_hex(4))
    payload = validate_refinement_token(token, consume=True)
    assert payload is not None
    with pytest.raises(ValueError, match="already been used"):
        validate_refinement_token(token)


def test_token_stats():
    audit_id = "AUDIT007_STATS_" + secrets.token_hex(4)
    create_refinement_token(audit_id)
    stats = get_token_stats(audit_id)
    assert stats["audit_id"] == audit_id
    assert stats["callbacks_used"] == 0
    assert stats["callbacks_remaining"] == 3


def test_rate_limiting_after_max_callbacks():
    audit_id = "AUDIT008_RATE_" + secrets.token_hex(4)
    for _ in range(3):  # max callbacks per audit
        token = create_refinement_token(audit_id)
        validate_refinement_token(token, consume=True)
    token4 = create_refinement_token(audit_id)
    with pytest.raises(ValueError, match="Maximum callbacks"):
        validate_refinement_token(token4)


def test_token_uniqueness_nonce():
    t1 = create_refinement_token("AUDIT009")
    t2 = create_refinement_token("AUDIT009")
    assert t1 != t2, "Tokens should be unique even for same audit_id"


def test_token_revocation():
    audit_id = "AUDIT010_REVOKE_" + secrets.token_hex(4)
    token = create_refinement_token(audit_id)
    revoke_token(token, audit_id)
    with pytest.raises(ValueError, match="already been used"):
        validate_refinement_token(token)


def test_extra_claims_carried_and_protected():
    """WP-2: extra_claims reach the payload; base claims cannot be overwritten."""
    token = create_refinement_token(
        "AUDIT011_CLAIMS", extra_claims={"created_by_key_hash": "abc", "purpose": "hijack"}
    )
    payload = validate_refinement_token(token)
    assert payload["created_by_key_hash"] == "abc"
    assert payload["purpose"] == "refinement_callback"  # base claim preserved
