"""
Toscanini Phase B1 — Token System Test Suite
Automated validation of all token security paths.
"""
import sys
sys.path.insert(0, '/app')

from tos.security.tokens import (
    create_refinement_token,
    validate_refinement_token,
    revoke_token,
    get_token_stats
)

passed = 0
failed = 0

def test(name, fn):
    global passed, failed
    try:
        fn()
        print(f"  ✓ {name}")
        passed += 1
    except Exception as e:
        print(f"  ✗ {name}: {e}")
        failed += 1

print("\n=== TOKEN SYSTEM TEST SUITE ===\n")

# TEST 1: Token generation
def t1():
    token = create_refinement_token("AUDIT001", "test@test.com")
    assert token is not None and len(token) > 50
test("Token generation", t1)

# TEST 2: Token validation
def t2():
    token = create_refinement_token("AUDIT002", "test@test.com")
    payload = validate_refinement_token(token)
    assert payload["audit_id"] == "AUDIT002"
    assert payload["user_email"] == "test@test.com"
    assert payload["purpose"] == "refinement_callback"
test("Token validation", t2)

# TEST 3: Token without email
def t3():
    token = create_refinement_token("AUDIT003")
    payload = validate_refinement_token(token)
    assert payload["audit_id"] == "AUDIT003"
    assert payload["user_email"] is None
test("Token without email", t3)

# TEST 4: Invalid token rejected
def t4():
    try:
        validate_refinement_token("invalid.token.string")
        assert False, "Should have rejected invalid token"
    except ValueError:
        pass
test("Invalid token rejected", t4)

# TEST 5: Tampered token rejected
def t5():
    token = create_refinement_token("AUDIT005")
    tampered = token[:-5] + "XXXXX"
    try:
        validate_refinement_token(tampered)
        assert False, "Should have rejected tampered token"
    except ValueError:
        pass
test("Tampered token rejected", t5)

# TEST 6: Single-use enforcement
def t6():
    token = create_refinement_token("AUDIT006_UNIQUE_" + __import__("secrets").token_hex(4))
    payload = validate_refinement_token(token, consume=True)
    assert payload is not None
    try:
        validate_refinement_token(token)
        assert False, "Should have rejected used token"
    except ValueError as e:
        assert "already been used" in str(e)
test("Single-use enforcement", t6)

# TEST 7: Token stats
def t7():
    audit_id = "AUDIT007_STATS_" + __import__("secrets").token_hex(4)
    token = create_refinement_token(audit_id)
    stats = get_token_stats(audit_id)
    assert stats["audit_id"] == audit_id
    assert stats["callbacks_used"] == 0
    assert stats["callbacks_remaining"] == 3
test("Token stats", t7)

# TEST 8: Rate limit after max callbacks
def t8():
    import secrets as s
    audit_id = "AUDIT008_RATE_" + s.token_hex(4)
    # Use 3 tokens (max)
    for i in range(3):
        token = create_refinement_token(audit_id)
        validate_refinement_token(token, consume=True)
    # 4th should be rate limited
    token4 = create_refinement_token(audit_id)
    try:
        validate_refinement_token(token4)
        assert False, "Should have rate limited"
    except ValueError as e:
        assert "Maximum callbacks" in str(e)
test("Rate limiting after max callbacks", t8)

# TEST 9: Nonce uniqueness (two tokens for same audit are different)
def t9():
    audit_id = "AUDIT009"
    t1 = create_refinement_token(audit_id)
    t2 = create_refinement_token(audit_id)
    assert t1 != t2, "Tokens should be unique even for same audit_id"
test("Token uniqueness (nonce)", t9)

# TEST 10: Revocation
def t10():
    import secrets as s
    audit_id = "AUDIT010_REVOKE_" + s.token_hex(4)
    token = create_refinement_token(audit_id)
    revoke_token(token, audit_id)
    try:
        validate_refinement_token(token)
        assert False, "Should have rejected revoked token"
    except ValueError as e:
        assert "already been used" in str(e)
test("Token revocation", t10)

print(f"\n=== RESULTS: {passed}/{passed+failed} passed ===\n")
if failed > 0:
    sys.exit(1)
