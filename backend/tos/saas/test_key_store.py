"""Quick smoke test for key_store."""
from key_store import create_key, get_key, check_quota, increment_usage, hash_key, revoke_key

print("Creating free tier key...")
key = create_key("free")
print(f"Key: {key[:20]}... (truncated)")

key_hash_val = hash_key(key)
print(f"Hash: {key_hash_val[:16]}...")

print("\nLookup key...")
api_key = get_key(key_hash_val)
print(f"Tier: {api_key.tier}, Active: {api_key.active}")

print("\nCheck quota...")
allowed, used, limit = check_quota(key_hash_val)
print(f"Allowed: {allowed}, Used: {used}/{limit}")

print("\nIncrement usage 3 times...")
for i in range(3):
    increment_usage(key_hash_val)
    allowed, used, limit = check_quota(key_hash_val)
    print(f"  After {i+1}: {used}/{limit}")

print("\nRevoke key...")
revoke_key(key_hash_val)
api_key = get_key(key_hash_val)
print(f"Active after revoke: {api_key.active}")

print("\n✅ All tests passed")
