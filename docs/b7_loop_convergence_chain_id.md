# LOOP CONVERGENCE SPEC — ADDENDUM: chain_id Decision

Recorded: 2026-03-01
Status: Design decision. No code. Freeze respected.

## chain_id

  Generated at: create_refinement_token()
  Storage:      JWT payload field "chain_id"
  Format:       UUID4 (secrets.token_hex(16))
  Derivation:   NOT content-addressed. Chains are user sessions,
                not structural hashes. Same structure submitted twice
                in different sessions = two different chains.
  Persistence:  Passed through validate_refinement_token() into
                comparison metadata on callback. Survives across
                callbacks for the same token-generation event.
  Scope:        One chain = one create_refinement_token() call
                and all callbacks that consume tokens from that chain.

## Implementation Note (for B.7, post-freeze)

  1. Add chain_id = secrets.token_hex(16) to create_refinement_token()
  2. Include in JWT payload alongside audit_id, jti
  3. Extract in validate_refinement_token(), return in payload dict
  4. Pass into store_comparison() as chain_id field
  5. No audit schema changes. Comparison metadata only.
