# OD-07 compliance fix — evidence (Round 14, direction 2)

Base: `feature/od07-audit-id` @ `08a9403a7dd376da08831d14f1fa4fa26b15430e`
(the pushed commit — preserved, not rewritten; this fix is a new commit
on top, per direction 1/2).

## Deviation being fixed (auditor finding, verified)

At `08a9403`, `governance.audit_id` = `sha256(coord_hash)` (the bare
WP-5 form) while the canon and matrix hashes sat beside it in a separate
`governance.audit_identity` object. Same structure therefore kept the
same ID under different canon/matrix versions — the opposite of the
ruling. The superseded assertion
(`test_audit_id_still_sha256_of_coord_hash`) asserted that deviation.

## The fix

- `backend/tos/governance/audit_identity.py` (new, single derivation
  point): `compute_audit_id(structure_digest, canon_hash, matrix_hash)`
  = SHA-256 over the fixed-order, domain-tagged canonical serialization
  `"od07-audit-id-v1|{digest}|{canon}|{matrix}"`; `audit_id_from_record`;
  `is_full_audit_id`.
- `backend/main.py`: `governance.audit_id` = `compute_audit_id(coord_hash,
  LAW_CANON_HASH, compute_matrix_hash())`. The separate `audit_identity`
  object is REMOVED (the binding values remain in
  `governance_fingerprint`). `legacy_audit_id` (deterministic pre-WP-5
  md5[:8].upper ref) kept for read-only legacy resolution.
- `backend/tos/schemas/response_v1.py`: the `AuditIdentity` model and
  `audit_identity` field removed with their payload source; `legacy_audit_id`
  kept (Optional).
- `backend/tos/storage/audit_store.py`: UNCHANGED — the pushed
  implementation's `resolve_audit_ref` is already read-only (no writes on
  GET; verified by reading every write call in the module and by the new
  `test_legacy_lookup_writes_nothing`), and its key scoping
  (`created_by_key_hash`, fail-closed for unattributed records) is
  preserved as-is.
- `scripts/migrate_audit_ids_od07.py` (new): the EXPLICIT migration
  script for stored-record rewrites (ruling: "A separate migration
  script handles any stored-record rewrites"). Dry-run by default;
  `--apply` re-keys WP-5-era bare-hash records and metadata-bearing
  pre-WP-5 records to the canonical ID (inner `governance.audit_id`
  updated, `legacy_audit_id` backfilled, `created_by_key_hash`
  preserved); metadata-less records are reported and never touched;
  old files retained unless `--remove-legacy`. Smoke-tested end-to-end
  (see below).
- Tests: the superseded assertions rewritten to the ruling; the ruling
  tests added.

## Ruling tests (direction 2, verbatim requirement)

- `test_same_coords_different_canon_different_id` — same coordinates,
  different canon hash, DIFFERENT ID (unit, on the shared derivation).
- `test_same_coords_different_matrix_different_id` — matrix component.
- `test_same_input_different_canon_changes_response_id` — integration:
  the same uploaded file, monkeypatched canon hash, different response
  `audit_id`.
- `test_audit_id_is_not_the_bare_coord_hash` — the superseded bare form
  is rejected.
- `test_audit_id_binds_coord_canon_and_matrix` — response ID equals
  `compute_audit_id(provenance.hash, fingerprint.canon_hash,
  fingerprint.matrix_hash)`; `audit_identity` object absent.
- `test_legacy_lookup_writes_nothing` — the audit-dir file set is
  byte-identical before/after a legacy GET (read-only reaffirmation).
- `test_identical_input_shares_identity` — WP-5 dedup kept.

## Commands run and output (Python 3.12.15; gemmi 0.7.5)

1. `cd backend && python3 -m pytest tests/test_od07_audit_identity.py tests/test_ids.py -q`
   -> 28 passed
2. `cd backend && python3 -m pytest tests -q --tb=line`
   -> 4 failed, 523 passed, 5 skipped (the 4 failures are the
   environmental TestUsageWidget cases, present at the science base;
   zero new failures)
3. Hash locks: canon `6a9cd4b4349b81de`, matrix `ed1731ace09a9a3b…`,
   law-semantic vs `governance/law_semantic.sha256` — ALL MATCH.
4. Migration smoke test (temp store): WP-5-era bare-hash record re-keyed
   to the canonical ID with inner ID updated, legacy ref backfilled
   (`014842D4` = md5("a"*64)[:8].upper()), owner hash preserved,
   metadata-less record skipped and intact; dry run wrote nothing.

## Verdict impact

NONE. Identity semantics only: no law, threshold, operator, scope, or
adjudication change; hash locks bit-identical. Audit IDs CHANGE VALUE
by design (the ruling binds law semantics into the identity); identical
input under identical semantics still shares one ID.
