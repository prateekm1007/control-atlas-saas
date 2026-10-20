# Preflight — Round-12 Direction 2: OD-07 audit identity (2026-10-20)

EC §2.5 / §14 artifact. Companion to `CONSTITUTION_READING_NOTE_2026-10-20.md`.

## Baseline

- Repository: `prateekm1007/control-atlas-saas`
- Branch: `feature/od07-audit-id` (created from `science/p0-correctness-2026-10-19` @ `ebc1fd79f20b9a6b6a262a70ce6942ba069d19e9`, per audit Direction 2)
- Working tree: clean (fresh worktree at branch creation)
- Confirmation required by the direction: `backend/main.py:182` ALREADY computes
  the deterministic audit ID — `hashlib.sha256(coord_hash.encode()).hexdigest()`
  (full 64-hex SHA-256 of the canonical coordinate hash, WP-5). It is kept
  unchanged; it is NOT replaced with a UUID. Covered by existing
  `tests/test_ids.py::TestAuditIds` and by the new
  `test_od07_audit_identity.py::test_audit_id_still_sha256_of_coord_hash`.
- `structure_digest`: NO field of that name exists anywhere at this base
  (`git grep -n "structure_digest" ebc1fd7` → no hits, whole tree). The OD-07
  clause "structure_digest is informational only, never an identifier" is
  therefore satisfied vacuously at this base; the nearest analog,
  `provenance.hash` (the canonical coordinate hash), is documented as
  provenance and is not used as a lookup identifier.

## Problem statement (audit finding)

OD-07 ruling not implemented on `main`: audit ID was `md5[:8]` (8 chars,
uppercase, collision-prone) at `backend/main.py:171`; `uuid4()[:8]` identifiers
were generated in the ingestion layer and could reach stored records; no
legacy lookup for 8-character references; no key scoping on audit-object
lookups.

## Change map (inputs → outputs)

- `backend/tos/storage/audit_store.py` — ADD `compute_legacy_audit_id()`
  (deterministic pre-WP-5 scheme `md5(coord_hash)[:8].upper()`), ADD
  `resolve_audit_ref()` (canonical 64-hex OR legacy 8-char resolution;
  fail-closed key scoping; unattributed records invisible to any key).
  No change to existing `store_audit_result` / `get_audit_result`.
- `backend/main.py` — `_run_physics_sync(..., key_hash="")`: payload gains
  `governance.audit_identity` (triple: audit_id + canon_hash + matrix_hash),
  `governance.legacy_audit_id`, and top-level `created_by_key_hash`.
  `/ingest` and `/v1/batch` attribute records to the calling key.
  `/audit/{ref}` resolves canonical or legacy refs, key-scoped (fail-closed).
  `/refinement/pdb/{ref}` resolves legacy refs to the canonical filename
  (auth model unchanged — documented in-code). `/refinement/callback` and
  `/refinement/internal-callback` inherit key attribution from token claims /
  job record respectively.
- `backend/tos/schemas/response_v1.py` — ADD `AuditIdentity` model; `Governance`
  gains optional `audit_identity` + `legacy_audit_id` (additive, backward
  compatible).
- `backend/tos/ingestion/processor.py` — REMOVE `uuid4()[:8]` generation and
  the `uuid` import; ingestion-level `Structure.audit_id` is now `""` (the
  audit identity lives only at the API layer per OD-07). The fixed
  empty-sample sentinel `AUD-empty` is NOT random and is left as-is.
- `backend/tests/test_od07_audit_identity.py` — NEW: 13 regression tests
  (identity triple; determinism; SHA-256-not-UUID; legacy scheme; no ingestion
  identifier; key attribution; canonical lookup; cross-key denial; anonymous
  fail-closed; legacy owner/deny; pre-WP-5 fail-closed; invalid refs; PDB
  legacy resolution; batch attribution).

## Risks and mitigations

- Response contract growth → additive optional fields only; existing consumers
  unaffected (verified by full suite).
- Legacy scan cost → Phase-1 JSON store; bounded directory scan; corrupted
  records skipped, never fatal.
- Key-scoping regressions for internal flows → internal callers use
  `key_hash=None` (no filtering); user-facing `/audit` is scoped; both tested.
- Test isolation → `test_api_key_guard.py` reloads `tos.security.auth` and
  resets `ALLOWED_KEYS`; the new test module re-adds its second key per test
  (fixture), documented in-code.

## Evidence (exact commands, this branch)

- Targeted: `cd backend && TOSCANINI_DATA_DIR=<tmp> python -m pytest
  tests/test_od07_audit_identity.py tests/test_ids.py -q`
  → **24 passed** (13 new + 11 WP-5).
- Full suite: `cd backend && TOSCANINI_DATA_DIR=<tmp> python -m pytest tests/
  -q -p no:warnings` → **514 passed, 9 failed, 5 skipped**.
- Clean-base comparison (`ebc1fd7`, fresh worktree, same interpreter):
  **497 passed, 9 failed, 5 skipped** — the SAME 9 failures
  (4× `test_afdb_acquisition` `ModuleNotFoundError: gemmi` — gemmi is not
  installable in this sandbox; 4× `test_b6_regression::TestUsageWidget`;
  1× `test_runner_hardening::test_measure_target_uses_content_detected_format`).
  **Zero new failures attributable to this branch.**

## Rollback

Single-commit revert on this branch (no history rewrite); no schema
migration, no data backfill required (legacy records remain readable by
their filenames).

## Final checklist

- [x] Base SHA recorded; tree clean at branch creation
- [x] Contracts/tests/callers read before edit
- [x] Regression + negative-path tests added (13)
- [x] Full suite delta vs clean base: zero new failures
- [x] No threshold/law/verdict change; defect-traceable only (OD-14)
- [x] No secrets introduced; no push to `main`; no merge performed
- [ ] CI on the pushed branch SHA (runs after push — result reported, not assumed)
