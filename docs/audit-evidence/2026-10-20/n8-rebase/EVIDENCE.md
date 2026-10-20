# N8 rebase onto the OD-07 fix — evidence (Round 15, direction 2)

Remote branch tip before: `b5afd6ca01a345432852d5865a79999e464636b4`
(descended from the NON-COMPLIANT `08a9403`; the Round-15 audit
correctly rejected that stack — merging it would have reverted the
OD-07 fix).

Rebase method: merge-based, append-only (a literal rebase + push would
require force, which is forbidden). The merge commit's parents are
`b5afd6c` (per_residue default + verdict-diff evidence) and `c02f774`
(OD-07 ruling implementation). The change-sets were disjoint, so the
merged tree is exactly their union; verified field-by-field:

- `processor.py`: `ALTLOC_POLICY = ALTLOC_POLICY_PER_RESIDUE` (R-4
  default preserved); `per_atom` opt-in and its pinned tests unchanged.
- `main.py`: `governance.audit_id = compute_audit_id(coord_hash,
  LAW_CANON_HASH, compute_matrix_hash())` — the canonical triple in ONE
  SHA-256 (OD-07).
- `tos/governance/audit_identity.py` present (single derivation point).
- The superseded `test_audit_id_still_sha256_of_coord_hash` is ABSENT
  from both test files; the ruling tests are present
  (`test_audit_id_is_not_the_bare_coord_hash`,
  `test_same_coords_different_canon_different_id`, matrix variant,
  integration variant, read-only legacy lookup).
- `scripts/migrate_audit_ids_od07.py` (explicit migration) carried.

## Commands run and output (Python 3.12.15; gemmi 0.7.5)

1. `cd backend && python3 -m pytest tests -q --tb=line`
   -> 4 failed, 537 passed, 5 skipped in 21.38s
   (the 4 failures are the environmental TestUsageWidget cases —
   streamlit absent locally; zero new failures vs `c02f774`'s run)
2. Hash locks: canon `6a9cd4b4349b81de`, matrix `ed1731ace09a9a3b…`,
   law-semantic vs `governance/law_semantic.sha256` — ALL MATCH.
3. Pinned-cohort verdict diff re-run on the merged tree (script mirrors
   `_run_physics_sync`; 40/40 experimental targets fetched from
   canonical RCSB bytes; both policies):
   -> verdict_changes: 0; deterministic-score changes: 0; failing-law
      changes: 0; fetch errors: 0; 13/40 targets carry altloc residues;
      the same 2 atom-set differences as before (1EJG, 3NIR — chimera
      prevention, both PASS/83 either way).
   Full JSON: `verdict_diff_after_od07_merge.json` (this directory).
   The verdict-diff evidence from `b5afd6c` therefore remains valid:
   the identity fix changes no verdict logic.

## Verdict impact

NONE. Altloc policy (unchanged from the accepted default) + identity
semantics (per the accepted OD-07 fix); hash locks bit-identical.
