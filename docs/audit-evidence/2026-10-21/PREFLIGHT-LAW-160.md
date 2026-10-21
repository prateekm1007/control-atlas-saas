# Preflight - LAW-160 sequence-pair adjacency correction

- Recorded: 2026-10-21T11:50:02Z (UTC)
- Base: `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (freshly fetched and verified before the branch was cut)
- Branch: `fix/law-160-sequence-adjacency` (new feature branch; PR #37's governance branch is
  NOT reused)
- Task: correct the LAW-160 measurement defect only.

## Problem statement

`LAW-160` (chain integrity, threshold 4.5 A, `backend/tos/engine/tier1_measurements.py`) measured
the Cα-Cα distance of every ADJACENT PAIR in the list returned by
`coord_math.get_sequential_residues()`. That helper sorts and de-duplicates residues; it does not
establish sequence adjacency. Two residues whose sequence numbers are not adjacent (residues in
between are missing) were therefore measured as if they were consecutive, and the long distance
across the gap was reported as a chain break (a false VETO).

## Expected behavior after the change

1. Only genuinely consecutive residue pairs (the adjacency rules already used by the other Tier-1
   laws, `Tier1Measurements._is_sequential`, including insertion codes) are measured.
2. `sample` counts only eligible pairs that were actually measured.
3. A structure with no eligible, measured Cα pair reports an ABSENT metric: the adjudicator's
   existing absent-measurement representation (`status="FAIL"`, `observed=0`, `sample=0`,
   `backend/tos/engine/adjudicator.py`), never `PASS`.
4. A genuinely excessive distance between consecutive residues is still detected (VETO).

## Files expected to change

- `backend/tos/engine/tier1_measurements.py` (the LAW-160 block only)
- `backend/tests/test_law_160.py` (regressions added; the two original tests unchanged)
- `docs/evidence/reading-note/READING_NOTE_LATEST.md` (this task's read record)
- `docs/audit-evidence/2026-10-21/` (this preflight and the verification record)
- `docs/evidence/push-log/PUSH_LOG.md` (two-phase push records)
- `SHA256SUMS` (manifest rows for the in-scope evidence files that changed)

## Explicitly NOT changed

The 4.5 A threshold, the law ID, the operator, the modality classification, the law canon, the
modality matrix, the calibration cohort, the dataset, the golden artifacts,
`coord_math.get_sequential_residues()` (other measurements depend on it), and every other law's
behavior.

## Risks

- A structure whose chain has no eligible Cα pair changes from `PASS` (sample 0) to `FAIL`
  (sample 0). This is intended (an absent metric is not a normal result) and follows the
  adjudicator's existing fail-closed mapping (deterministic FAIL -> VETO). It is a VERDICT-AFFECTING
  change for such structures and is recorded as such.
- Mis-pairing the adjacency rules would either drop genuine pairs or keep false ones; the
  regressions pin both directions.

## Regression tests

`backend/tests/test_law_160.py`: sequence gap not measured; only eligible measured pairs counted;
insertion-code adjacency measured; insertion-code transition (plain residue -> next residue with an
insertion code) NOT adjacent; zero sample reported as FAIL, not PASS; zero sample with two
nonconsecutive residues; excessive distance between consecutive residues still vetoes. The two
pre-existing tests are unchanged.

## Rollback

Revert the single commit that changes `backend/tos/engine/tier1_measurements.py`; the branch is
feature-only and nothing is merged.
