# N8 altloc-default verdict diff — pinned v23 cohort (2026-10-20)

Round-12 Direction 3 evidence. Diagnostic runs of the SAME pinned
cohort (dataset v23, 50 targets) before and after the N8
altloc-default flip. **No threshold change** (OD-14: change traces
to the altloc fixed defect; recorded diff only).

- BEFORE: `08a9403a7dd3` (feature/od07-audit-id tip, per-atom behavior) @ 2026-10-20T05:09:30.748291+00:00
- AFTER : `eb7e56523f24` (feature/n8-altlocs-per-residue tip, per-residue default) @ 2026-10-20T05:11:27.616091+00:00
- Input identity: 48 targets with matching SHA-256 coordinate digests, 0 mismatches (verified between runs).

## Aggregate verdicts

| Verdict | before | after |
|---|---|---|
| INDETERMINATE | 1 | 1 |
| PASS | 34 | 34 |
| VETO | 5 | 5 |
| not:MEASUREMENT_FAILED | 8 | 8 |
| not:NOT_MEASURED | 2 | 2 |

- expected_targets: 50 / acquired: 48 / measured: 40 (identical in both runs)
- acquisition_failures: 10 (recorded, identical in both runs; dominant cause: gemmi absent in this sandbox for mmCIF measurement)

## Verdict changes (per target)

**None.** No target changed its verdict across the altloc-default flip on this cohort.

## Law-status changes (per target, measured rows only)

**None.** No law row changed status across the flip.

## Interpretation

- The altloc policy only changes which alternate conformers are
  retained for residues carrying multiple conformers; every target
  in this cohort produced identical measurements under both
  policies here, so the recorded diff is empty. An empty diff is a
  RECORDED RESULT of these two runs over these pinned inputs — it
  is not a claim that the policies are equivalent in general.
- The per-residue selection logic itself is covered by
  backend/tests/test_altloc_per_residue.py (11 passing in this
  sandbox; the 3 CIF-path tests are gemmi-limited, pre-existing).
- JSON form: n8_verdict_diff_2026-10-20.json (same directory).
