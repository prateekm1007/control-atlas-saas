# Item 8 evidence — honest statuses completed

**Base SHA:** dfb9735 · **Executed:** 2026-10-19T11:37:53Z

## What changed
- Adjudicator malformed-row default (present-but-statusless) aligned from FAIL to NOT_MEASURED
  (absent-law default was already NOT_MEASURED from WP-SCI-1 — verified by test).
- Per-law inventory completed: evidence/item-8/LAW_INVENTORY.md (all 15 enforced laws:
  measurement, sample unit, empty-sample status, atom scope, pairing rule; contained + retired sets).
- LAW-110 counts observed peptide bonds (was already true on science head; now pinned by
  tests in test_residue_ordering.py: sample 3 and 2 on 4-residue fixtures).

## Test commands and results
```
$ cd backend && python -m pytest tests/test_honest_statuses_coverage.py -q
9 passed in 0.30s

$ cd backend && python -m pytest tests -q
481 passed, 5 skipped, 1209 warnings in 8.79s
```

## Verdict changes
- Malformed measurement rows no longer register as measured-FAIL; they block PASS via
  the indeterminate path (still fail-closed, now honestly labeled).
