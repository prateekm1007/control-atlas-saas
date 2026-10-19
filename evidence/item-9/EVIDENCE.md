# Item 9 evidence — coverage denominator from source residues

**Base SHA:** dfb9735 · **Executed:** 2026-10-19T11:37:53Z

## What changed
- processor.Structure gains source_residue_count: distinct polymer residue keys seen in
  the SOURCE records (ATOM + retained HETATM polymer), counted before per-atom parse
  filtering — unparseable coordinate lines can no longer shrink the denominator.
- run_full_audit: coverage = len(core) / source_residue_count (predicted structures);
  LAW-105 row gains sample_source_residues alongside the parsed-residue sample.

## Demonstrations (tests/test_honest_statuses_coverage.py)
- 6 source residues, 1 with garbage coordinates -> coverage 83.33% (5/6), not 100%
- MSE HETATM residue counts in the denominator
- 2/6 residues below pLDDT 70 -> 66.67% -> LAW-105 FAIL

## Test commands and results
```
$ cd backend && python -m pytest tests/test_honest_statuses_coverage.py -q
9 passed in 0.30s
$ cd backend && python -m pytest tests -q
481 passed, 5 skipped, 1209 warnings in 8.79s
```

## Verdict changes
- Predicted structures with dropped/unparseable residues now score LOWER, honest coverage.
