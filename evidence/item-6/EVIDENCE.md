# Item 6 evidence — altloc selection (one conformer per atom)

**Base SHA:** 490ab36 · **Executed:** 2026-10-19T11:25:34Z

## Rule (applied everywhere)
Highest occupancy wins; ties keep the FIRST record in file order.

Applied in four places:
1. IngestionProcessor (parse-time dedup — LAW-130/155/200 consume structure.atoms directly)
2. coord_math.get_atom_pos (was: first match)
3. coord_math.get_residue_atoms (was: last match — disagreed with get_atom_pos)
4. Tier1Measurements._extract (was: last record; _conf now averages RETAINED conformers only)

## Test on a structure with altlocs (real bytes)
3ERD: 66 altloc atom records; ASP A 321 backbone in conformers A/B at 0.50/0.50 —
a genuine tie-break case. Pinned fixture sha256:
7a192385c93b90670a52e8fc478b441152a49811b93e6ca202c9a487cb4e7418

## Commands and results
```
$ cd backend && python -m pytest tests/test_altloc_selection.py -q
8 passed in 0.70s

$ cd backend && python -m pytest tests -q
467 passed, 5 skipped, 1209 warnings in 8.66s
```

## Verdict changes
- LAW-130/155/200 no longer double-count altloc atoms (denominators shrink on
  altloc-bearing structures); atom lookups are consistent across all paths.
- No law threshold or meaning changed.
