# Item 10 evidence — pLDDT scope of LAW-130/155/200

**Base SHA:** d5fa87a · **Executed:** 2026-10-19T11:41:14Z

## Verification (code)
LAW-130/155/200 consume structure.atoms (ALL atoms); LAW-100..160 consume core
(pLDDT >= 70 for predicted). Cited lines in evidence/item-10/DECISION_REQUIRED.md.

## Verification (demonstrated divergence)
scope_demo.py: compact pLDDT-90 core + extended pLDDT-40 tail ->
LAW-200: 477.01 A^3/atom VETO (all atoms) vs 68.09 PASS (core only) — verdict flips.
Output: evidence/item-10/scope_demo_output.txt

## Disposition
Verdict-changing scientific-semantics choice -> OWNER GATE. No law code changed.
Current scope pinned by backend/tests/test_plddt_scope_pin.py (3 tests) so it
cannot drift silently; the pin tests are the ones to update in the same PR if
the owner selects core-filtering. See DECISION_REQUIRED.md for both positions.

## Test commands and results
```
$ cd backend && python -m pytest tests/test_plddt_scope_pin.py -q
3 passed in 0.31s
$ cd backend && python -m pytest tests -q
484 passed, 5 skipped, 1209 warnings in 8.78s
```
