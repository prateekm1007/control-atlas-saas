# Item 7 evidence — dead LAW-155 removed; one definition per LAW ID

**Base SHA:** a256980 · **Executed:** 2026-10-19T11:33:25Z

## Dead code verification (before deletion)
```
$ grep -rn "tier1_laws\|Tier1Laws" backend/ --include="*.py" | grep -v nested-copy
backend/tos/engine/tier1_laws.py:2:class Tier1Laws:      # zero importers

$ grep -rn "def _law_" backend/tos/ --include="*.py"
backend/tos/engine/tier1_laws.py:3:    def _law_155(self, atoms):   # the only _law_ def; no callers
```
The dead method implemented 'Steric Clash Prohibition' while the live LAW-155 is
Voxel Occupancy — a silent competing source of truth. Module deleted
(git rm); zero references outside the nested copy removed separately in item 17.

## New invariant tests (backend/tests/test_law_definition_uniqueness.py)
1. every LAW ID has exactly one measurement block (banner scan)
2. engine emits exactly the canon set (no orphans, no gaps)
3. no 'def _law_*' helper anywhere under backend/tos (pattern banned)
4. engine/tier1_laws.py stays deleted
5. registry enforced set == canon; retired == {190,210,220,230}

## Test-isolation defect found and fixed during this item
Full-suite runs exposed a cross-test/cross-run state leak:
tos_credits/credits.py resolved CREDITS_DIR at IMPORT time, and main.py
imports it lazily inside /refinement/submit — the first test to reach that
endpoint pinned the credits path for the whole session, defeating the
per-test TOSCANINI_DATA_DIR isolation and accumulating u@x.com's balance
across runs (proof: /tmp/toscanini-test-data/credits/u_at_x_com.json with
credits_used counting submits from earlier sessions; symptom:
test_ids.py::test_job_id_is_full_uuid_hex 402s only in full-suite runs).
Fix: _credits_dir() resolves the path per CALL from the environment;
module import no longer creates directories as a side effect.
This is the isolation the test fixtures already assumed — no test was
weakened, skipped, or deleted.

## Commands and results (three consecutive full runs, determinism proof)
```
$ cd backend && python -m pytest tests -q --tb=line   # run 1
472 passed, 5 skipped, 1209 warnings in 8.69s
$ cd backend && python -m pytest tests -q --tb=line   # run 2
472 passed, 5 skipped, 1209 warnings in 8.82s
$ cd backend && python -m pytest tests -q --tb=line   # run 3
472 passed, 5 skipped, 1209 warnings in 8.81s
```

## Files changed
- deleted: backend/tos/engine/tier1_laws.py
- added: backend/tests/test_law_definition_uniqueness.py (5 tests)
- fixed: backend/tos_credits/credits.py (per-call dir resolution)

## Verdict changes
- None: the deleted code had no callers and no runtime path.
