# Item 15 evidence — runner hardening

**Base SHA:** 9de4fff · **Executed:** 2026-10-19T11:57:25Z

## Changes (backend/calibration/runner.py)
1. out_dir.exists() -> SystemExit 'refusing to overwrite' (evidence runs are
   append-only); mkdir no longer uses exist_ok.
2. Hardcoded 'origin.pdb' filename replaced by content-detected format
   (structure.cif / structure.pdb via detect_structure_format).
3. Provenance records the commit of the repo containing the IMPORTED tos
   engine package + working_tree_clean + tree_status. The CLI (main)
   refuses to start an evidence run on a dirty tree (exit 2) — the exact
   v24 failure mode (post-fix numbers under pre-fix SHA) can no longer
   occur silently.
4. run_status 'complete' only when measured == expected_target_count AND
   zero failures (e.g. 50 of 50).
5. Exit non-zero on any target failure (retained; verified by tests).

## Commands and results
```
$ cd backend && python -m pytest tests/test_runner_hardening.py tests/test_calibration_runner.py -q
11 passed   (6 new hardening + 5 existing runner tests, no regressions)

$ cd backend && python -m pytest tests -q
506 passed, 5 skipped, 1209 warnings in 8.78s
```

## Verdict changes
- None from hardening alone; AF- targets that now acquire (v6 cif) will be
  measured in v25 instead of NOT_MEASURED.
