# C5 diagnostic - exact base vs OD-scope candidate

Recorded: 2026-10-22T13:37:28Z
Base engine: `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
Candidate engine: `3029e21196dedd7cecd3927598639c9bd1a22e6c`
Branch: `review/od-scope-v25`
Status: diagnostic evidence for owner review. Not a threshold validation.

## Scope and provenance

- Engine checkouts were materialized byte-exactly from git objects with `git archive` and evaluated in isolated subprocesses.
- Coordinate bytes were taken from the pinned v24 manifest and verified against the pinned byte counts and SHA-256 values.
- The v24 coordinate manifest pin is `9722c3e45d316ff58ed326df074eec138b81ceb520191a492088cc0ddef63762`; the run manifest pin is `f602351f4844531d409cfc256939508103e3c2d4f57fe21d7f77c34494aa4036`.
- This archive uses base vs candidate only. No prior `c78` comparison is used as evidence.

## Acquisition and analysis scope

- v24 manifest: 50 targets; 40 acquired and digest-verified; 10 historical AlphaFold `model_v4` targets unavailable.
- v24 primary set: 38 targets. `1I0T` and `1DDL` are excluded from primary counts and retained separately per OD-18.
- Supplementary current-AlphaFold set: 9 acquired; `AF-P0DTC2-F1` is excluded under OD-11; `AF-P62988-F1` is unavailable; 8 targets are analyzed. Acquisition records are in `supplementary-acquisition.json`.
- Total analyzed structures: 46 (38 v24-pinned primary + 8 supplementary).

## Verdict comparison

Primary v24 base: {'PASS': 34, 'VETO': 4}.
Primary v24 candidate: {'PASS': 33, 'VETO': 4, 'INDETERMINATE': 1}.
Primary verdict changes: 1 of 38.
Supplementary base: {'PASS': 5, 'INDETERMINATE': 3}.
Supplementary candidate: {'PASS': 5, 'INDETERMINATE': 3}.
Supplementary verdict changes: 0 of 8.

| Target | Base verdict | Candidate verdict |
|---|---|---|
| 5T4O | PASS | INDETERMINATE |

The sole primary verdict change is `5T4O`: `PASS` to `INDETERMINATE`. Its candidate LAW-150 row is `NOT_MEASURED` because no eligible residue had measurable chi1 coordinates; the former zero-sample `PASS` is removed.

## Law status transitions (analyzed sets)

Total changed law rows: 102 (79 primary, 23 supplementary).

| Law | Base status | Candidate status | Rows |
|---|---|---|---|
| LAW-100 | PASS | FAIL (Advisory) | 2 |
| LAW-150 | PASS | NOT_MEASURED | 1 |
| LAW-155 | PASS | RETIRED | 46 |
| LAW-182 | VETO | FAIL (Advisory) | 16 |
| LAW-195 | PASS | NOT_APPLICABLE | 33 |
| LAW-200 | VETO | FAIL (Advisory) | 4 |

These transitions are the intended truthfulness and policy effects: retirement of LAW-155, composition-driven `NOT_APPLICABLE` for LAW-195, advisory display for LAW-182/LAW-200, the canonical LAW-100 percentage path, and the removed zero-sample LAW-150 pass. Full per-target observations and samples are in `law-status-deltas-3029e21196de.json`.

## Separate audit set (OD-18)

`1DDL` and `1I0T` remain VETO in both engines and are not included in the primary counts. Their law rows are retained in `v24-pinned-separate-audit-*.json` and in the JSON delta file.

## Interpretation constraints

- Diagnostic only. No validation of v24 calibration, thresholds, cohorts, splits, or the historical verdict-change table (OD-9, OD-14, OD-16).
- No coordinate bytes were substituted or fabricated. Unavailable targets are disclosed, not replaced.
- No v25 calibration, held-out evaluation, threshold/cohort change, merge, or deployment was performed.
- The archive commit that carries this file is not the engine SHA named above; the engine SHA is the exact code candidate materialized for this comparison.
