# OD-scope iteration 4 - constitution-first reading note

Recorded: 2026-10-22T17:55:25Z (hash/decision scan); 2026-10-22T17:56:40Z (remote verification)
Branch: `review/od-scope-v25`
Base: `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (`origin/main`)
Engine candidate: `3029e21196dedd7cecd3927598639c9bd1a22e6c`
Evidence/archive candidate: `e71cec36d8189125b74568343a80fff1eb9637cd`
Working tree at note start: clean, except no unreviewed modifications
PAT: owner rotates on the owner's schedule; no PAT material is read, stored, or asserted here

## Governing-file hashes and blobs

| File | SHA-256 | Git blob at `e71cec3` |
|---|---|---|
| `AGENTS.md` | `5e320939a11135a0d8aed5caa231df6aeaf633d724fcb4dba97b5909e51420b4` | `85767237399ffb30f3e45894d9dfbe6f4240584f` |
| `constitution.md` | `cfdd5ba3babda31513b55894a49597fbfb18c2e870535a068c8b0ea423c53ce8` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` |
| `ENGINEERING_CONSTITUTION.md` | `e0f256f78f2c1256228c4286311cba15800dcd5a6c8005583508ee66d27983f5` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` |
| `docs/coder/BRIEF.md` | `8f07edfbba9d3e57de623e242125e942369bc5d60df5805343b26b9048c52e07` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` |
| `docs/decisions/OWNER_DECISIONS.md` | `59ccd5e84176a8a610c8c234cbb95fec4734478cc1f35668599ad2e3aef3172a` | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` |
| `docs/P0_STATUS.md` | `3f129dd18b9e5cf8888b763158cda42de270a8429a48ed7cc9a3aa675e810e0f` | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` |

All six files were re-read in full. No governing file is edited by this iteration.

## Live remote state

Verified by `git ls-remote --heads origin` at `2026-10-22T17:56:40Z` (full capture in `remote-heads-pre-fix.txt`):

- `main` = `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
- `review/od-scope-v25` = `e71cec36d8189125b74568343a80fff1eb9637cd`
- `review/od-scope-v25-2026-10-22` = `e2cd17d0f9c5231c489c42b8dd90d74055d1f648` (do not modify)

The live remote ref for `review/od-scope-v25` is exactly the archive candidate B.
Owner-reported CI state at B: Stability Gate `39825159703` success (459 passed, 5 skipped, 1,262 warnings, Python 3.12), Governance Acceptance `39825159704` failure (package manifest stale `PUSH_LOG.md` digest; push-log remote verification fails closed on pruned historical refs).

## Decision map applied to this iteration

- **OD-16 remains active and central**: no merge to `main`, no deployment, no v25 calibration run, no held-out evaluation, no threshold/cohort change, no golden-artifact regeneration. Only normal feature-branch pushes are authorized.
- **OD-02/OD-03/OD-05**: coverage-first fail-closed verdicts; pLDDT advisory-only; retired laws (190/210/220/230) stay out of the enforced active set.
- **OD-09/OD-10**: no calibration/threshold tuning, no editable cohort construction, no threshold derived from the contaminated historical cohort.
- **OD-11**: `AF-P0DTC2-F1` (identity `AF-0000000365840314`) stays out of the authoritative cohort; no silent ID substitution.
- **OD-13**: pin exact model/version, entity, HTTPS source, format/size, SHA-256, and run-manifest provenance; no `latestVersion` movement.
- **OD-14**: the historical 18-degree LAW-120 threshold and the historical verdict-change table remain diagnostic-only and unapproved.
- **OD-17**: keep `backend/calibration/dataset_v23.py` and all historical v23 results byte-for-byte; corrected metadata must be a new version with row-level sources.
- **OD-18**: `1I0T` (DNA-only) and `1DDL` (protein+RNA) stay out of the protein-only core; they remain a separately classified audit set.
- **EC-17/EC-22/EC-25** are active policy (live automated enforcement not yet on `main`); **EC-16/18/19/20/21/23/24/26** remain proposed/inert until exact-candidate checks and activation are recorded. Do not weaken or bypass any gate.
- Standing rule from `AGENTS.md` / `P0_STATUS.md`: no direct push to `main`, no force-push, no secrets, "evidence or it did not happen", decision IDs are never reused.

## What this iteration may do

- Normal follow-up commits on `review/od-scope-v25` only (no force, no rebase of published history).
- Focused engine/geometry/confidence/evaluator/runner corrections and their tests.
- Reconciliation of the push-log and package manifest with truthful, append-only records.
- Reproducible C5 regeneration on the next exact candidate, subsets kept separate per OD-18/OD-11.
- Stop for review; no merge, deploy, calibration, held-out run, threshold/cohort change, or golden-artifact regeneration.

## What this iteration may not do

- Edit `AGENTS.md`, `constitution.md`, `ENGINEERING_CONSTITUTION.md`, `docs/coder/BRIEF.md`, `docs/decisions/OWNER_DECISIONS.md`, or `docs/P0_STATUS.md`.
- Invent the absent evaluation protocol `backend/calibration/protocol/evaluation_protocol_v1.md`; it stays diagnostic-only until recovered and reviewed.
- Recreate pruned branches, fabricate intents, rewrite historical push-log entries, or weaken the two-phase schema.
- Rotate or read the PAT; the owner handles rotation independently.
