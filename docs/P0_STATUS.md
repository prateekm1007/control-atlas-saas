# P0_STATUS — Implementation Status, Blockers, Test Evidence, Candidate SHAs

**Status date:** 2026-10-20 (Round 14)
**Authority:** OD-01 target layout; maintained per the Engineering
Constitution (EC §6, §12, §14). Decision status lives in
`docs/decisions/OWNER_DECISIONS.md` (canonical copy on the
`docs/owner-decisions` branch); this file tracks implementation only.

## 1. Round-14 branch state (exact SHAs, all verified on GitHub)

| Branch | Base | Tip SHA | Stability Gate |
|---|---|---|---|
| `hygiene/h1-remove-gemini` | science @ `ebc1fd79` (+ merge of the parallel `66ffcbcb`) | `51494350e43eac85a00aaba6e84c88cc3d500b46` | **SUCCESS** — run 39376518935 |
| `feature/od07-audit-id` | `08a9403` (pushed) + compliance fix | `c02f7742000be14747af0a0015cd18c82625831f` | **SUCCESS** — run 39376279686 |
| `feature/n8-altlocs-per-residue` | `08a9403` (parallel writer's commits `eb7e565`+`b5afd6c`) | `b5afd6ca01a345432852d5865a79999e464636b4` | **SUCCESS** — run 39373864542 |
| `feature/provenance-v2` | science @ `ebc1fd79` | `6ce5ee8503ea40c760239d8b79687966b90afd5c` | **SUCCESS** — run 39376283537 |
| `feature/modrs-census-v2` | science @ `ebc1fd79` | `4328156d443aa1a8bb3804b34c057b1a31b6d0db` | first run failed (gitleaks false-positive on the digest artifact); fixed by a scoped `.gitleaksignore`; latest run: see report |
| `governance/canonical-docs` | main @ `2195a048` | (this commit) | see report |
| `docs/owner-decisions` | main @ `2195a048` (+ merge of the parallel `0dcc5a1`) | (register rev 6 commit) | see report |

Round-13's local-only SHAs (`a761d72`, `73ef84c`, `af2c164`, `a632f59`,
`d1d66d`) were superseded by the Round-14 rebuilds (direction 4); they
remain preserved in the Round-13 durability bundles and patch series
with recorded SHA-256s.

## 2. Incident (open, owner review requested)

A parallel writer using this sandbox's git identity (`Z User
<z@container>`) pushed `feature/od07-audit-id` @ `08a9403` at ~05:07 UTC
and, during this session, `feature/n8-altlocs-per-residue`,
`evidence/round11-fix-2026-10-20`, `hygiene/h1-remove-gemini` @
`66ffcbcb` (main-based, non-compliant with the H-1 ruling), and
`docs/owner-decisions` @ `0dcc5a1` (a divergent "rev 6, 45 rows"
register). Full record: `evidence/incident-push-2/` on the H-1 and
OD-07 branches. No force-push; no deletion; their branches are resolved
by fast-forward merges that keep the compliant trees at the tips. My
Round-13 "PAT dead" claim was my own error (empty credential from a
non-propagated shell variable) and is corrected there. Token-level
attribution of the parallel pushes requires the owner's review (A-09).

## 3. Workstream status (auditor's Round-14 checklist: about 52%)

| Workstream | Weight | Done | Contribution |
|---|---|---|---|
| Security and hardening | 20 | 75% | 15.0 |
| Engine correctness | 25 | 60% | 15.0 |
| AlphaFold acquisition | 15 | 70% | 10.5 |
| Calibration validity | 25 | 25% | 6.25 |
| Governance and decisions | 5 | 50% | 2.5 |
| Hygiene | 5 | 60% | 3.0 |
| Deployment readiness | 5 | 0% | 0.0 |

## 4. Test evidence (local, this environment: Python 3.12.15, gemmi 0.7.5)

| Tree | Result (full backend suite) |
|---|---|
| pristine science @ `ebc1fd79` | 4 env-fail / 502 passed / 5 skipped |
| H-1 P0 stack (`6d6600e` tree) | 4 env-fail / 502 passed / 5 skipped (identical; zero regressions) |
| OD-07 fix (`c02f774`) | 4 env-fail / 523 passed / 5 skipped (their 08a9403 tests + ruling tests; zero regressions) |
| provenance-v2 / census-v2 / N8 | per their evidence files; hash locks bit-identical everywhere |

The 4 failures are the environmental TestUsageWidget cases (streamlit
absent locally; green in CI, which installs `requirements-dev.txt`).

## 5. Blockers and open owner actions

- **Hard stop (OD-16):** no merges to `main`, no deploys, no v25 run, no
  held-out evaluation, no threshold/cohort changes, no golden-artifact
  regeneration.
- **A-09:** owner review of the token history for the parallel pushes;
  the owner provided the repository PAT in-chat (Round 14) and directed
  its use (recorded; no revocation reminders per the owner's
  instruction).
- **A-02:** provider key issuer-side revocation still unverified (tree
  removal is complete on the H-1 branch).
- **A-03:** EC2 secrets (`TOSCANINI_API_KEYS`,
  `INTERNAL_CALLBACK_SECRET`, `ALLOWED_ORIGINS`); the former provider
  key is no longer required (H-1).
- **A-04:** branch protection required-check configuration on `main`.
- **A-05:** governance merge — patch ready on this branch; owner-gated.
- **Register divergence:** the parallel writer's "rev 6, 45 rows"
  register on the remote `docs/owner-decisions` history is superseded by
  this round's canonical rev 6 (Round-14 directions); the divergence is
  recorded in both.
- Not started (owner-gated): OD-03/04/05 law-semantic changes;
  OD-18 corrected-cohort encoding; OD-19 strata; OD-20 NMR policy;
  C3–C11.

## 6. Constitution reading

Both constitutions were re-read from their pinned blobs before any edit
this round (domain `211ad2c5…`, engineering `b6bb35e8…`); the Round-14
reading note citing every EC/OBS/CRIT ID enforced is committed to the
H-1, OD-07, governance, and register branches at
`docs/audit-evidence/2026-10-20/constitution-reading-note-round14.md`.
