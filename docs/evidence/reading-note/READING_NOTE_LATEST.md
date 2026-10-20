# Constitution reading note — 2026-10-20, Round-22 acceptance-defect repair

- **Directive series:** SERIES-2026-10-20-R22 (auditor ruling "Round-21 is NOT
  CLOSED / PR #30 is NOT APPROVED FOR MERGE" + directions 1–8; see
  `docs/directives/INDEX.md` and register Part 14).
- **Discipline:** this read was performed BEFORE any repository edit in this
  session (auditor direction 1). The previous worklog and prior audit
  conclusions were not trusted: the live remote state was re-verified by fresh
  ls-remote BEFORE reading and again after the parallel-writer observation
  (below). No retroactive evidence was manufactured.
- **Live state at preflight (2026-10-20T19:06:52Z):** `main` =
  `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` (unchanged); the branch under
  repair was observed at `75f6a8b…` and, minutes later, at
  `908de6dbb9bb183c3c078cbd6be8df1a3475fa32` — two ordinary fast-forward
  pushes (R21B-3 `47022be`, R21B-4 `908de6d`, committed 19:06:36Z/19:12:40Z by
  the same sandbox identity as the documented parallel-writer incident class)
  carrying push-log bookkeeping only (P-069-CONFIRM, P-070-INTENT/CONFIRM,
  P-071 void of P-068, P-072 merge-intent, P-073-INTENT; evidence-table and
  manifest regeneration). No checker, workflow, or documentation content was
  touched by those pushes; this cycle builds forward from `908de6d` with the
  single-writer pre-check re-verified before any push of its own.

## 1. Sources read (verified before reading; Round-22 cycle)

| Document | Pinned blob (full SHA-1) | Verified | Source commit at read time | Read window (UTC, real) |
|---|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES — `git rev-parse` at the candidates under repair returned the identical blob ID | `repair/round19-mainline-acceptance` @ `75f6a8b…` then `908de6d…` | start `2026-10-20T19:06:52Z` (12 lines; OBS-153 and CRIT-154 both ACTIVE) |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES — same method, identical ID; read in full (178 lines, EC §1–§15 + Final Rule) | same tree | same window, immediately following (`19:07:10Z`–`19:07:45Z`) |
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | YES | same tree | read in full (75 lines) |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | YES | same tree | read in full |
| `docs/decisions/OWNER_DECISIONS.md` (register rev 10, Parts 1–13) | `5cbabd6d5e7f41beeb1fabd10c7763e7a2aa37e8` | YES | same tree | read in full; Part 14 appended this cycle |
| `docs/P0_STATUS.md` (sections 1–11) | `d11a86be8ff774ab92ef62a6f59def354bcfb1b1` | YES | same tree | read in full; §12 appended this cycle |
| `docs/directives/INDEX.md` (pre-R22-row) | `a3eaac242128ec9eb10f92bb14c10cb3a208d10c` | YES | same tree | read in full (inventory parsed) |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` (schema v4) | `58b66544c2ffb80d56fe3b51ab429beb3094d0d2` | YES | same tree | read in full; schema v5 defined this cycle |
| `docs/evidence/push-log/PUSH_LOG.md` (P-001..P-073-INTENT at `908de6d`) | `cae6f09d9ac6e59b09d8ec72fb7bd8e02e1c1745` (at `75f6a8b`; extended by the R21B-3/4 appends read in full) | YES | `908de6d…` | read in full |
| `docs/evidence/push-log/RECONCILIATION-2026-10-20-R19-FOLLOWUP-PUSHES.md` | `b1a3df7b25ecfb57c301b7fccbc202575181c10d` | YES | same tree | read in full |
| `docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md` | `a3a6680120c9e95d6e9c13eddddfd477b39b0ccb` | YES | same tree | read in full |
| `docs/evidence/ec16-ec26/R21_EVIDENCE_TABLE.md` + `RECEIPT.md` (R13/R14/R15) | current at `908de6d` | YES | same tree | read in full |

The governed-document reads completed by `2026-10-20T19:25:00Z`; the first
repository writes of this cycle (the checker rewrites) began after that. Every
remote tip above was re-verified by fresh ls-remote before reading and before
writing.

## 2. Rule → enforcement → acceptance mapping (EC-26; the inventory in INDEX.md is normative)

Legend: `workflow:FILE#JOB~events=…~governs=…` = a LIVE GitHub Actions job under
`.github/workflows/` together with the events its enforcement claims and the
changes it governs (Round-20 trigger-liveness rule; Round-21 path-coverage rule;
Round-22 conservative glob-coverage rule — the live workflow's path/branch
filters must cover EVERY structural witness of the governed globs under GitHub's
ORDERED filter semantics, with `!` negations evaluated in order; empty filter
lists and malformed event configuration types FAIL — enforced by
`scripts/ci/reading_note_map_check.py`, AT-23/AT-29/AT-34); `script:PATH` =
checker script; `procedural:PATH` = a tracked document stating the procedure.
`AT-nn` = acceptance test in `scripts/ci/tests/ACCEPTANCE_TESTS.json`. Policy
status and enforcement status are recorded separately: a rule in force is not
proof its checks run, and a check running on a PR branch is not a rule activated
on `main`.

### 2a. ACTIVE rules (in force; live enforcement points after this repair)

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | Purpose and authority; fail closed on ambiguity | procedural:AGENTS.md | NONE |
| EC-§2 | Pre-coding protocol (constitutions first) | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request~governs=docs/evidence/reading-note/**+docs/directives/INDEX.md+scripts/ci/reading_note_map_check.py; procedural:AGENTS.md | AT-06, AT-08, AT-17, AT-23, AT-29, AT-34 |
| EC-§3 | Security is a release requirement (secrets, fail closed; logs MUST not expose credentials) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request~governs=**; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request~governs=**; pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request~governs=** |
| EC-§4 | Product and scientific integrity | procedural:AGENTS.md; pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request~governs=** |
| EC-§5 | Architecture and maintainability | procedural:AGENTS.md | NONE |
| EC-§6 | Testing and evidence (deterministic, exact-SHA results) | script:scripts/ci/tests/run_all.py; workflow:.github/workflows/governance-acceptance.yml#acceptance-suite~events=push+pull_request~governs=scripts/ci/**; procedural:AGENTS.md | AT-16 |
| EC-§7 | CI/CD and supply chain (exact-candidate checks; pinned tooling; PR review before merge) | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request~governs=.github/workflows/**+scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification~events=pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py; procedural:.github/workflows/toscanini-gate.yml | AT-03, AT-14, AT-33 |
| EC-§8 | Reliability and operations | procedural:AGENTS.md | NONE |
| EC-§9 | Customer experience and cost control | procedural:AGENTS.md | NONE |
| EC-§10 | LLM and external-provider boundaries | procedural:AGENTS.md | NONE |
| EC-§11 | Change, release, incident discipline | procedural:AGENTS.md; procedural:docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md; procedural:docs/evidence/push-log/RECONCILIATION-2026-10-20-R19-FOLLOWUP-PUSHES.md | NONE |
| EC-§12 | Definition of done ("not run" never "passed") | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request~governs=docs/evidence/**+scripts/ci/status_vocabulary_lint.py; procedural:AGENTS.md | AT-11, AT-19 |
| EC-§13 | Prohibited shortcuts (no gate weakening; no secrets in logs, reports, or artifacts) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request~governs=** | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request~governs=** |
| EC-§14 | Required task artifacts | procedural:AGENTS.md; procedural:SHA256SUMS (package manifest, machine-checked) | NONE |
| EC-§15 | External engineering references | procedural:ENGINEERING_CONSTITUTION.md | NONE |
| OBS-153 | Avidity quantification (contact density) — ACTIVE domain rule | procedural:constitution.md | NONE |
| CRIT-154 | Structural admissibility (warhead-RMSD, steric veto) — ACTIVE domain rule | procedural:constitution.md | NONE |
| EC-17 | Push on completion (durable storage before subsequent work) | procedural:docs/evidence/push-log/PUSH_LOG.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28, AT-31, AT-32, AT-33 |
| EC-22 | Status vocabulary (RUN_PASS needs a real run ID) | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request~governs=docs/evidence/**+scripts/ci/status_vocabulary_lint.py | AT-04, AT-11, AT-19 |
| EC-25 | No ungated deploy on push | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request~governs=.github/workflows/**+scripts/ci/workflow_deploy_trigger_lint.py; procedural:proposals/remediation/README.md | AT-03, AT-14 |

`NONE` rows are declared absences (honest, not hidden). The live jobs above run
on this branch's pushes/PRs and are NOT live on `main` until PR #30 merges
through the reviewed path and the post-merge verification succeeds.

### 2b. PROPOSED rules (not in force; per the Round-22 direction-8 adjudications)

EC-16, EC-18, EC-19, EC-20, EC-21, EC-23, EC-24, EC-26 remain PROPOSED. Their
checks run live in the governance-acceptance workflow, but RUNNING a check is
not ACTIVATING a rule: each becomes in force only when its activation criteria
pass on the exact final candidate and the auditor records it in the register
(register Part 14.3 carries the Round-22 adjudication table). EC-21's PR-event
failures on #31–#36 are expected while the declared repair base is not on
`main` and must not be weakened or waived.

| Rule ID | Enforcement point(s) (live job on this branch + inert draft) | Acceptance test(s) |
|---|---|---|
| EC-16 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-18 | script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request; workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13, AT-24 |
| EC-19 | workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#pr-fields; procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09, AT-18 |
| EC-21 | script:scripts/ci/check_declared_base.py; workflow:.github/workflows/governance-acceptance.yml#declared-base~events=pull_request; workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-23 | script:scripts/ci/secret_baseline_check.py; script:scripts/ci/secret_control_check.py; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request~governs=**; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks-all-refs | AT-15, AT-21, AT-22, AT-25, AT-30 |
| EC-24 | procedural:docs/evidence/push-log/PUSH_LOG.md; procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification~events=pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28, AT-31, AT-32, AT-33 |
| EC-26 | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request~governs=docs/evidence/reading-note/**+docs/directives/INDEX.md+scripts/ci/reading_note_map_check.py; workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08, AT-17, AT-23, AT-29, AT-34 |

## 3. Round-22 context recorded at read time

- The auditor's independently verified state at ruling time: `main` =
  `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` unchanged; PR #30 open at head
  `0fe4840d1f51bf2c375115116f7e46928fd7ee78` (15 commits ahead, 0 behind; no
  scientific-backend or calibration changes in the diff); Stability Gate runs
  `39473614585` (push) / `39473669119` (PR) and Governance Acceptance runs
  `39473614599` (push, 7+2 jobs) / `39473669306` (PR, 9/9 jobs) green on that
  head; the terminal `P-060-INTENT` and `P-061` merge-intent bound to the
  previous candidate; the ONLY review recorded on PR #30 an automated
  COMMENTED review on the older `2dcbbd781a` head — no current-head approval.
  Merge NOT approved; green CI explicitly ruled necessary-but-not-sufficient.
- Provenance note (direction 3 in practice): the events-API search (full
  retained window, 300 events, 2026-10-18T16:00:59Z..2026-10-20T19:13:33Z)
  found no PushEvent whose `after` equals either P-011's or P-013's SHA, and
  zero Actions runs exist on both heads (the pre-gate-trigger era). For
  P-013's SHA the search found a retained PushEvent (24612279632,
  2026-10-20T10:32:21Z) on the claimed ref whose `before` field IS that SHA
  — a later push's own record that the branch pointed there (prior-tip
  citation, schema v5). For P-011's SHA no ref-transition evidence remains
  in the retained window: it is classified NOT PROVEN by an explicit
  reclassification, never counted as evidence.
- Feed-reality adaptation (disclosed): live verification on 2026-10-20 showed
  this repository's events feed ELIDES PushEvent.after (only `before` is
  reliable) and that ref updates made through the git-data API generate NO
  PushEvent (the last four repair pushes are proven only by their push-event
  Actions runs). The v5 binding therefore accepts, in preference order: a
  full before+after PushEvent match; a before-matching PushEvent plus a
  push-event run on the exact candidate; or a push-event run on the exact
  candidate alone (disclosed). Every decidable check is enforced; no
  PR-event run ever substitutes for push-event provenance.
- The Round-22 blockers repaired in this cycle: EC-23 report-derived metadata
  disclosure (RuleID/File/StartLine interpolation into diagnostics); EC-24
  commit-existence-vs-push-provenance confusion (P-062/P-064 prove only
  existence; P-063/P-065 cite failed push-event runs); missing event/
  head_branch validation; unscoped/duplicate intent reconciliations and the
  branch-local unmatched-intent invariant; expected_remote_tip_before_push
  not bound to push-event before/after state; the merge transition's missing
  `pull-requests: read` permission, missing merge_commit_sha/source-ref/
  merged_at/fresh-approval verification, and absent pre-merge gate; EC-26
  representative-sample coverage and unordered `!` handling. See
  `docs/P0_STATUS.md` §12 and register Part 14.
- The pre-merge-verification job introduced by this cycle is EXPECTED TO FAIL
  on PR events while no independent approval exists on the exact head — that
  is the functioning pre-merge gate per auditor direction 5, not a defect; it
  passes once the required non-author approval is recorded and the check
  re-runs.
- Hard stop (OD-16) remains in force for everything outside this repair's
  authorized scope: no deploys, no v25 calibration, no held-out evaluation,
  no threshold/cohort changes, no law-semantic changes, no golden-artifact
  regeneration. A-01 remains ACCEPTED RISK (no rotation reminder, no rotation
  initiation, no credential operations; the token value is never written
  anywhere). The 40-row cohort metadata table, the independent r2 archive
  verification, A-02 issuer-side revocation, and the EC2 deployment secrets
  remain OPEN. A-04 (the exact branch-protection required-check
  configuration) remains an open blocker until independently verified.
