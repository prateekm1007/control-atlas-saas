# Constitution reading note — 2026-10-20, Round-21 acceptance-defect repair

- **Directive series:** SERIES-2026-10-20-R21 (auditor ruling "PR #30 is not approved
  for merge yet" + directions 1–10; see `docs/directives/INDEX.md`).
- **Discipline:** this read was performed BEFORE any repository edit in this session
  (auditor direction 2). The Round-20 reading note is NOT evidence for this cycle's
  preflight; both constitutions were re-read in full from the exact pinned blobs at
  the SHAs below, with `date -u` timestamps captured around the read window. The
  previous worklog and prior audit conclusions were not trusted: the live remote
  state was re-verified by fresh ls-remote before reading (`main`
  `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`, PR #30 head
  `62fab8a5ad0c07c71a508cb34b93e1a4af24e73e`, matching the auditor's Round-21 table),
  and the candidate worktree was verified clean at `62fab8a…`. No retroactive
  evidence was manufactured.

## 1. Sources read (verified before reading; Round-21 cycle)

| Document | Pinned blob (full SHA-1) | Verified | Source commit at read time | Read window (UTC, real) |
|---|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES — `git rev-parse 62fab8a:constitution.md` returned the identical blob ID (both pinned blobs also re-verified against `origin/main` by the preflight script) | candidate `repair/round19-mainline-acceptance` @ `62fab8a5ad0c07c71a508cb34b93e1a4af24e73e` (remote-verified 2026-10-20T17:39:55Z) | start `2026-10-20T17:40:00Z` (11 lines; OBS-153 and CRIT-154 both ACTIVE) |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES — same method, identical ID; read in full (178 lines, EC §1–§15 + Final Rule) | same tree | same window, immediately following (`17:40:00Z`–`17:40:06Z`) |
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | YES | same tree | read in full (75 lines) |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | YES | same tree | read in full |
| `docs/decisions/OWNER_DECISIONS.md` (register rev 9, pre-Part-13) | `d05bf86f33e74de95b9570522f026263d690b47c` | YES | same tree | read in full (Parts 1–12; Part 13 appended this cycle as rev 10) |
| `docs/P0_STATUS.md` (pre-section-10) | `e646a9e6972764bea2568a239799894a8bd4978c` | YES | same tree | read in full (sections 1–9) |
| `docs/directives/INDEX.md` (pre-R21-row) | `f95dc52e939cf419f2a2c5c8c17b9962cf511c6f` | YES | same tree | read in full (inventory parsed) |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` (schema v3) | `88f3511f267a62e961277786ad788a7ef318af78` | YES | same tree | read in full; schema v4 (identifier pairing, evidence-back, intent-reconcile, merge-intent) defined this cycle |
| `docs/evidence/push-log/PUSH_LOG.md` (pre-R21 appends) | `b1a61144222bee65832a86d3efad37ab66b916c1` | YES | same tree | read in full (33 entries: P-001..P-011 legacy, P-011..P-015-CONFIRM, P-017, P-050..P-059) |
| `docs/evidence/push-log/RECONCILIATION-2026-10-20-R19-FOLLOWUP-PUSHES.md` | `b1a3df7b25ecfb57c301b7fccbc202575181c10d` | YES | same tree | read in full |
| `docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md` | `a3a6680120c9e95d6e9c13eddddfd477b39b0ccb` | YES | same tree | read in full |

The governed-document reads completed by `2026-10-20T17:45:00Z`; the first
repository writes of this cycle (the checker rewrites and records) began after
`17:47Z`. Every remote tip in the auditor's Round-21 table was re-verified
before editing.

## 2. Rule → enforcement → acceptance mapping (EC-26; the inventory in INDEX.md is normative)

Legend: `workflow:FILE#JOB~events=…~governs=…` = a LIVE GitHub Actions job under
`.github/workflows/` together with the events its enforcement claims and the
changes it governs (Round-20 trigger-liveness rule: the live workflow's `on:`
block must actually contain each declared event with branch coverage of `main`;
Round-21 path-coverage rule: the workflow's `paths`/`paths-ignore` filters must
not exclude any governed path, and filters of an invalid type FAIL — enforced by
`scripts/ci/reading_note_map_check.py`, AT-23/AT-29); `script:PATH` = checker
script; `procedural:PATH` = a tracked document stating the procedure. `AT-nn` =
acceptance test in `scripts/ci/tests/ACCEPTANCE_TESTS.json`. Policy status and
enforcement status are recorded separately: a rule in force is not proof its
checks run, and a check running on a PR branch is not a rule activated on
`main`.

### 2a. ACTIVE rules (in force; live enforcement points after this repair)

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | Purpose and authority; fail closed on ambiguity | procedural:AGENTS.md | NONE |
| EC-§2 | Pre-coding protocol (constitutions first) | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request~governs=docs/evidence/reading-note/**+docs/directives/INDEX.md+scripts/ci/reading_note_map_check.py; procedural:AGENTS.md | AT-06, AT-08, AT-17, AT-23, AT-29 |
| EC-§3 | Security is a release requirement (secrets, fail closed; logs MUST not expose credentials) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request~governs=**; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request~governs=**; pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request~governs=** |
| EC-§4 | Product and scientific integrity | procedural:AGENTS.md; pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request~governs=** |
| EC-§5 | Architecture and maintainability | procedural:AGENTS.md | NONE |
| EC-§6 | Testing and evidence (deterministic, exact-SHA results) | script:scripts/ci/tests/run_all.py; workflow:.github/workflows/governance-acceptance.yml#acceptance-suite~events=push+pull_request~governs=scripts/ci/**; procedural:AGENTS.md | AT-16 |
| EC-§7 | CI/CD and supply chain (exact-candidate checks; pinned tooling) | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request~governs=.github/workflows/**+scripts/ci/workflow_deploy_trigger_lint.py; procedural:.github/workflows/toscanini-gate.yml | AT-03, AT-14 |
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
| EC-17 | Push on completion (durable storage before subsequent work) | procedural:docs/evidence/push-log/PUSH_LOG.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28 |
| EC-22 | Status vocabulary (RUN_PASS needs a real run ID) | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request~governs=docs/evidence/**+scripts/ci/status_vocabulary_lint.py | AT-04, AT-11, AT-19 |
| EC-25 | No ungated deploy on push | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request~governs=.github/workflows/**+scripts/ci/workflow_deploy_trigger_lint.py; procedural:proposals/remediation/README.md | AT-03, AT-14 |

`NONE` rows are declared absences (honest, not hidden). The live jobs above run
on this branch's pushes/PRs and are NOT live on `main` until PR #30 merges
through the reviewed path and the post-merge verification succeeds.

### 2b. PROPOSED rules (not in force; per the Round-21 direction-10 adjudications)

EC-16, EC-18, EC-19, EC-20, EC-21, EC-23, EC-24, EC-26 remain PROPOSED. Their
checks run live in the governance-acceptance workflow, but RUNNING a check is
not ACTIVATING a rule: each becomes in force only when its activation criteria
pass on the exact final candidate and the auditor records it in the register
(register rev 10 Part 13.3 carries the Round-21 adjudication table). EC-21's
PR-event failures on #31–#36 are expected while the declared repair base is
not on `main` and must not be weakened or waived.

| Rule ID | Enforcement point(s) (live job on this branch + inert draft) | Acceptance test(s) |
|---|---|---|
| EC-16 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-18 | script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request; workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13, AT-24, AT-26, AT-27, AT-28 |
| EC-19 | workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#pr-fields; procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09, AT-18 |
| EC-21 | script:scripts/ci/check_declared_base.py; workflow:.github/workflows/governance-acceptance.yml#declared-base~events=pull_request; workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-23 | script:scripts/ci/secret_baseline_check.py; script:scripts/ci/secret_control_check.py; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request~governs=**; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks | AT-15, AT-21, AT-22, AT-25 |
| EC-24 | procedural:docs/evidence/push-log/PUSH_LOG.md; procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28 |
| EC-26 | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request~governs=docs/evidence/reading-note/**+docs/directives/INDEX.md+scripts/ci/reading_note_map_check.py; workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08, AT-17, AT-23, AT-29 |

## 3. Round-21 context recorded at read time

- The auditor's independently verified state at ruling time: `main` =
  `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` unchanged; PR #30 open at head
  `62fab8a5ad0c07c71a508cb34b93e1a4af24e73e` with Governance Acceptance run
  `39459283617` and Stability Gate run `39459283515` green on that exact head;
  the required PR approval NOT satisfied (the recorded Codex review is
  COMMENTED, not APPROVED); merge blocked.
- The three critical defects (EC-23 secret-value disclosure; EC-24 unbound
  terminal intents and unproven historical classifications), the unproven
  protected-main merge transition, the EC-26 path-filter gap, and the status
  misalignments are the work items of this cycle (fixed in this candidate; see
  `docs/P0_STATUS.md` §10 and register rev 10 Part 13).
- Hard stop (OD-16) remains in force for everything outside this repair's
  authorized scope: no deploys, no v25 calibration, no held-out evaluation,
  no threshold/cohort changes, no law-semantic changes, no golden-artifact
  regeneration. A-01 remains ACCEPTED RISK (no rotation reminder, no rotation
  initiation, no credential operations; the token value is never written
  anywhere). The 40-row cohort metadata table, the independent r2 archive
  verification, A-02 issuer-side revocation, and the EC2 deployment secrets
  remain OPEN (direction 10).

---

# Round-21-followup reading note (SERIES-2026-10-20-R21B, 2026-10-20)

- **Session:** continuation cycle on `repair/round19-mainline-acceptance`; the
  candidate moved from `62fab8a…` to `0fe4840…` (the Round-21 push) during this
  session, and the preflight below was performed for BOTH candidates (the
  pinned blobs are content-addressed and byte-identical on both).
- **Fetch/read evidence:** read-only fetch of `origin` at 2026-10-20T~17:44Z;
  worktree at `62fab8a` created 17:46:57Z (worktree file mtime); BOTH
  constitutions read in full 17:47–17:50Z (Read tool; probe extraction
  17:50:27Z); `git hash-object` re-verified both pinned blobs on the `0fe4840`
  checkout at ~18:39Z: `constitution.md` =
  `211ad2c5d4b40059d9f3f03de12e0701f9235bc7`, `ENGINEERING_CONSTITUTION.md` =
  `b6bb35e84956d84a8c9d9e65734c212baf3ed075` (both match the pinned blob SHAs
  recorded in the INDEX inventory; timestamps are session-log approximations,
  disclosed as such). First repository edit of this cycle: the conftest probe
  placement 17:50:27Z (probe), first substantive edit ~18:2xZ (checker).
- **Files read in full before editing:** `AGENTS.md`
  (`85767237399ffb30f3e45894d9dfbe6f4240584f`), `docs/coder/BRIEF.md`
  (`052a84ca829d321dc0841a4b1c8702145dcbf6b5`),
  `docs/decisions/OWNER_DECISIONS.md` (rev 9 → rev 10 Part 13 on `0fe4840`),
  `docs/P0_STATUS.md` (`e646a9e…` at `62fab8a`; §10 read at `0fe4840`),
  `docs/directives/INDEX.md` (`f95dc52e…` at `62fab8a`; R21 row read at
  `0fe4840`), `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` (v3 and v4),
  `docs/evidence/push-log/PUSH_LOG.md` (P-050..P-065 + schema records),
  `RECONCILIATION-2026-10-20-R19-FOLLOWUP-PUSHES.md`, the Round-21 evidence
  table, and the R21 state of every checker/test touched.
- **EC/OBS/CRIT IDs enforced this cycle:** EC §1 (purpose/authority; conflicts
  recorded, not edited away), EC §2 (mandatory pre-coding protocol; this note),
  EC §3 (fail-closed secret handling; no credential values anywhere; the
  OWNER-PAT-1 value never written to any file), EC §4 (integrity: synthetic
  fixtures stay test-only; measured vs inferred distinguished), EC §5 (smallest
  complete change; duplicate-entry-point search before the draft repair), EC §6
  (testing and evidence: deterministic tests; exact commands/SHAs/counts; the
  unsound-green defect class eliminated), EC §7 (feature branch + reviewed PR;
  required checks on the exact candidate), EC §11 (change/incident discipline;
  the sandbox read anomaly disclosed, not smoothed), EC §12 (definition of
  done: this cycle's evidence below), EC §13 (no weakened/deleted checks — one
  check ADDED at AT-23), EC §14 (task artifacts: P0_STATUS §11, RECEIPT R15,
  this note), EC-17 (two-phase push records P-060-CONFIRM/P-066-INTENT),
  EC-18/EC-24 (append-only log; remote verification; evidence-field branch
  binding), EC-19 (this INDEX-anchored status), EC-20 (executed acceptance
  tests recorded with counts), EC-22 (status vocabulary respected), EC-25 (no
  deploy-trigger changes), EC-26 (rule inventory consistent; AT IDs valid),
  OBS-153 and CRIT-154 (untouched; no scientific claim made or altered).
- **Verdict impact: none.** This cycle repairs the acceptance suite's
  soundness and the evidence bookkeeping only; no scientific artifact, law,
  threshold, cohort, calibration, or engine verdict is touched. Hash locks
  (canon/matrix/law-semantic) are byte-identical: the files are not referenced
  by this change (verified by the changed-file manifest in the commit).
- **Verification performed on the exact candidate** (pristine `0fe4840`
  checkout + this cycle's changes): pytest suite with the conftest fixture
  active `29 passed` × 3 consecutive runs; `run_all.py` green; reading-note +
  inventory check PASS (20 required rules mapped); EC-25 trigger lint PASS;
  status vocabulary PASS on the four CI-scanned evidence files; the hardened
  evidence-sequence checker PASSES the real push log with the CI-equivalent
  invocation (`--repo` + API token + `--branch`/`--expect-tip`); parent-prefix
  append-only verified at commit time.
