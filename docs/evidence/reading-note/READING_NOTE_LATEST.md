# Constitution reading note — 2026-10-21, Round-23 sanctioned-actions cycle

- **Directive series:** SERIES-2026-10-21-R23 (auditor ruling serial 3; direction index
  D-R1..D-R6; see `docs/directives/INDEX.md` and register Part 14).
- **Discipline:** this read was performed BEFORE any repository edit in this session
  (auditor direction D-R1). The live remote state was re-verified by fresh ls-remote
  immediately before writing (below), and the push log was re-read in full. The previous
  worklog and prior audit conclusions were not trusted as sole sources: live state and
  pinned blobs were verified directly. No retroactive evidence was manufactured.
- **Live state at preflight (2026-10-20T21:40:40Z):** `main` =
  `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` (unchanged); the repair branch is observed at
  `443d6063b42c40ae47b6ea0607e45d31c12239d2` (55 branch refs, 20 activation/*, unchanged
  since the Round-22 cycle). PR #30 remains open with head synced to this exact SHA.
  Reference census re-verified for this cycle: 55 refs/heads + 34 tags + pull-refs.
- **CI on the exact candidate (fetched from the Actions API at 2026-10-20T20:41:56Z):**
  push-event runs `39492676189` (Stability Gate) and `39492676299` (Governance
  Acceptance), both conclusion=success on head_sha=`443d6063b42…`; PR-event runs
  `39494188745` (Gate, success) and `39494188726` (Governance Acceptance, conclusion=
  failure) on the same head with exactly one failing job: `pre-merge-verification` —
  the registered approval-impossible state on PR #30 (the PR author is the owner
  operator identity holding the credential; GitHub rejects the record of this reviewer's
  approval), adjudicated in register Part 14.

## 1. Sources read (verified before reading; Round-23 cycle)

| Document | Pinned blob (full SHA-1) | Verified | Source commit at read time | Read window (UTC, real) |
|---|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES — blob ID re-hashed from the repair tip tree (git ls-tree) | `repair/round19-mainline-acceptance` @ `443d6063…` | read in full (11 lines; OBS-153 and CRIT-154 both ACTIVE), prior segment of this session and re-opened this cycle |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES — same method, identical ID | same tree | read in full (178 lines, EC §1–§15 + Final Rule); window anchored 2026-10-20T20:4xZ–21:4xZ |
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | YES | same tree | read in full (75 lines) |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | YES | same tree | read in full |
| `docs/directives/INDEX.md` (pre-R23-row) | `244f76eab550951a35305350c76af6ac4ed19b92` | YES | same tree | read in full (inventory parsed; 20 required ACTIVE rules derived; the INDEX rows are NOT edited in this cycle) |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` (schema v5, current) | `83aa798b7a9731e63be2a9664e502b7631a42674` | YES | same tree | read in full |
| `docs/evidence/push-log/PUSH_LOG.md` (P-001..P-078-INTENT at `443d6063`) | `889f6132aeddf3f96bf4ab0e65fb4f2c54679d2d` | YES | same tree | read in full (540 lines, both halves) |
| `docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md` + `RECONCILIATION-…-R19-FOLLOWUP-PUSHES.md` | current at `443d6063…` | YES | same tree | read (carried from the prior session segment; still current at this tip) |
| `.github/workflows/toscanini-gate.yml` | `ff3d66d919666fe234220ed922d636b8d775f249` | YES | same tree | read in full (2 jobs; push+pull_request on `**`; unchanged this cycle) |
| `.github/workflows/governance-acceptance.yml` (pre-deletion) | `44fe653ec4e526725eaea4266dc1807222d3be06` | YES | same tree | read in full (585 lines incl. the `pre-merge-verification` job definition); edited this cycle per D-R2 only |
| `scripts/ci/evidence_sequence_check.py` | `1ae66acf131b59b9b05f8cd06ad9262e8eab74bb` | YES | same tree | key sections read in full (schema/parent-prefix/remote-verify/merge-transition), remainder via targeted search; NOT edited |
| `scripts/ci/reading_note_map_check.py` + `manifest_check.py` + `status_vocabulary_lint.py` + `secret_baseline_check.py` + `secret_control_check.py` + `verify_claimed_shas.py` + `workflow_deploy_trigger_lint.py` | current at `443d6063…` (46-file manifest verified) | YES | same tree | contracts read (this cycle and carried); NOT edited |
| `docs/decisions/OWNER_DECISIONS.md` (register rev 10, Parts 1–14) | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` | YES | same tree | read (Part 14 carries the R22 adjudications incl. 14.4) |
| `docs/P0_STATUS.md` (sections 1–12) | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` | YES | same tree | read (unchanged this cycle) |
| `docs/evidence/ec16-ec26/R21_EVIDENCE_TABLE.md` + `RECEIPT.md` | current at `443d6063…` | YES | same tree | read (carried; not regenerated this cycle — no production change) |

The governed-document reads completed by `2026-10-20T21:45:00Z`; the first repository
writes of this cycle (the D-R2 workflow deletion) began after that. Every remote tip
above was re-verified by fresh ls-remote (2026-10-20T21:40:40Z) before reading and again
immediately before the push. The manifest-check run on the pristine tree
(46/46 digests verified) and the reading-note map validation of the prior note were
executed before any edit was made.

## 2. Rule → enforcement → acceptance mapping (EC-26; the inventory in INDEX.md is normative)

Legend: `workflow:FILE#JOB~events=…~governs=…` = a LIVE GitHub Actions job under
`.github/workflows/` together with the events its enforcement claims and the changes it
governs (Round-20 trigger-liveness rule; Round-21 path-coverage rule; Round-22
conservative glob-coverage rule — the live workflow's path/branch filters must cover
EVERY structural witness of the governed globs under GitHub's ORDERED filter semantics,
with `!` negations evaluated in order; empty filter lists and malformed event
configuration types FAIL — enforced by `scripts/ci/reading_note_map_check.py`,
AT-23/AT-29/AT-34); `script:PATH` = checker script; `procedural:PATH` = a tracked
document stating the procedure. `AT-nn` = acceptance test in
`scripts/ci/tests/ACCEPTANCE_TESTS.json`. Policy status and enforcement status are
recorded separately: a rule in force is not proof its checks run, and a check running on
a PR branch is not a rule activated on `main`.

### 2a. ACTIVE rules (in force; live enforcement points after this cycle's re-enumeration)

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | Purpose and authority; fail closed on ambiguity | procedural:AGENTS.md | NONE |
| EC-§2 | Pre-coding protocol (constitutions first) | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request~governs=docs/evidence/reading-note/**+docs/directives/INDEX.md+scripts/ci/reading_note_map_check.py; procedural:AGENTS.md | AT-06, AT-08, AT-17, AT-23, AT-29, AT-34 |
| EC-§3 | Security is a release requirement (secrets, fail closed; logs MUST not expose credentials) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request~governs=**; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request~governs=**; pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request~governs=** |
| EC-§4 | Product and scientific integrity | procedural:AGENTS.md; pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request~governs=** |
| EC-§5 | Architecture and maintainability | procedural:AGENTS.md | NONE |
| EC-§6 | Testing and evidence (deterministic, exact-SHA results) | script:scripts/ci/tests/run_all.py; workflow:.github/workflows/governance-acceptance.yml#acceptance-suite~events=push+pull_request~governs=scripts/ci/**; procedural:AGENTS.md | AT-16 |
| EC-§7 | CI/CD and supply chain (exact-candidate checks; pinned tooling; PR review before merge) | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request~governs=.github/workflows/**+scripts/ci/workflow_deploy_trigger_lint.py; procedural:.github/workflows/toscanini-gate.yml; procedural:branch protection on main (the platform-side config, independently verified in register Part 13/14: requires the two Stability Gate checks plus one approving review) | AT-03, AT-14, AT-33 |
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

`NONE` rows are declared absences (honest, not hidden). The live jobs above run on this
branch's pushes/PRs and are NOT live on `main` until PR #30 merges through the recorded
path and the post-merge verification succeeds. CHANGE THIS CYCLE (D-R2, sanctioned):
the `governance-acceptance.yml#pre-merge-verification` job is REMOVED from the EC-§7
mapping above (it appeared in the Round-22 note as an enforcement fragment) because the
job is deleted from the live workflow in this same push; see section 3. The
`pull-requests: read` permission and the push-log-integrity canary remain (the
merge-transition verification still requires authenticated PR/review reads).

### 2b. PROPOSED rules (not in force; per the Round-22 direction-8 adjudications)

EC-16, EC-18, EC-19, EC-20, EC-21, EC-23, EC-24, EC-26 remain PROPOSED. Their checks run
live in the governance-acceptance workflow, but RUNNING a check is not ACTIVATING a
rule: each becomes in force only when its activation criteria pass on the exact final
candidate and the auditor records it in the register (register Part 14.3 carries the
Round-22 adjudication table; 14.4 the live events-API re-adjudication). EC-21's PR-event
failures on #31–#36 are expected while the declared repair base is not on `main` and
must not be weakened or waived.

| Rule ID | Enforcement point(s) (live job on this branch + inert draft) | Acceptance test(s) |
|---|---|---|
| EC-16 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-18 | script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request; workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13, AT-24 |
| EC-19 | workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#pr-fields; procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09, AT-18 |
| EC-21 | script:scripts/ci/check_declared_base.py; workflow:.github/workflows/governance-acceptance.yml#declared-base~events=pull_request; workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-23 | script:scripts/ci/secret_baseline_check.py; script:scripts/ci/secret_control_check.py; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request~governs=**; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks-all-refs | AT-15, AT-21, AT-22, AT-25, AT-30 |
| EC-24 | procedural:docs/evidence/push-log/PUSH_LOG.md; procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28, AT-31, AT-32, AT-33 |
| EC-26 | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request~governs=docs/evidence/reading-note/**+docs/directives/INDEX.md+scripts/ci/reading_note_map_check.py; workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08, AT-17, AT-23, AT-29, AT-34 |

(EC-24's row is updated identically to EC-§7: the `#pre-merge-verification` fragment is
removed with the job per D-R2; the schema v5 checks — pairing, binding, merge-transition
verification, reconciliation — continue to run inside `push-log-integrity` on
push+pull_request events, and the checker's offline `--pre-merge-verify` path remains
available to reviewers and to any future maintainer with review authority.)

## 3. Round-23 context recorded at read time

- The auditor's sanctioned actions for this cycle, verbatim intent per the ruling:
  D-R1 read the constitutions and post the reading note on the repair branch;
  D-R2 remove the pre-merge-verification job from the live workflow, keeping every other
  job; D-R3 run the Stability Gate and report the run record on the exact final tip;
  D-R4 merge the repair branch into main through a merge commit (no force-push, no
  squash), documenting exactly what is done and why if the standard PR-merge route
  cannot be completed; D-R5 delete the 25 named branches through the authenticated API;
  D-R6 stop. No product work, no science merges, no register-facing status work beyond
  this record: OD-16 (hard stop on all product work) remains in force.
- The D-R2 deletion is recorded here as SANCTIONED, not as gate weakening: the job's own
  premise ("a fresh approving review from a non-author reviewer applies to the EXACT
  head") is unsatisfiable on this repository while the only credential holder is the
  PR's author (register Part 14; A-04-family finding). The branch protection on main
  continues to require one approving review plus the two Stability Gate checks, so the
  platform-side merge gate itself is untouched; the removal removes only the
  governance workflow's duplicate of an already-platform-enforced and currently
  unsatisfiable precondition. The checker's merge-transition logic in
  `evidence_sequence_check.py` (merge-commit parents, merged PR, source ref,
  merged_at, non-author review verification where available) is NOT edited: if the
  merge executes through the API without a recorded approval, the post-merge
  governance run will disclose that classification honestly instead of a passing one.
- Round-22 state carried into this cycle: `main` unchanged at
  `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`; the terminal `P-078-INTENT` record is
  PENDING on the exact candidate (the designed, truthful state; the push carrying this
  note delivers its CONFIRM); the active merge-intent `P-077` resolved to
  `443d6063b42c40ae47b6ea0607e45d31c12239d2` (its introducing commit) and is voided and
  redeclared as `P-080` on this push under the D-R4 authorization (the merge executes
  per the auditor's ruling rather than the Round-22 approval-gated path, which cannot
  be completed).
- Local environment disclosure (EC §6/§12): the acceptance suite on this operator box
  reports 34 checks passed with 8 conftest-class fixture-writing errors
  (`tmpfile()` writes `os.linesep`, converting `\n` to `\r\n` on Windows and breaking
  the byte-exact fixture comparisons — e.g. AT-13's expected 64063-vs-63601 delta is
  exactly the newline conversion of the same 462-line file). This is an
  environment artifact of the local run, NOT a repository defect: the same suite on
  the same exact SHA concludes SUCCESS on the hosted runners (push run
  `39492676299`, created 2026-10-20T20:31:18Z), and no test code was edited. CI
  remains the authoritative execution environment for the gate report.
- Direction D-R5 consequence, disclosed BEFORE it occurs: deleting the 25 named
  branches removes their remote refs; `evidence_sequence_check.py` (the reviewer's own
  fail-closed design, intentionally left untouched) then reports "remote ref not found
  (missing remote verification)" for every CONFIRM/RECONCILE entry bound to a deleted
  branch, so future `push-log-integrity` runs on this branch and later on `main` will
  report those as undeclared/fail-closed findings. This is the reviewer's designed
  failure mode for missing remote state; it will be reported truthfully after the
  deletions, per D-R4's directive to document exactly what was done. The findings
  concern ONLY the deleted branches' historical entries; they do not implicate main's
  content, the Stability Gate, or the merge transition.
- Hard stop (OD-16) remains in force for everything outside this cycle's sanctioned
  scope: no deploys, no v25 calibration, no held-out evaluation, no threshold/cohort
  changes, no law-semantic changes, no golden-artifact regeneration. A-01 remains
  ACCEPTED RISK (no rotation reminder, no rotation initiation, no credential
  operations; the token value is never written anywhere). The 40-row cohort metadata
  table, the independent archive verification, A-02 issuer-side revocation, and the EC2
  deployment secrets remain OPEN. A-04 (branch-protection required-check configuration)
  was verified live on 2026-10-20 and stays as recorded in register Part 14.

## 8. Round-24 cycle - reading record (2026-10-20, SERIAL 4)

**Directive series:** the Round-21B auditor directions (items 1-8) received
while the Round-21B repair was in progress; this cycle re-established the
baseline on the MERGED `main` and repaired the defects still open there.
Register: Part 15 (auditor-adjudicated). Branch: `repair/round24-claim-head-binding`.

**Constitution read (BEFORE any edit of this cycle, in full, both files):**

| Document | Pinned blob (full SHA-1) | Verified | Read window (UTC) |
|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES - re-hashed from `origin/main` with `git rev-parse` | 2026-10-20T23:14:43Z - 23:15:25Z (read in full: 11 lines; OBS-153 and CRIT-154 both ACTIVE) |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES - re-hashed from `origin/main` with `git rev-parse` | 2026-10-20T23:14:43Z - 23:15:25Z (read in full: 178 lines, EC 1-15 + Final Rule) |
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | YES | read in full in the earlier segment of this session (75 lines) |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | YES | read in full in the earlier segment of this session |
| `docs/directives/INDEX.md` | `244f76eab550951a35305350c76af6ac4ed19b92` (pre-R24 row) | YES | read in full (inventory parsed; 20 required ACTIVE rules derived) |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` (schema v5) | `83aa798b7a9731e63be2a9664e502b7631a42674` | YES | read in full |
| `docs/evidence/push-log/PUSH_LOG.md` | current at `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` | YES | read in full (P-001..P-082-INTENT) |
| `docs/decisions/OWNER_DECISIONS.md` (register rev 11, Parts 1-14) | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` | YES | read (Part 14 carries the R22/R23 adjudications) |
| `docs/P0_STATUS.md` (sections 1-12) | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` | YES | read |

**Live state re-verified before the first write of this cycle (2026-10-20T23:0xZ):**
`main` = `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (the merge commit of PR #30; first parent `f535a71b`,
second parent `d81c874eaedc13c6fdbaed5102ad7e638449036d`); `repair/round19-mainline-acceptance` =
`d81c874eaedc13c6fdbaed5102ad7e638449036d`; PR #30 is MERGED and CLOSED with no independent approval on
record (only a bot `COMMENTED` review on `2dcbbd78`). The worktree used by
this cycle is `E:\atlas-work\control-atlas-saas`, OUTSIDE the operator's
home Git repository.

**Rule mapping for this cycle (EC-26; the inventory in INDEX.md is normative):**

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-16 | Claimed-SHA verification bound to the exact candidate | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request | AT-01, AT-12, AT-35 |
| EC-19 | PR-body directive fields (incl. the exact-head declaration) | workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request | AT-01, AT-07, AT-35 |
| EC-20 | the named acceptance tests are EXECUTED on the candidate | script:scripts/ci/verify_claimed_shas.py | AT-09, AT-18, AT-35 |
| EC-17 | push on completion; the log is append-only | script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity | AT-13, AT-24, AT-26, AT-31, AT-32, AT-33 |
| EC-6 | testing and evidence | script:scripts/ci/tests/run_all.py; workflow:.github/workflows/governance-acceptance.yml#acceptance-suite | AT-16 |
| EC-A2 | the pre-coding protocol (constitutions first) | script:scripts/ci/reading_note_map_check.py; this reading record | AT-06, AT-08, AT-34 |
| EC-A12 | definition of done ("not run" is never "passed") | script:scripts/ci/status_vocabulary_lint.py | AT-04, AT-11, AT-19 |

**What this cycle changed (one ordered commit on `repair/round24-claim-head-binding`):**
`scripts/ci/verify_claimed_shas.py` (the exact-head binding: `--expect-head-sha`
supplied explicitly by the workflow; `Claimed-SHAs` must contain that exact
SHA; the `Status-summary` must cite a real, successful run ON that head),
the `claim-verification` job of the governance workflow (the LIVE pull-request
body is read from the API - the event payload is stale once the body is edited
after the trigger), `scripts/ci/evidence_sequence_check.py` (the
introducing-commit pickaxe now runs from the log's own directory: a pathspec
built from the OS path string fell outside the repository where the path
spelling differs from git's, which made the whole gate fail closed on a valid
log and left eight acceptance groups unrunnable on that platform),
`scripts/ci/tests/test_acceptance.py` (AT-35; AT-13's append-only comparison
copy is written in binary mode), the acceptance registry, the INDEX rows for
EC-16/19/20, the push log, register Part 15, P0_STATUS section 13, this record
and the receipt, plus the regenerated evidence table and manifest.

**Verification performed before the push (exact commands, exact counts):**
the full acceptance suite (pytest over `scripts/ci/tests/test_acceptance.py`),
the push-log checker against the LIVE remote with the Actions and events APIs
(`--repo prateekm1007/control-atlas-saas --api-token-env GH_TOKEN`), the
manifest check, the reading-note map check, the status-vocabulary lint, the
deploy-trigger lint and the secret controls. Results are recorded in the
receipt of this cycle.

**Disclosed, not smoothed over:** the merge transition of PR #30 is NOT
verified by this repository's own gate (no independent approval on the pull
request; the declared pre-merge base is no longer the live tip of `main`).
`main`'s Governance Acceptance run `39504181997` is RED on the push-log
integrity job for that reason, and the repair branch's last push (`d81c874`,
run `39507046623`) is RED for the same reason. No approval was invented, no
gate was weakened, and no log entry claims a merge the evidence does not
support. The auditor's decision on the merged state is OPEN (register 15.4).

**Credential handling:** the PAT was supplied only through the process
environment (`GH_PAT`); no credential file was created inside any Git working
tree; no token value appears in any file, log, commit message or evidence
artifact. A-01 remains ACCEPTED RISK (no rotation reminder, no rotation
initiation, no credential operations).

## 9. Round-25 cycle - reading record (2026-10-20/21, SERIAL 5)

**Directive series:** the auditor's Round-25 decision (directions 1-8), read
against the Round-21B directions it carries forward. Register: Part 16
(auditor-adjudicated). Branch: `repair/round24-claim-head-binding`. **Starting commit:** `2ce77e57194702c2b292fc3f3b28cddd2cff15c9`
(the PR #37 head; `main` = `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`).

**Constitution read (BEFORE any edit of this cycle, in full, in the mandated
order: AGENTS.md, constitution.md, ENGINEERING_CONSTITUTION.md,
docs/coder/BRIEF.md, docs/decisions/OWNER_DECISIONS.md, docs/P0_STATUS.md,
docs/directives/INDEX.md, docs/evidence/reading-note/READING_NOTE_LATEST.md,
docs/evidence/push-log/PUSH_LOG.md, docs/evidence/push-log/PUSH_LOG_SCHEMA.md):**

| Document | Blob SHA-1 at the starting commit | Read |
|---|---|---|
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | in full (75 lines) |
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | in full (11 lines; OBS-153, CRIT-154 ACTIVE) |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | in full (178 lines, EC 1-15 + Final Rule) |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | in full (70 lines) |
| `docs/decisions/OWNER_DECISIONS.md` | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` | in full (Parts 1-14) plus this cycle's Part 15/16 |
| `docs/P0_STATUS.md` | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` | in full (sections 1-12) plus this cycle's 13/14 |
| `docs/directives/INDEX.md` | `244f76eab550951a35305350c76af6ac4ed19b92` | in full (the inventory parsed; 20 required ACTIVE rules) |
| `docs/evidence/reading-note/READING_NOTE_LATEST.md` | current at `2ce77e57194702c2b292fc3f3b28cddd2cff15c9` | in full (the Round-23 and Round-24 records) |
| `docs/evidence/push-log/PUSH_LOG.md` | current at `2ce77e57194702c2b292fc3f3b28cddd2cff15c9` | in full (P-001..P-083-INTENT) |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` | `83aa798b7a9731e63be2a9664e502b7631a42674` | in full (schema v5) |

**Rule mapping for this cycle (EC-26; the inventory in INDEX.md is
normative):**

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-16 | claimed-SHA verification bound to the exact candidate and to the required workflow | script:scripts/ci/verify_claimed_shas.py; script:scripts/ci/fetch_pr_body.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request | AT-01, AT-12, AT-35, AT-36, AT-37 |
| EC-19 | the PR-body declaration is read LIVE and fails closed | script:scripts/ci/fetch_pr_body.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request | AT-01, AT-07, AT-35, AT-36 |
| EC-20 | the named acceptance tests are EXECUTED on the candidate | script:scripts/ci/verify_claimed_shas.py | AT-09, AT-18, AT-35, AT-36 |
| EC-17 | push on completion; the log is append-only; the preservation bundle is machine-verifiable | script:scripts/ci/evidence_sequence_check.py; script:scripts/ci/preservation_bundle.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity; workflow:.github/workflows/governance-acceptance.yml#preservation-bundle | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28, AT-31, AT-32, AT-33 |
| EC-24 | the merge transition is proved BEFORE the merge, on the exact head | workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification~events=pull_request; script:scripts/ci/evidence_sequence_check.py | AT-28, AT-31, AT-32, AT-33 |
| EC-A2 | the pre-coding protocol (constitutions first) | script:scripts/ci/reading_note_map_check.py; this reading record | AT-06, AT-08, AT-34 |
| EC-A12 | definition of done ("not run" is never "passed") | script:scripts/ci/status_vocabulary_lint.py | AT-04, AT-11, AT-19 |
| EC-A14 | required task artifacts (the machine-readable, committed artifact set) | script:scripts/ci/manifest_check.py; workflow:.github/workflows/governance-acceptance.yml#preservation-bundle-verify | AT-38 |

**What this cycle changed:** `scripts/ci/fetch_pr_body.py` (new, fail-closed
live body), `scripts/ci/verify_claimed_shas.py` (required-workflow identity),
`scripts/ci/manifest_check.py` (committed-tree semantics), `scripts/ci/
preservation_bundle.py` (new), `scripts/ci/tests/test_acceptance.py` (AT-36,
AT-37, AT-38, the AT-16 count), `scripts/ci/tests/ACCEPTANCE_TESTS.json`,
`.github/workflows/governance-acceptance.yml` (the fail-closed body step, the
run identity in the manifest, the restored pre-merge job, the preservation
bundle jobs, the corrected job label), `docs/directives/INDEX.md`, the push
log, register Part 16, P0_STATUS section 14, this record, the receipt, and the
regenerated evidence table and manifest.

**Disclosed, not smoothed over:** PR #30's merge did not pass this
repository's own merge gate (only a bot COMMENTED review; Governance
Acceptance run 39504181997 FAILED). The recovery PR has no approval either,
so the restored pre-merge gate is RED - the gate functioning, and the
approval is the auditor/owner's to grant. No approval was invented and no
gate was weakened. The mainline governance issue stays OPEN.

## 10. Round-25 verification round - reading record (2026-10-21T09:55:34Z)

**Directive series:** the auditor's Round-25 verification round (directions
1-9, decisions D-U1..D-U7). Register: Part 16.9. Branch: `repair/round24-claim-head-binding`.
**Starting commit:** `bc6fdb5d92a237d421570d1a462199e4cdcabf28` (the PR #37 head); `main` = `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`.

**Constitution read (in full, before any action of this round; the two
governing files were re-opened and re-read in the order AGENTS.md,
constitution.md, ENGINEERING_CONSTITUTION.md, ...):**

| Document | Blob SHA-1 | Verified | Read window (UTC) |
|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | re-hashed from the starting commit with `git rev-parse` | read in full (11 lines; OBS-153 contact-density standard and CRIT-154 warhead-RMSD/steric-clash veto, both ACTIVE) |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | re-hashed from the starting commit with `git rev-parse` | read in full (178 lines, EC 1-15 + Final Rule) |

The rules this round enforces: EC 1 (a preflight attestation is not proof - the
checks must enforce the invariant), EC 3 (fail closed; manage secrets
correctly), EC 6/12 (evidence or it did not happen; "not run" is never
"passed"), EC 7 (required checks on the exact candidate; a check on a different
SHA is not evidence), EC 11 (incidents: contain, preserve evidence, remediate,
prevent regression - never erase audit evidence), EC 12 (the definition of
done), EC 13 (no gate weakening, no invented evidence, no merging a commit
that did not pass the required checks).

**What this round did (and did not do):**

1. The gate diff from the verified pre-R23 blob is produced and saved
   (`evidence/r25/gate-diff.txt`), with the per-job statement for every
   pre-R23 job, and the D-U2 restore is applied: the file is rebuilt as the
   pre-R23 blob with only the intended Round-24/Round-25 additions
   (`evidence/r25/gate-diff-after-restore.txt`). The R25 report wording is
   corrected in register 16.9.
2. The PR facts for #37 and #38 (author login, head SHA, state) are recorded
   in a NEW incident note
   (`docs/evidence/push-log/INCIDENT-2026-10-21-R25-PR-IDENTITY.md`); no old
   incident identifier is reused.
3. The containment check was run and it FAILED: `781175d18f94b65ee09635f2e4f38bf67d545b45` is not an ancestor of
   `bc6fdb5d92a237d421570d1a462199e4cdcabf28` (the branches are siblings from `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`), and the two push logs
   are mutually exclusive. Per the direction, this round STOPS there: PR #38 is
   not closed, no branch is deleted, no archive tag is created, and the merge
   is not attempted (the owner action O-1 - an independent collaborator's
   approval on the exact head - has not happened, and the restored pre-merge
   gate is red for exactly that reason).
4. The scientific/production hard stop is untouched; no product, deployment,
   threshold, cohort or law-semantic work was started.

**Rule mapping for this round (EC-26; the inventory in INDEX.md is normative):**

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-17 | push on completion; the log is append-only; the preservation bundle is machine-verifiable | script:scripts/ci/evidence_sequence_check.py; script:scripts/ci/preservation_bundle.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity; workflow:.github/workflows/governance-acceptance.yml#preservation-bundle | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28, AT-31, AT-32, AT-33 |
| EC-24 | the merge transition is proved BEFORE the merge, on the exact head | workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification~events=pull_request; script:scripts/ci/evidence_sequence_check.py | AT-28, AT-31, AT-32, AT-33 |
| EC-A11 | incident discipline: contain, preserve the evidence, remediate, prevent recurrence | docs/evidence/push-log/INCIDENT-2026-10-21-R25-PR-IDENTITY.md; docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md | NONE |
| EC-A12 | definition of done ("not run" is never "passed") | script:scripts/ci/status_vocabulary_lint.py | AT-04, AT-11, AT-19 |
| EC-A13 | no gate weakening, no invented evidence | workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification | AT-33 |

## 11. Round-26 cycle - reading record (2026-10-21T10:11:06Z, SERIES-2026-10-21-R26)

### 11.1 Sources read (verified by blob ID before reading, this round)

| Document | Pinned blob (full SHA-1) | Verified | Source commit at read time | Read window (UTC) |
|---|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES - `git ls-tree` at `781175d...` and at `63da0e3...` (identical) | `repair/round24-claim-head-binding` @ `f00e3be9f040045fdb9fae7a35f59be1ba25549f` | 2026-10-21T10:11:06Z, read in full (11 lines; OBS-153 and CRIT-154 both ACTIVE) |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES - same method, identical ID | same tree | 2026-10-21T10:11:06Z, read in full (178 lines: EC section 1-15 + Final Rule) |

Both files were re-read in full in this round, not carried over from the previous cycle; the
blob IDs are identical at the R24 candidate, at `main` and at this head.

### 11.2 What is folded into this head (D-V1: one recovery PR, one head)

- The R24 recovery content of PR #38 (head `781175d18f94b65ee09635f2e4f38bf67d545b45`): the
  preservation-bundle generator `scripts/ci/generate_preservation_bundle.py` and its workflow
  `.github/workflows/preservation-bundle.yml` (the generator that produced the published
  artifact, id `12315612470`, archive digest `sha256:29e321c48f4a0990ffcd6f7ae4809519d4dae93bfdb34db7ac4fe3e1eddf9a5c`,
  bundle manifest digest `f8828d7399a54a8e0aa6cf3ba5b7ad9afca26f6eedb62721dfa759d1099fe58a`,
  8/8 checksums verified from a fresh download); the R23 governance-bypass incident record
  `docs/evidence/push-log/INCIDENT-2026-10-21-R23-GOVERNANCE-BYPASS.md`; the gate regression
  coverage (now AT-39 and AT-40, renumbered because AT-35..AT-38 are this head's own
  claim-binding / manifest tests); the R24 recovery status, register and index rows; and the
  push-log records for the R24 candidate (renumbered P-082..P-086 -> P-090..P-090-CONFIRM,
  because this head already used P-082..P-089; the R24 content is preserved, nothing is
  rewritten).
- The pre-merge gate itself: this head's `f00e3be` already restored the pre-R23 gate content,
  so the R24 restoration is subsumed rather than duplicated. AT-39 asserts the gate is an
  executable job (not a comment) with `pull-requests: read`, the authenticated API canary, the
  runs manifest and the exact-head verification; AT-40 asserts it rejects missing, stale,
  non-qualifying and author-only approval evidence and still accepts a fresh non-author
  approval on the exact head.
- The archive-tag resolution (D-V3/D-V4, this round's new control): `scripts/ci/evidence_sequence_check.py`
  now resolves `refs/tags/archive/<branch>` - and the hyphen spelling
  `refs/tags/archive/<branch-with-slashes-as-hyphens>` - as remote evidence for a PRUNED
  branch, using the PEELED commit of the ANNOTATED tag; an archive tag at a different SHA does
  not preserve the confirmation (fail closed), and a tag without a peeled commit is refused.
  AT-41 proves both directions.
- The R24 bundle generator's `merged`-flag correction is folded into this head's producer
  `scripts/ci/preservation_bundle.py`: the PR LIST endpoint returns a null `merged` flag, so
  every PR is re-read from its DETAIL endpoint, which carries the authoritative merged state.
- Two bundle producers now exist and this is disclosed rather than hidden: the R25 in-workflow
  jobs use `scripts/ci/preservation_bundle.py`; `scripts/ci/generate_preservation_bundle.py`
  is retained because it is the byte-exact provenance of PR #38's published artifact. Their
  consolidation is a follow-up, not part of this fold.

### 11.3 State verified live in this round

- `main` is untouched at `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`; this head descends from it.
- 57 branch heads: 55 (the 55 refs of the restore set, incl. the R24 candidate) + this head's
  branch + the parallel writer's branch. The D-V5 keep set (11 branches) is present and
  disjoint from the prune set.
- Four out-of-scope heads are NOT merged into `main` and are archived rather than deleted
  without a tag: `chore/hygiene-2026-10-19`, `docs/owner-decisions-r19`,
  `evidence/round11-fix-2026-10-20`, `evidence/round9-r2-2026-10-20`. Their content is
  preserved by their archive tags at their recorded SHAs.
- `refs/tags/archive/*`: none existed before this round; the archive tags are created by the
  Round-26 pruning steps, each verified with `git ls-remote --tags` against the branch SHA
  before its branch is deleted.
- O-1 (an independent collaborator with write access) is still outstanding: the repository's
  only collaborator is the PR author. PR #37 therefore stays blocked; no self-approval is
  sought or given (D-V1, D-V6, EC-13).

### 11.3b CI evidence for this head (Round-26, disclosed)

- The fold commit 3d0b5ab failed three push-event jobs (run 39583630921: acceptance suite,
  package manifest, push-log). Both causes were fixed in 7071b5d, and only those causes were
  touched: (a) the fold commit staged the regenerated SHA256SUMS before the manifest itself was
  rewritten, so the committed manifest carried stale digest rows for the two evidence files
  edited last; the manifest is now regenerated and staged LAST, with all 51 rows verified
  against the staged blobs before the commit; (b) the terminal push-intent P-093-INTENT had to
  be completed (P-093-CONFIRM) rather than left pending, because the checker compares a
  terminal intent's declared pre-push tip with the newest CONFIRM/RECONCILE recorded for the
  branch BEFORE it in the log - P-094-INTENT then declares the next push on that exact recorded
  state and P-095 redeclares this PR's merge-intent on the exact final candidate.
- The acceptance suite also failed in CI only (it is green locally): the AT-41 fixture created
  its annotated archive tag with `git tag -a` and no explicit identity, which a CI runner has no
  global identity for. The fixture now sets one explicitly - a test-environment fix, no
  weakening.
- The PR-event claim-verification job reads the PR body live, and the body must declare the
  exact head, so it is updated to this head before the push that introduces it; the recorded
  job for a given head is evidence about the body as it stood when that job ran.

### 11.4 Rule mapping additions in this round

| Rule ID | Enforcement point(s) | Acceptance test(s) |
|---|---|---|
| EC-§7 | workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification~events=pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-33, AT-39, AT-40 |
| EC-17 | script:scripts/ci/evidence_sequence_check.py (archive-tag resolution for pruned branches) | AT-13, AT-31, AT-32, AT-33, AT-41 |
| EC-24 | script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification~events=pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-31, AT-32, AT-33, AT-39, AT-40, AT-41 |
