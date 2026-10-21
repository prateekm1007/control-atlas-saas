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

## LAW-130 branch read record (2026-10-21T18:43:54Z, fix/law-130-clash-boundary)

The constitution was read before any change on this branch. The read is recorded with the exact
blob hashes and the read timestamp, as AGENTS.md and EC-2 require.

| Document | Blob (full SHA-1) | Read (UTC) |
|---|---|---|
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | 2026-10-21T18:43:54Z |
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | 2026-10-21T18:43:54Z |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | 2026-10-21T18:43:54Z |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | 2026-10-21T18:43:54Z |

Branch `fix/law-130-clash-boundary`, based on `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
(a fresh fetch; the branch is isolated from the LAW-160 work so that the two engine corrections
never share a patch). Rules enforced here: EC-2 (constitutions first, preflight, machine-readable
evidence), EC-4 (a failed measurement is never a numerical zero and never a scientific failure),
EC-6 (regression tests first, exact-SHA evidence, "not run" is never "passed"), EC-7 (feature
branch, checks on the exact candidate), EC-11 (append-only records, no rewriting of published
history), EC-12 (definition of done), EC-13 (no weakened gate, no invented status, no
fabricated measurement), EC-14 (artifacts), OD-02/OD-03/OD-05 (the law canon, its thresholds
and the active law set are frozen; this change corrects the implementation of the EXISTING
contract and alters no threshold, operator, scope or law id), and OD-16 (the scientific/production
hard stop: no merge to `main`, no deployment, no v25 run, no held-out evaluation, no threshold or
cohort change, no golden-artifact regeneration).

Scope of this branch: the LAW-130 clashscore measurement only. A failed calculation becomes
NOT_MEASURED (never a PASS with a zero value) and the comparison is brought into exact agreement
with the declared inclusive operator (`<=` at 20.0). The canon, the modality matrix, the
calibration cohort and the golden artifacts are not touched.
