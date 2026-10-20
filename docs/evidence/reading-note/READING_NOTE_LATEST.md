# Constitution reading note — 2026-10-20, Round-20 acceptance-defect repair

- **Directive series:** SERIES-2026-10-20-R20 (auditor ruling "Round-19 is not closed" +
  directions 1–8; see `docs/directives/INDEX.md`).
- **Discipline:** this read was performed BEFORE any repository edit in this session. The
  Round-19 reading note is NOT evidence for this cycle's preflight (auditor direction 2:
  "The Round-19 report records an earlier constitution preflight; it does not replace the
  required preflight for this next edit"); both constitutions were re-read in full from the
  exact pinned blobs at the SHAs below, with `date -u` timestamps captured around the read
  window. No retroactive evidence was manufactured.

## 1. Sources read (verified before reading; Round-20 cycle)

| Document | Pinned blob (full SHA-1) | Verified | Source commit at read time | Read window (UTC, real) |
|---|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES — `git hash-object` on the candidate worktree returned the identical blob ID (both pinned blobs also re-verified against `origin/main`) | candidate `repair/round19-mainline-acceptance` @ `f44bf4d93f8d69a4bd651396f383e8f3ad2b707a` (remote-verified by fresh ls-remote at session start; `main` still `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`) | start `2026-10-20T15:37:43Z` |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES — same method, identical ID; read in full (178 lines) | same tree | same window, immediately following |
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | YES | same tree | read in full (75 lines) |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | YES | same tree | read in full |
| `docs/decisions/OWNER_DECISIONS.md` (register rev 7, pre-port) | `e4f7808d490411e257a523ff8a336515c45028ac` | YES | same tree | read in full (Parts 1–10); the rev-8 register on `docs/owner-decisions-r19` @ `df98f5e…` was then read in full before the content-only port (rev 9, Part 12) |
| `docs/P0_STATUS.md` | `c845034a669cd6c0e7d06c43c6625cbffb8e8e6a` | YES | same tree | read in full (sections 1–8) |
| `docs/directives/INDEX.md` | `f815f51cf3e0391a3f86fca840569088ba8decff` | YES | same tree | read in full (inventory parsed) |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` | schema v2 (two-phase; pre-port) | YES | same tree | read in full; schema v3 (RECONCILE + remote verification) defined this cycle |

The first repository writes of this cycle (the EC-23/EC-24/EC-26 checker
changes, the baseline, and the records) happened after the read window above.
Every remote tip in the auditor's Round-20 table was re-verified by fresh
ls-remote before editing (`main` `f535a71…`, PR #30 head `f44bf4d…`, PRs
#31–#36 heads `bd8d86a…`/`32d272d…`/`151b8c1…`/`97067a0…`/`fa854f8…`/
`28ee983…`, register carrier `df98f5e…`).

## 2. Rule → enforcement → acceptance mapping (EC-26; the inventory in INDEX.md is normative)

Legend: `workflow:FILE#JOB~events=…` = a LIVE GitHub Actions job under
`.github/workflows/` together with the events its enforcement claims (Round-20
trigger-liveness rule: the live workflow's `on:` block must actually contain
each declared event with branch coverage of `main` — enforced by
`scripts/ci/reading_note_map_check.py`, AT-23); `script:PATH` = checker
script; `procedural:PATH` = a tracked document stating the procedure. `AT-nn` =
acceptance test in `scripts/ci/tests/ACCEPTANCE_TESTS.json`. Policy status and
enforcement status are recorded separately: a rule in force is not proof its
checks run, and a check running on a PR branch is not a rule activated on
`main`.

### 2a. ACTIVE rules (in force; live enforcement points after this repair)

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | Purpose and authority; fail closed on ambiguity | procedural:AGENTS.md | NONE |
| EC-§2 | Pre-coding protocol (constitutions first) | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request; procedural:AGENTS.md | AT-06, AT-08, AT-17, AT-23 |
| EC-§3 | Security is a release requirement (secrets, fail closed) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request; pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request |
| EC-§4 | Product and scientific integrity | procedural:AGENTS.md; pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request |
| EC-§5 | Architecture and maintainability | procedural:AGENTS.md | NONE |
| EC-§6 | Testing and evidence | script:scripts/ci/tests/run_all.py; workflow:.github/workflows/governance-acceptance.yml#acceptance-suite~events=push+pull_request; procedural:AGENTS.md | AT-16 |
| EC-§7 | CI/CD and supply chain (exact-candidate checks; pinned tooling) | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request; procedural:.github/workflows/toscanini-gate.yml | AT-03, AT-14 |
| EC-§8 | Reliability and operations | procedural:AGENTS.md | NONE |
| EC-§9 | Customer experience and cost control | procedural:AGENTS.md | NONE |
| EC-§10 | LLM and external-provider boundaries | procedural:AGENTS.md | NONE |
| EC-§11 | Change, release, incident discipline | procedural:AGENTS.md; procedural:docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md; procedural:docs/evidence/push-log/RECONCILIATION-2026-10-20-R19-FOLLOWUP-PUSHES.md | NONE |
| EC-§12 | Definition of done ("not run" never "passed") | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request; procedural:AGENTS.md | AT-11, AT-19 |
| EC-§13 | Prohibited shortcuts (no gate weakening; no secrets) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request |
| EC-§14 | Required task artifacts | procedural:AGENTS.md; procedural:SHA256SUMS (package manifest, machine-checked) | NONE |
| EC-§15 | External engineering references | procedural:ENGINEERING_CONSTITUTION.md | NONE |
| OBS-153 | Avidity quantification (contact density) — ACTIVE domain rule | procedural:constitution.md | NONE |
| CRIT-154 | Structural admissibility (warhead-RMSD, steric veto) — ACTIVE domain rule | procedural:constitution.md | NONE |
| EC-17 | Push on completion (durable storage before subsequent work) | procedural:docs/evidence/push-log/PUSH_LOG.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request | AT-13, AT-20, AT-24 |
| EC-22 | Status vocabulary (RUN_PASS needs a real run ID) | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request | AT-04, AT-11, AT-19 |
| EC-25 | No ungated deploy on push | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request; procedural:proposals/remediation/README.md | AT-03, AT-14 |

`NONE` rows are declared absences (honest, not hidden). Before the Round-19
repair, EC-17/EC-22/EC-25 were in force as POLICY with no live automated
enforcement; the live jobs above run on this branch's pushes/PRs and are NOT
live on `main` until PR #30 merges through the reviewed path.

### 2b. PROPOSED rules (not in force; per the Round-20 direction-7 adjudications)

EC-16, EC-18, EC-19, EC-20, EC-21, EC-23, EC-24, EC-26 remain PROPOSED. Their
checks run live in the governance-acceptance workflow (claim-verification,
push-log-integrity, declared-base, secret-scan-all-refs,
reading-note-inventory jobs), but RUNNING a check is not ACTIVATING a rule:
each becomes in force only when its exact final candidate meets its activation
criteria and the auditor records it in the register (register rev 9 Part 12.6
carries the Round-20 adjudication table verbatim in substance). EC-21's
PR-event failures on #31–#36 are expected while the declared repair base is
not on `main` and must not be weakened or waived.

| Rule ID | Enforcement point(s) (live job on this branch + inert draft) | Acceptance test(s) |
|---|---|---|
| EC-16 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-18 | script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request; workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13, AT-24 |
| EC-19 | workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#pr-fields; procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09, AT-18 |
| EC-21 | script:scripts/ci/check_declared_base.py; workflow:.github/workflows/governance-acceptance.yml#declared-base~events=pull_request; workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-23 | script:scripts/ci/secret_baseline_check.py; script:scripts/ci/secret_control_check.py; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks | AT-15, AT-21, AT-22 |
| EC-24 | procedural:docs/evidence/push-log/PUSH_LOG.md; procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request | AT-13, AT-20, AT-24 |
| EC-26 | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request; workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08, AT-17, AT-23 |

## 3. Round-20 context recorded at read time

- The auditor's independently verified state at ruling time: `main` =
  `f535a71…` unchanged; PR #30 head `f44bf4d…` with gate and Governance
  Acceptance green (runs `39443515998`/`39443516042`) but NOT merge-approved
  (latent acceptance defects; no approving review visible); PRs #31–#36
  blocked on expected EC-21 PR-event failures; register rev 8 prepared on
  `docs/owner-decisions-r19` but outside the reviewed path; no merge approved.
- The two critical deficiencies (EC-23 category exemption + nonzero-exit
  positive control; EC-24 no live-remote-ref verification with the R19
  follow-up pushes unreconciled), the EC-26 trigger-liveness gap, and the
  register gap are the work items of this cycle (fixed in this candidate;
  see `docs/P0_STATUS.md` §9 and register rev 9 Part 12).
- Hard stop (OD-16) remains in force for everything outside this repair's
  authorized scope: no deploys, no v25 calibration, no held-out evaluation,
  no threshold/cohort changes, no law-semantic changes, no golden-artifact
  regeneration. A-01 remains ACCEPTED RISK (no rotation reminder, per the
  owner's standing instruction). The 40-row cohort metadata table and the
  independent r2 archive acceptance remain OPEN (direction 8).
