# Constitution reading note — 2026-10-20, EC-16..EC-26 correction cycle

- **Directive series:** SERIES-2026-10-20-EC (correction phase; see `docs/directives/INDEX.md`).
- **Discipline:** this read was performed BEFORE any repository edit in this session (worklog-recorded;
  scratch timestamps captured outside the repo at read time). No retroactive evidence was manufactured.

## 1. Sources read (verified before reading)

| Document | Pinned blob (full SHA-1) | Verified | Source commit at read time | Read window (UTC, real) |
|---|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES — `git cat-file blob <id>` then `git hash-object` returned the identical ID | `main` @ `2195a048e87dcec4aa63e8dc7981b3bb6d237ca0` | start `2026-10-20T09:05:44Z` → end `2026-10-20T09:06:36Z` |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES — same method, identical ID | `policy/engineering-constitution` @ `4993cb3b6954ffb8ecad281817b2c59f1f8fc041` | same window; read in full (178 lines / 20,572 bytes) |
| `AGENTS.md` (entry point) | `85767237399ffb30f3e45894d9dfbe6f4240584f` | YES — `git ls-tree` + `git show` from the governance branch | `governance/canonical-docs` @ `f368bd2ad6bca7b72360d3c390771f0496d36aca` | same window; read in full (75 lines / 4,021 bytes) |

All three documents were read completely in this session, in full, before the first edit
(`date -u` captured at read start and read end into an out-of-repo evidence file; the session's
first repository write happened after the read end above).

## 2. Rule → enforcement → acceptance mapping (EC-26)

Legend — enforcement values: `workflow:FILE#JOB` = GitHub Actions job in a tracked workflow file;
`script:PATH` = checker script (inert until the owning amendment is adopted); `procedural:PATH` =
a tracked document that states the procedure (no automated execution today — declared, not hidden).
Acceptance values: `AT-nn` = test in `scripts/ci/tests/ACCEPTANCE_TESTS.json` (existence enforced);
`workflow:FILE#JOB` = a CI job that must exist; `NONE` = declared absence of an automated test today.

### 2a. ACTIVE rules (in force now)

| Rule ID | Rule (verbatim anchor) | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | Purpose and authority; conflict handling; fail closed on ambiguity | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§2 | Mandatory pre-coding protocol (constitutions read in full first) | script:scripts/ci/reading_note_map_check.py; procedural:AGENTS.md@origin/governance/canonical-docs | AT-06, AT-08 |
| EC-§3 | Security is a release requirement (fail closed, deny by default, secrets) | workflow:.github/workflows/toscanini-gate.yml@origin/science/p0-correctness-2026-10-19#Secret Scan (gitleaks) (main tip's gate currently carries only the unit-tests job — recorded honestly); pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock |
| EC-§4 | Product and scientific integrity (measured vs inferred; law registry) | procedural:AGENTS.md; pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock |
| EC-§5 | Architecture and maintainability (contracts, dependencies, migrations) | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§6 | Testing and evidence (exact command, SHA, counts; never invent) | script:scripts/ci/tests/run_all.py; procedural:AGENTS.md@origin/governance/canonical-docs | AT-16 |
| EC-§7 | CI/CD and supply chain (required checks on exact candidate; pinned actions) | script:scripts/ci/workflow_deploy_trigger_lint.py; procedural:.github/workflows/toscanini-gate.yml | AT-03, AT-14 |
| EC-§8 | Reliability, observability, operations (real health signals; no fake 200) | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§9 | Customer experience, accessibility, cost control | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§10 | LLM and external-provider boundaries (advisory narrative; untrusted input) | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§11 | Change, release, incident discipline (containment → evidence → RCA) | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§12 | Definition of done ("Not run" is never reported as "passed") | script:scripts/ci/status_vocabulary_lint.py | AT-11 |
| EC-§13 | Prohibited shortcuts (no gate-weakening, no secret commits, no fake claims) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks) |
| EC-§14 | Required task artifacts (preflight, manifest, machine-readable results) | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§15 | External engineering references (SSDF / ASVS / SLSA as guidance) | procedural:ENGINEERING_CONSTITUTION.md@origin/policy/engineering-constitution | NONE |
| OBS-153 | Avidity quantification (contact-density standard) — ACTIVE domain rule | procedural:constitution.md | NONE |
| CRIT-154 | Structural admissibility (warhead-RMSD, steric-clash veto) — ACTIVE domain rule | procedural:constitution.md | NONE |

`NONE` rows are declared absences recorded under the proposed rule-authoring standard (R-AS2):
they are honest today and are ADVISORY candidates if and when the owner adopts that standard —
no reclassification is applied while it is still a proposal.

### 2b. EC-16…EC-26 package rules (status per auditor ruling D-5, 2026-10-20)

> **Status update 2026-10-20T11:23:35Z (post-read; source: auditor ruling D-5 of the same day, quoted):** "EC-17 (push on completion), EC-22 (status vocabulary), and EC-25 (no deploy on push) are in force on merge of D-2." The D-2 merge is `b851e01` (gate run 39412416212 = success), so those three are shown `ACTIVE` below. "EC-16, EC-18, EC-19, EC-20, EC-21, EC-23, EC-24, and EC-26 remain proposals. Each becomes in force only when its activation workflow has run green on a feature branch, and I record it in the register" (Direction 6). The read-time content of sections 1–3 is preserved unchanged (EC-18); only status tokens and this notice were added.

| Rule ID | Status | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-16 | PROPOSED | script:scripts/ci/verify_claimed_shas.py; workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-17 | ACTIVE | procedural:docs/evidence/push-log/PUSH_LOG.md; script:scripts/ci/evidence_sequence_check.py | AT-13 |
| EC-18 | PROPOSED | script:scripts/ci/evidence_sequence_check.py; workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13 |
| EC-19 | PROPOSED | workflow:proposals/ci/verify-claimed-shas.yml#pr-fields; procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | PROPOSED | script:scripts/ci/verify_claimed_shas.py; workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09 |
| EC-21 | PROPOSED | script:scripts/ci/check_declared_base.py; workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-22 | ACTIVE | script:scripts/ci/status_vocabulary_lint.py; workflow:proposals/ci/status-vocabulary-lint.yml#vocabulary | AT-04, AT-11 |
| EC-23 | PROPOSED | workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks-all-refs | AT-15 |
| EC-24 | PROPOSED | procedural:docs/evidence/push-log/PUSH_LOG.md; procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md | AT-13 |
| EC-25 | ACTIVE | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:proposals/ci/no-deploy-on-push-lint.yml#lint; procedural:proposals/remediation/README.md | AT-03, AT-14 |
| EC-26 | PROPOSED | script:scripts/ci/reading_note_map_check.py; workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08 |

## 3. What this session's work applies

- The ACTIVE set above governs this session's edits (EC §2 read-first satisfied; EC §6 evidence
  with exact commands/SHAs/counts; EC §12 "not run" never reported as "passed"; EC §13 no secrets).
- The PROPOSED set is implemented as inert corrections under `proposals/` and `scripts/ci/`,
  per the auditor's 2026-10-20 binding directions (register rev 7 records the adjudications).
  Nothing here activates a check or amends a constitution file.
