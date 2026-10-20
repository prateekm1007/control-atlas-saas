# Directive index — single source of series status (EC-19)

- **Started:** 2026-10-20 by the EC-16..EC-26 proposal package (branch `proposal/ec16-ec26-enforcement`, base `main` @ `2195a048e87dcec4aa63e8dc7981b3bb6d237ca0`).
- **Status of this file:** rows marked RECONSTRUCTED were backfilled at index creation from the worklog and remote git metadata, not written contemporaneously (EC-18 discipline). Going forward, rows are appended at series open/close time.
- This index is the machine-readable source for: the ACTIVE series (R-DL1), in-force rule IDs and their enforcement points (EC-26), and superseded directives.

## Series register

| Series ID | Source | Opened (UTC) | Status | Items closed | Notes |
|---|---|---|---|---|---|
| SERIES-2026-10-20-R12 | Round-12 auditor directions 1–6 | 2026-10-20 | COMPLETE (RECONSTRUCTED) | 6/6 | Worklog Task 35-round12: reading note; OD-07 (`feature/od07-audit-id`); N8 (`feature/n8-altlocs-per-residue`); H-1 (`hygiene/h1-remove-gemini`); register (`docs/owner-decisions`); evidence-README retraction (`evidence/round11-fix-2026-10-20`). Tips later advanced by SERIES-2026-10-20-R13R14. |
| SERIES-2026-10-20-R13R14 | Round-13/14 audit remediation (commits only; directive text not present in this environment) | 2026-10-20 (commit times 05:35–05:49Z) | COMPLETE on remote (RECONSTRUCTED) | OD-07 canonical triple (`c02f774`); H-1 re-based onto P0 stack (`5149435`); register rev 6 (`0773f0b`) | **Worklog gap:** no worklog entry exists for this cycle; recorded 2026-10-20 from remote evidence (ls-remote + fetched commits). CI: `c02f774` success (runs 118475032870/118475032461), `5149435` success (runs 118475776494/118475776358), `0773f0b` no runs (main trigger defect). |
| SERIES-2026-10-20-EC | Auditor message, 2026-10-20: EC-16…EC-26 amendment proposals + two meta-amendments + adoption order | 2026-10-20 | **ACTIVE** — packaging; PROPOSALS PENDING OWNER APPROVAL | — | Proposal bundle = this branch (`proposal/ec16-ec26-enforcement`). Nothing in force until owner approval. Adoption order (auditor): Batch 1 {EC-17, EC-25}; Batch 2 {EC-22, EC-23}; Batch 3 {rest}. |
| SERIES-2026-10-20-EC-FIX | Auditor direction, 2026-10-20: EC-16…EC-26 adjudications (ACCEPT / ACCEPT WITH REVISIONS, register rev 7) + 9 defect-fix items (3.1–3.9) on this branch | 2026-10-20 | ACTIVE — corrections on this branch; amendments stay PROPOSED/inert until defects fixed AND validation gates pass | 3.1–3.9 | Adjudications recorded in `docs/decisions/OWNER_DECISIONS.md` rev 7 (branch `docs/owner-decisions`). This cycle does NOT activate any workflow and does NOT merge anything. |

## In-force rule set (EC-26)

Until the owner approves the EC-16..EC-26 amendments, the in-force set is unchanged: the EC sections enforced by the pinned constitutions (`ENGINEERING_CONSTITUTION.md` blob `b6bb35e8…`, `constitution.md` blob `211ad2c5…`) plus OBS-153 / CRIT-154 as explicitly out-of-scope — as mapped in the Round-12 reading note (worklog Task 35, first commit of `feature/od07-audit-id`).

Upon approval, each adopted EC-16…EC-26 rule is appended here with its enforcement point (file/command/CI job) per the rule-authoring standard.

## Superseded directives

None recorded in this index. Historical supersessions prior to 2026-10-20 (register revs 1/3/4 rows, D-cycle, etc.) are recorded in `docs/decisions/OWNER_DECISIONS.md` (register rev 6, branch `docs/owner-decisions`) and are not duplicated here.

## Rule inventory (EC-26, machine-readable)

<!-- BEGIN RULE INVENTORY (parsed by scripts/ci/reading_note_map_check.py --index; format below is normative) -->
reading-note-path: docs/evidence/reading-note/READING_NOTE_LATEST.md
pinned-blobs: b6bb35e84956d84a8c9d9e65734c212baf3ed075,211ad2c5d4b40059d9f3f03de12e0701f9235bc7
| Rule ID | Status | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§2 | ACTIVE | script:scripts/ci/reading_note_map_check.py, procedural:AGENTS.md@origin/governance/canonical-docs | AT-06, AT-08 |
| EC-§3 | ACTIVE | workflow:.github/workflows/toscanini-gate.yml@origin/science/p0-correctness-2026-10-19#Secret Scan (gitleaks), pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock |
| EC-§4 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs, pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock |
| EC-§5 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§6 | ACTIVE | script:scripts/ci/tests/run_all.py, procedural:AGENTS.md@origin/governance/canonical-docs | AT-16 |
| EC-§7 | ACTIVE | script:scripts/ci/workflow_deploy_trigger_lint.py, procedural:.github/workflows/toscanini-gate.yml@origin/main | AT-03, AT-14 |
| EC-§8 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§9 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§10 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§11 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§12 | ACTIVE | script:scripts/ci/status_vocabulary_lint.py, procedural:AGENTS.md@origin/governance/canonical-docs | AT-11 |
| EC-§13 | ACTIVE | workflow:.github/workflows/toscanini-gate.yml@origin/science/p0-correctness-2026-10-19#Secret Scan (gitleaks) | workflow:.github/workflows/toscanini-gate.yml@origin/science/p0-correctness-2026-10-19#Secret Scan (gitleaks) |
| EC-§14 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§15 | ACTIVE | procedural:ENGINEERING_CONSTITUTION.md@origin/policy/engineering-constitution | NONE |
| OBS-153 | ACTIVE | procedural:constitution.md | NONE |
| CRIT-154 | ACTIVE | procedural:constitution.md | NONE |
| EC-16 | PROPOSED | script:scripts/ci/verify_claimed_shas.py, workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-17 | PROPOSED | procedural:docs/evidence/push-log/PUSH_LOG.md, script:scripts/ci/evidence_sequence_check.py | AT-13 |
| EC-18 | PROPOSED | script:scripts/ci/evidence_sequence_check.py, workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13 |
| EC-19 | PROPOSED | workflow:proposals/ci/verify-claimed-shas.yml#pr-fields, procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | PROPOSED | script:scripts/ci/verify_claimed_shas.py, workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09 |
| EC-21 | PROPOSED | script:scripts/ci/check_declared_base.py, workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-22 | PROPOSED | script:scripts/ci/status_vocabulary_lint.py, workflow:proposals/ci/status-vocabulary-lint.yml#vocabulary | AT-04, AT-11 |
| EC-23 | PROPOSED | workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks, workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks-all-refs | AT-15 |
| EC-24 | PROPOSED | procedural:docs/evidence/push-log/PUSH_LOG.md, procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md | AT-13 |
| EC-25 | PROPOSED | script:scripts/ci/workflow_deploy_trigger_lint.py, workflow:proposals/ci/no-deploy-on-push-lint.yml#lint, procedural:proposals/remediation/README.md | AT-03, AT-14 |
| EC-26 | PROPOSED | script:scripts/ci/reading_note_map_check.py, workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08 |
<!-- END RULE INVENTORY -->

Status vocabulary: `ACTIVE` = in force now (pinned constitutions / standing decisions); `PROPOSED` =
amendment pending owner approval (register rev 7 records the adjudications); no other tokens are
valid in this table. Enforcement grammar: `script:PATH` / `procedural:PATH` (PATH must exist);
`workflow:FILE#JOB` (JOB must exist in FILE); `pytest:NODE` (the file part of NODE must exist).
Acceptance grammar: `AT-nn` (must exist in `scripts/ci/tests/ACCEPTANCE_TESTS.json`),
`pytest:NODE`, `workflow:FILE#JOB`, or `NONE` (a declared absence — honest, never hidden).
The required set for the reading-note check = all ACTIVE rows (plus PROPOSED rows when the
checker runs with `--include-proposed`, as this proposal branch does for self-testing).
