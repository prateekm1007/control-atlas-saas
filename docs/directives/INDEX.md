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
| SERIES-2026-10-20-EC-FIX | Auditor direction, 2026-10-20: EC-16…EC-26 adjudications (ACCEPT / ACCEPT WITH REVISIONS, register rev 7) + 9 defect-fix items (3.1–3.9) on this branch | 2026-10-20 | COMPLETE — corrections delivered on `proposal/ec16-ec26-enforcement` @ `a0e321c` | 3.1–3.9 | Adjudications recorded in `docs/decisions/OWNER_DECISIONS.md` rev 7 (branch `docs/owner-decisions`). That branch was left 18 commits behind main after the release; rebuilt per D-4 as SERIES-2026-10-20-R17 below. |
| SERIES-2026-10-20-R17 | Auditor verdict "partially accepted" + D-1..D-9 + Directions 1–7 (2026-10-20) | 2026-10-20 | ACTIVE — Directions 1–3 verified complete on remote; Direction 4 = this rebuild; Directions 5–6 in progress; Direction 7 = stop | D-1, D-2, D-3 | Release `main` @ `d292a38`, gate run 39411457483 success. Governance tip `d51ad7c` gate 39412154545 success, merged as `b851e01` (gate 39412416212). Register tip `afd445c` gate 39412159331 success, merged as `eb6d6a9` (gate 39412574505). Per D-4 the proposal is rebuilt onto `eb6d6a9` (which contains `d292a38`) as branch `proposal/ec16-ec26-enforcement-rebased`; original branch untouched (no force-push); `.github/` byte-identical to main. D-5: EC-17/EC-22/EC-25 in force upon the D-2 merge; remaining eight PROPOSED pending Direction-6 activation runs. |

## In-force rule set (EC-26)

Per the auditor's 2026-10-20 ruling (D-5, under full delegation): **EC-17 (push on completion), EC-22 (status vocabulary), and EC-25 (no deploy on push) are IN FORCE** upon merge of D-2 (`governance/canonical-docs` merged to main as `b851e01`, gate run 39412416212 = success). They are recorded `ACTIVE` in the rule inventory below; the auditor records this activation in the decision register.

**EC-16, EC-18, EC-19, EC-20, EC-21, EC-23, EC-24, and EC-26 remain PROPOSED**: each becomes in force only when its activation workflow has run green on a feature branch AND the auditor records it in the register (Direction 6 of the same ruling).

The base in-force set is unchanged: the EC sections enforced by the pinned constitutions (`ENGINEERING_CONSTITUTION.md` blob `b6bb35e8…`, `constitution.md` blob `211ad2c5…`) plus OBS-153 / CRIT-154 — as mapped in the Round-12 reading note (worklog Task 35) and in `docs/evidence/reading-note/READING_NOTE_LATEST.md`.

| SERIES-2026-10-20-R19 | Auditor Round-19 directions 1-8 (2026-10-20): mainline repair-forward from `f535a71`, push-log collision/deviation reconciliation, live end-to-end enforcement, activation adjudications | 2026-10-20 | COMPLETE (repair + six activation candidates + register rev 8 carrier pushed; PRs #30-#36 open) — NOT closed by the Round-20 audit; superseded by SERIES-2026-10-20-R20 for the acceptance defects below | D2 (adjudication applied), D1, D3a, D3b, D4 (records), D5 (root cause + fix), D6 (EC-26 fix) | Branch `repair/round19-mainline-acceptance` from exact `main` @ `f535a71`. SDK tokens removed; README corrected; manifest scoped + regenerated; `INCIDENT-2026-10-20-P015-COLLISION.md` recorded; EC-20 collection fix (runs 39417485854/39417488957/39417492355); EC-26 ACTIVE-liveness fix; live `governance-acceptance.yml`. Auditor Round-20 findings on this series: (a) the all-refs classification carried a broad docs//evidence/ category exemption; (b) the positive control accepted any nonzero scanner exit; (c) the push-log checker did not verify recorded remote_sha values against live remote refs, and the R19 follow-up pushes (tips `f44bf4d`, `bd8d86a`, `32d272d`, `151b8c1`, `97067a0`, `fa854f8`, `28ee983`, `df98f5e`) were themselves undeclared by INTENT records — reconciled truthfully (RECONCILE records, never backdated) in the R20 series. |
| SERIES-2026-10-20-R20 | Auditor Round-20 ruling "Round-19 is not closed" + directions 1-8 (2026-10-20): EC-23 fingerprint-allowlist + genuine-detection controls, EC-24 remote-ref verification, EC-26 trigger liveness, register rev 8 into the reviewed path, truthful status records | 2026-10-20 | ACTIVE — fixes implemented on `repair/round19-mainline-acceptance`; merge of PR #30 only after the complete gate is green on the exact final candidate AND the required independent review (GitHub reports no approving review; do not bypass branch protection) | EC-23 (AT-21, AT-22), EC-24 (AT-24), EC-26 (AT-23), register rev 8 port, P0_STATUS §9, RECONCILE records P-050..P-057 | Two critical acceptance defects fixed: the docs//evidence/ category exemption replaced by the reviewed fingerprint allowlist (`docs/evidence/ec16-ec26/secret_scan_baseline.json`, 47 documented entries; H-1 class separate under A-02; unknown fingerprints fail wherever they live); the positive control now verifies the expected finding, the fired rule (`aws-access-token`), and key redaction, failing closed on scanner execution errors (`--exit-code 101` disambiguation). EC-24: every CONFIRM/RECONCILE verified against the live remote ref (CURRENT/HISTORICAL classification; undeclared branch-under-test tips FAIL; missing confirmation = PENDING, disclosed, never backdated; cited run IDs checked against the Actions runs manifest). EC-26: liveness now requires the declared events (`~events=` suffix) to actually be present in the live workflow's `on:` block with branch coverage of `main` (manual-only / branch-disabled / untriggered workflows FAIL). Register rev 8 carried into the reviewed path content-only (no history merge of the divergent register branch). Nothing activated; no merge approved by this series. |

## Superseded directives

None recorded in this index. Historical supersessions prior to 2026-10-20 (register revs 1/3/4 rows, D-cycle, etc.) are recorded in `docs/decisions/OWNER_DECISIONS.md` (register rev 6, branch `docs/owner-decisions`) and are not duplicated here.

## Rule inventory (EC-26, machine-readable)

<!-- BEGIN RULE INVENTORY (parsed by scripts/ci/reading_note_map_check.py --index; format below is normative) -->
reading-note-path: docs/evidence/reading-note/READING_NOTE_LATEST.md
pinned-blobs: b6bb35e84956d84a8c9d9e65734c212baf3ed075,211ad2c5d4b40059d9f3f03de12e0701f9235bc7
| Rule ID | Status | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§2 | ACTIVE | script:scripts/ci/reading_note_map_check.py, workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request, procedural:AGENTS.md | AT-06, AT-08, AT-17, AT-23 |
| EC-§3 | ACTIVE | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request, workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request, pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request |
| EC-§4 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs, pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request |
| EC-§5 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§6 | ACTIVE | script:scripts/ci/tests/run_all.py, workflow:.github/workflows/governance-acceptance.yml#acceptance-suite~events=push+pull_request, procedural:AGENTS.md | AT-16 |
| EC-§7 | ACTIVE | script:scripts/ci/workflow_deploy_trigger_lint.py, workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request, procedural:.github/workflows/toscanini-gate.yml | AT-03, AT-14 |
| EC-§8 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§9 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§10 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§11 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§12 | ACTIVE | script:scripts/ci/status_vocabulary_lint.py, workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request, procedural:AGENTS.md | AT-11, AT-19 |
| EC-§13 | ACTIVE | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request |
| EC-§14 | ACTIVE | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§15 | ACTIVE | procedural:ENGINEERING_CONSTITUTION.md@origin/policy/engineering-constitution | NONE |
| OBS-153 | ACTIVE | procedural:constitution.md | NONE |
| CRIT-154 | ACTIVE | procedural:constitution.md | NONE |
| EC-16 | PROPOSED | script:scripts/ci/verify_claimed_shas.py, workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request, workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-17 | ACTIVE | procedural:docs/evidence/push-log/PUSH_LOG.md, script:scripts/ci/evidence_sequence_check.py, workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request | AT-13, AT-20, AT-24 |
| EC-18 | PROPOSED | script:scripts/ci/evidence_sequence_check.py, workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request, workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13, AT-24 |
| EC-19 | PROPOSED | workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request, workflow:proposals/ci/verify-claimed-shas.yml#pr-fields, procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | PROPOSED | script:scripts/ci/verify_claimed_shas.py, workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request, workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09, AT-18 |
| EC-21 | PROPOSED | script:scripts/ci/check_declared_base.py, workflow:.github/workflows/governance-acceptance.yml#declared-base~events=pull_request, workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-22 | ACTIVE | script:scripts/ci/status_vocabulary_lint.py, workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request | AT-04, AT-11, AT-19 |
| EC-23 | PROPOSED | script:scripts/ci/secret_baseline_check.py, script:scripts/ci/secret_control_check.py, workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request, workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks, workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks-all-refs | AT-15, AT-21, AT-22 |
| EC-24 | PROPOSED | procedural:docs/evidence/push-log/PUSH_LOG.md, procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md, script:scripts/ci/evidence_sequence_check.py, workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request | AT-13, AT-20, AT-24 |
| EC-25 | ACTIVE | script:scripts/ci/workflow_deploy_trigger_lint.py, workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request, procedural:proposals/remediation/README.md | AT-03, AT-14 |
| EC-26 | PROPOSED | script:scripts/ci/reading_note_map_check.py, workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request, workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08, AT-17, AT-23 |
<!-- END RULE INVENTORY -->

Status vocabulary: `ACTIVE` = in force now (pinned constitutions / standing decisions, and
EC-17/EC-22/EC-25 per auditor ruling D-5 of 2026-10-20, in force upon the D-2 merge `b851e01`);
`PROPOSED` = amendment not yet in force (EC-16, EC-18..EC-21, EC-23, EC-24, EC-26 — activation
workflows must run green on feature branches first, Direction 6); no other tokens are valid in
this table. Round-19 liveness rule (EC-26 fix): for `ACTIVE` rows, a `workflow:` reference must
be unpinned, under `.github/workflows/`, and present in the working tree — enforced by
`scripts/ci/reading_note_map_check.py` (AT-17) and run live by the governance-acceptance
workflow. Round-20 trigger-liveness rule (auditor direction 5): an ACTIVE row's `workflow:`
reference must also declare the events its enforcement claims via the `~events=` suffix
(e.g. `~events=push+pull_request`), and the live workflow's `on:` block must actually contain
each declared event with branch filters covering the enforcement branch (`main`);
manual-only, branch-disabled, or otherwise untriggered workflows FAIL (AT-23). PROPOSED rows
cite the inert `proposals/ci/` drafts (activation path) and, where recorded, the live jobs that
already run their checks on this branch — running a check is not activating a rule. A rule being in
force (policy) is recorded separately from whether its automated
checks run (enforcement); the live jobs are named above per row. Enforcement grammar: `script:PATH` / `procedural:PATH` (PATH must exist);
`workflow:FILE#JOB` (JOB must exist in FILE; optional `~events=E1+E2` declares the required
triggers — REQUIRED and validated for ACTIVE rows); `pytest:NODE` (the file part of NODE must exist).
Acceptance grammar: `AT-nn` (must exist in `scripts/ci/tests/ACCEPTANCE_TESTS.json`),
`pytest:NODE`, `workflow:FILE#JOB`, or `NONE` (a declared absence — honest, never hidden).
The required set for the reading-note check = all ACTIVE rows (plus PROPOSED rows when the
checker runs with `--include-proposed`, as this proposal branch does for self-testing).
