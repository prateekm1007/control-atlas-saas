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
| SERIES-2026-10-20-R18 | Auditor Round-18 directions 1–7 (2026-10-20): fix-forward from `f535a71`, activation fixes EC-16/19/20 + EC-18/21/23/24/26, incident + collision records, manifest completion | 2026-10-20 | ACTIVE — all seven exact candidates GREEN; awaiting auditor activation rulings + required review (no merge authorized by this cycle) | 1–6 (delivered; 7 = stop honored) | Direction 1: Round-17 treated as historical evidence; live refs re-verified (`main @ f535a71`, PRs #18–#22 open, activation branches at their reported SHAs). Direction 2: constitution preflight re-executed fresh (pinned blobs `211ad2c5…` / `b6bb35e8…` + AGENTS.md + 5 governance files; recorded with real read windows). Direction 3: fix-forward branch from exactly `f535a71` (provider-SDK refs removed, remediation docs corrected, manifest completed + checker added, acceptance gate added to active CI). Direction 4: PR #17 review-bypass recorded as INC-1; P-015-CONFIRM collision recorded as INC-2; `--pushlog-schema` extended (pairing + remote verification; pending ≠ complete). Direction 5: EC-20 collection root cause proven (`gpu_worker/tests/test_credits.py` raises `SystemExit` at import → pytest INTERNALERROR; run 39417492355); checker collection scoped to the acceptance suite; duplicate activation PRs #18–#20 superseded by one consolidated candidate once green. Direction 6: EC-18/21/23/24/26 fixes rebuilt on new branches (no force-push); auditor adjudications recorded in the register. Direction 7: PAT A-01 unchanged; hard stop maintained. **Round outcomes (verified, exact-SHA CI):** PR #23 fix-forward green (acceptance gate + ordinary gate, runs 39443533044/39443544992); PR #24 CONSOLIDATED EC-16/19/20 candidate green (all three activation jobs + gates, runs 39446141415/39446076597) — PRs #18/#19/#20 marked superseded per direction 5; PR #25 EC-18 rebuild green (run 39437335099); PR #26 EC-23 rebuild green (runs 39441644268/39441653042 — all three jobs, incl. the verified synthetic positive control and the classified 46-fingerprint disposition); PR #27 EC-24 rebuild green (run 39438198824 — pairing + remote verification + PENDING semantics live); PR #28 EC-26 rebuild green (run 39438290430 — computed Enforcement column validated); PR #29 EC-21 rebuild green (run 39443598543) — remains PROPOSED per the auditor ruling. Disclosed process defects this round (all caught by the fail-closed gates and repaired append-only): two skipped CONFIRM entries, one manifest bytecode-scope defect, one parser dot-stripping defect, hand-typed SHA tails in P-022/P-023/P-024/P-025/P-026/P-027/P-029/P-032 (corrected via P-035/P-036/P-037 records), two self-referential disposition-doc matches, one structural FP ('commits API:' + public id), and three disclosed no-force-push-protected rebuilds of the never-merged consolidated activation candidate. |

## In-force rule set (EC-26)

Per the auditor's 2026-10-20 ruling (D-5, under full delegation): **EC-17 (push on completion), EC-22 (status vocabulary), and EC-25 (no deploy on push) are IN FORCE** upon merge of D-2 (`governance/canonical-docs` merged to main as `b851e01`, gate run 39412416212 = success). They are recorded `ACTIVE` in the rule inventory below; the auditor records this activation in the decision register.

**EC-16, EC-18, EC-19, EC-20, EC-21, EC-23, EC-24, and EC-26 remain PROPOSED**: each becomes in force only when its activation workflow has run green on a feature branch AND the auditor records it in the register (Direction 6 of the same ruling).

The base in-force set is unchanged: the EC sections enforced by the pinned constitutions (`ENGINEERING_CONSTITUTION.md` blob `b6bb35e8…`, `constitution.md` blob `211ad2c5…`) plus OBS-153 / CRIT-154 — as mapped in the Round-12 reading note (worklog Task 35) and in `docs/evidence/reading-note/READING_NOTE_LATEST.md`.

## Superseded directives

None recorded in this index. Historical supersessions prior to 2026-10-20 (register revs 1/3/4 rows, D-cycle, etc.) are recorded in `docs/decisions/OWNER_DECISIONS.md` (register rev 6, branch `docs/owner-decisions`) and are not duplicated here.

## Rule inventory (EC-26, machine-readable)

<!-- BEGIN RULE INVENTORY (parsed by scripts/ci/reading_note_map_check.py --index; format below is normative) -->
reading-note-path: docs/evidence/reading-note/READING_NOTE_LATEST.md
pinned-blobs: b6bb35e84956d84a8c9d9e65734c212baf3ed075,211ad2c5d4b40059d9f3f03de12e0701f9235bc7
| Rule ID | Status | Enforcement | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|---|
| EC-§1 | ACTIVE | PROCEDURAL | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§2 | ACTIVE | LIVE | script:scripts/ci/reading_note_map_check.py, procedural:AGENTS.md@origin/governance/canonical-docs | AT-06, AT-08 |
| EC-§3 | ACTIVE | LIVE | workflow:.github/workflows/toscanini-gate.yml@origin/science/p0-correctness-2026-10-19#Secret Scan (gitleaks), pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock |
| EC-§4 | ACTIVE | PROCEDURAL | procedural:AGENTS.md@origin/governance/canonical-docs, pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock |
| EC-§5 | ACTIVE | PROCEDURAL | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§6 | ACTIVE | LIVE | script:scripts/ci/tests/run_all.py, procedural:AGENTS.md@origin/governance/canonical-docs | AT-16 |
| EC-§7 | ACTIVE | LIVE | script:scripts/ci/workflow_deploy_trigger_lint.py, procedural:.github/workflows/toscanini-gate.yml@origin/main | AT-03, AT-14 |
| EC-§8 | ACTIVE | PROCEDURAL | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§9 | ACTIVE | PROCEDURAL | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§10 | ACTIVE | PROCEDURAL | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§11 | ACTIVE | PROCEDURAL | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§12 | ACTIVE | LIVE | script:scripts/ci/status_vocabulary_lint.py, procedural:AGENTS.md@origin/governance/canonical-docs | AT-11 |
| EC-§13 | ACTIVE | LIVE | workflow:.github/workflows/toscanini-gate.yml@origin/science/p0-correctness-2026-10-19#Secret Scan (gitleaks) | workflow:.github/workflows/toscanini-gate.yml@origin/science/p0-correctness-2026-10-19#Secret Scan (gitleaks) |
| EC-§14 | ACTIVE | PROCEDURAL | procedural:AGENTS.md@origin/governance/canonical-docs | NONE |
| EC-§15 | ACTIVE | PROCEDURAL | procedural:ENGINEERING_CONSTITUTION.md@origin/policy/engineering-constitution | NONE |
| OBS-153 | ACTIVE | PROCEDURAL | procedural:constitution.md | NONE |
| CRIT-154 | ACTIVE | PROCEDURAL | procedural:constitution.md | NONE |
| EC-16 | PROPOSED | INERT | script:scripts/ci/verify_claimed_shas.py, workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-17 | ACTIVE | LIVE | procedural:docs/evidence/push-log/PUSH_LOG.md, script:scripts/ci/evidence_sequence_check.py | AT-13 |
| EC-18 | PROPOSED | LIVE | script:scripts/ci/evidence_sequence_check.py, workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13 |
| EC-19 | PROPOSED | INERT | workflow:proposals/ci/verify-claimed-shas.yml#pr-fields, procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | PROPOSED | INERT | script:scripts/ci/verify_claimed_shas.py, workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09 |
| EC-21 | PROPOSED | INERT | script:scripts/ci/check_declared_base.py, workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-22 | ACTIVE | LIVE | script:scripts/ci/status_vocabulary_lint.py, workflow:proposals/ci/status-vocabulary-lint.yml#vocabulary | AT-04, AT-11 |
| EC-23 | PROPOSED | INERT | workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks, workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks-all-refs | AT-15 |
| EC-24 | PROPOSED | PROCEDURAL | procedural:docs/evidence/push-log/PUSH_LOG.md, procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md | AT-13 |
| EC-25 | ACTIVE | LIVE | script:scripts/ci/workflow_deploy_trigger_lint.py, workflow:proposals/ci/no-deploy-on-push-lint.yml#lint, procedural:proposals/remediation/README.md | AT-03, AT-14 |
| EC-26 | PROPOSED | LIVE | script:scripts/ci/reading_note_map_check.py, workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08 |
<!-- END RULE INVENTORY -->

Status vocabulary: `ACTIVE` = in force now (pinned constitutions / standing decisions, and
EC-17/EC-22/EC-25 per auditor ruling D-5 of 2026-10-20, in force upon the D-2 merge `b851e01`);
`PROPOSED` = amendment not yet in force (EC-16, EC-18..EC-21, EC-23, EC-24, EC-26 — activation
workflows must run green on feature branches first, Direction 6); no other tokens are valid in
this table.

Enforcement vocabulary (Round-18, EC-26 revision — SEPARATE from policy status): `LIVE` =
verified running automation (a `.github/workflows/` file with a non-empty trigger and the
named job present, or a script invoked by such a live workflow — computed, not claimed);
`INERT` = the mapped automation exists only as drafts outside `.github/workflows/` (e.g.
`proposals/ci/*.yml`) and runs NOWHERE; `PROCEDURAL` = documented procedure and/or a checker
script that no live workflow invokes; `NONE` = nothing beyond the normative text. An `ACTIVE`
policy status NEVER implies `LIVE` enforcement — read the Enforcement column. The checker
recomputes every row and fails closed on a declared label that does not match the computed
class (scripts/ci/reading_note_map_check.py). Labels above were computed on the Round-18
fix-forward tree (`.github/workflows/acceptance-gate.yml` is what makes EC-§2/§6/§7/§12 and
EC-17/EC-22/EC-25's script points live).

Enforcement grammar: `script:PATH` / `procedural:PATH` (PATH must exist);
`workflow:FILE#JOB` (JOB must exist in FILE; LIVE additionally requires FILE under
`.github/workflows/`, a non-empty `on:` block, and the JOB present in the live tree);
`pytest:NODE` (the file part of NODE must exist).
Acceptance grammar: `AT-nn` (must exist in `scripts/ci/tests/ACCEPTANCE_TESTS.json`),
`pytest:NODE`, `workflow:FILE#JOB`, or `NONE` (a declared absence — honest, never hidden).
The required set for the reading-note check = all ACTIVE rows (plus PROPOSED rows when the
checker runs with `--include-proposed`, as this proposal branch does for self-testing).
