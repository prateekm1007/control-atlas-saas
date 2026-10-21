# Constitution reading note — 2026-10-21, Round-24 governance recovery (series SERIES-2026-10-21-R24)

Recovery branch `repair/r23-governance-recovery`, created from the exact current `main` tip
`63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (the merge commit of PR #30). No history was
rewritten, no ref on `main` was moved, and no prior push-log entry was edited or deleted by
this series. The Round-23 authorization (gate removal + merge without a qualifying independent
approval) is **superseded and non-compliant**; this note records the recovery, not an
endorsement of those actions.

## 1. Sources read (verified before reading; Round-24 recovery session)

| Document | Pinned blob (full SHA-1) | Verified | Source commit at read time | Read window (UTC) |
|---|---|---|---|---|
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | YES — blob ID re-hashed from the base tree (`git ls-tree 63da0e3…`) | `main` @ `63da0e3…` | 2026-10-20T22:0xZ–22:1xZ, read in full (75 lines) |
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES — same method, identical ID | same tree | read in full (11 lines; OBS-153 and CRIT-154 both ACTIVE) |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES — same method, identical ID | same tree | read in full (178 lines, EC §1–§15 + Final Rule) |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | YES | same tree | read in full |
| `docs/decisions/OWNER_DECISIONS.md` (register rev 10, Parts 1–14) | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` | YES | same tree | read (Part 14 carries the Round-22 adjudications; Part 15 appended by this series) |
| `docs/P0_STATUS.md` (sections 1–12) | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` | YES | same tree | read in full; section 13 appended by this series |
| `docs/directives/INDEX.md` | `244f76eab550951a35305350c76af6ac4ed19b92` | YES | same tree | read in full (series register + rule inventory parsed; 20 ACTIVE rules derived); edited by this series (R23/R24 series rows, superseded directives, live-status corrections, the restored gate as an EC-§7 enforcement point) |
| `docs/evidence/reading-note/READING_NOTE_LATEST.md` (Round-23 note, superseded by this note) | `b8468317dc6b23beb4406ec23569fd12d81b4823` | YES | same tree | read in full; superseded — its "sanctioned" characterization of the Round-23 actions is corrected here |
| `docs/evidence/push-log/PUSH_LOG.md` (entries P-001…P-081-INTENT) | `45a924456eddfac857cdcd5aec4f6c8c8ee56c9b` | YES | same tree | read in full, including P-078-CONFIRM, P-079, P-080 and P-081-INTENT |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` (schema v5) | `83aa798b7a9731e63be2a9664e502b7631a42674` | YES | same tree | read in full (record phases, reconciliation verbs, remote verification) |
| `.github/workflows/governance-acceptance.yml` (on `main`, gate removed) | `939261fbbfb200f0f4e83553b7840642f826aaea` | YES | same tree | read in full — the `pre-merge-verification` text survives only as a comment saying the job was removed; it is **not** an executable job |
| `.github/workflows/governance-acceptance.yml` (verified pre-R23 version, commit `443d6063b42c40ae47b6ea0607e45d31c12239d2`) | `44fe653ec4e526725eaea4266dc1807222d3be06` | YES | `443d6063…` | read in full (585 lines, including the `pre-merge-verification` job definition at lines 352–458) |
| `.github/workflows/toscanini-gate.yml` | `ff3d66d919666fe234220ed922d636b8d775f249` | YES | same tree | read in full (2 jobs; push + pull_request on `**`; unchanged by this series) |
| `scripts/ci/evidence_sequence_check.py` | `1ae66acf131b59b9b05f8cd06ad9262e8eab74bb` | YES | same tree | key sections read in full (schema, parent-prefix, remote verification, merge transition, `_valid_approvals`); **not edited** by this series |
| `scripts/ci/reading_note_map_check.py` | `467b97a86c1edc88da4e0399efc3a19ef8434004` | YES | same tree | contract read; not edited |
| `scripts/ci/status_vocabulary_lint.py` | `5b407803f728388069104e8cabfe1face84ffe25` | YES | same tree | contract read; not edited |
| `scripts/ci/manifest_check.py` | `b2c0b5a44e9dc03f18aebc88f7374b4a4146c144` | YES | same tree | scope policy read (8 scope globs; `.github/workflows/`, `docs/P0_STATUS.md` and `docs/decisions/` are deliberately out of scope); not edited |
| `scripts/ci/tests/test_acceptance.py` | `615b796ff03d8f97d2d8d49a237c855ad3e2c202` | YES | same tree | AT-33 read in full (the pre-merge gate fixture mechanism); AT-35/AT-36 appended by this series |
| `scripts/ci/tests/ACCEPTANCE_TESTS.json` | `a37b8bc6963aa17b02cff70e1c8ac9262eab9193` | YES | same tree | read; AT-35/AT-36 registered by this series |
| `scripts/ci/generate_preservation_bundle.py` | (new in this series) | n/a — authored here | `repair/r23-governance-recovery` | new file; generates the preservation bundle; not a modification of any prior checker |
| `.github/workflows/preservation-bundle.yml` | (new in this series) | n/a — authored here | `repair/r23-governance-recovery` | new file; publishes the bundle as a downloadable CI artifact; read-only with respect to the repository |
| `scripts/ci/tests/conftest.py`, `scripts/ci/tests/run_all.py` | `4869c1d5f47ef390766d67f25c9b789ea83c06c0`, `17ce9a1c45a0b40bcebb2454c17e0bdc70881402` | YES | same tree | read (suite soundness fixture) |

## 2. Rule → enforcement → acceptance mapping (EC-26; the inventory in `docs/directives/INDEX.md` is normative)

### 2a. ACTIVE rules (in force; live enforcement points after this series)

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | Purpose and authority; fail closed on ambiguity | procedural:AGENTS.md | NONE |
| EC-§2 | Pre-coding protocol (constitutions first) | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request~governs=docs/evidence/reading-note/**+docs/directives/INDEX.md+scripts/ci/reading_note_map_check.py; procedural:AGENTS.md | AT-06, AT-08, AT-17, AT-23, AT-29, AT-34 |
| EC-§3 | Security is a release requirement (secrets, fail closed; logs MUST not expose credentials) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request~governs=**; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request~governs=**; pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request~governs=** |
| EC-§4 | Product and scientific integrity | procedural:AGENTS.md; pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock~events=push+pull_request~governs=** |
| EC-§5 | Architecture and maintainability | procedural:AGENTS.md | NONE |
| EC-§6 | Testing and evidence (deterministic, exact-SHA results) | script:scripts/ci/tests/run_all.py; workflow:.github/workflows/governance-acceptance.yml#acceptance-suite~events=push+pull_request~governs=scripts/ci/**; procedural:AGENTS.md | AT-16 |
| EC-§7 | CI/CD and supply chain (exact-candidate checks; pinned tooling; PR review before merge) | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request~governs=.github/workflows/**+scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification~events=pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py; procedural:.github/workflows/toscanini-gate.yml; procedural:branch protection on main | AT-03, AT-14, AT-33, AT-35, AT-36 |
| EC-§8 | Reliability and operations | procedural:AGENTS.md | NONE |
| EC-§9 | Customer experience and cost control | procedural:AGENTS.md | NONE |
| EC-§10 | LLM and external-provider boundaries | procedural:AGENTS.md | NONE |
| EC-§11 | Change, release, incident discipline | procedural:AGENTS.md; procedural:docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md; procedural:docs/evidence/push-log/RECONCILIATION-2026-10-20-R19-FOLLOWUP-PUSHES.md; procedural:docs/evidence/push-log/INCIDENT-2026-10-21-R23-GOVERNANCE-BYPASS.md | NONE |
| EC-§12 | Definition of done ("not run" never "passed") | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request~governs=docs/evidence/**+scripts/ci/status_vocabulary_lint.py; procedural:AGENTS.md | AT-11, AT-19 |
| EC-§13 | Prohibited shortcuts (no gate weakening; no secrets in logs, reports, or artifacts) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request~governs=** | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks)~events=push+pull_request~governs=** |
| EC-§14 | Required task artifacts | procedural:AGENTS.md; procedural:SHA256SUMS (package manifest, machine-checked) | NONE |
| EC-§15 | External engineering references | procedural:ENGINEERING_CONSTITUTION.md | NONE |
| OBS-153 | Avidity quantification (contact density) — ACTIVE domain rule | procedural:constitution.md | NONE |
| CRIT-154 | Structural admissibility (warhead-RMSD, steric veto) — ACTIVE domain rule | procedural:constitution.md | NONE |
| EC-17 | Push on completion (durable storage before subsequent work) | procedural:docs/evidence/push-log/PUSH_LOG.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28, AT-31, AT-32, AT-33, AT-36 |
| EC-22 | Status vocabulary (RUN_PASS needs a real run ID) | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary~events=push+pull_request~governs=docs/evidence/**+scripts/ci/status_vocabulary_lint.py | AT-04, AT-11, AT-19 |
| EC-25 | No ungated deploy on push | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint~events=push+pull_request~governs=.github/workflows/**+scripts/ci/workflow_deploy_trigger_lint.py; procedural:proposals/remediation/README.md | AT-03, AT-14 |

### 2b. PROPOSED rules (not in force; unchanged by this series)

| Rule ID | Enforcement point(s) | Acceptance test(s) |
|---|---|---|
| EC-16 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#claimed-shas | AT-01, AT-12 |
| EC-18 | script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request; workflow:proposals/ci/evidence-sequence-check.yml#sequence | AT-05, AT-13, AT-24 |
| EC-19 | workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#pr-fields; procedural:docs/directives/INDEX.md | AT-01, AT-07 |
| EC-20 | script:scripts/ci/verify_claimed_shas.py; workflow:.github/workflows/governance-acceptance.yml#claim-verification~events=pull_request; workflow:proposals/ci/verify-claimed-shas.yml#acceptance-tests-executed | AT-09, AT-18 |
| EC-21 | script:scripts/ci/check_declared_base.py; workflow:.github/workflows/governance-acceptance.yml#declared-base~events=pull_request; workflow:proposals/ci/check-declared-base.yml#declared-base | AT-02, AT-10 |
| EC-23 | script:scripts/ci/secret_baseline_check.py; script:scripts/ci/secret_control_check.py; workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs~events=push+pull_request~governs=**; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks; workflow:proposals/ci/gitleaks-all-branches.yml#gitleaks-all-refs | AT-15, AT-21, AT-22, AT-25, AT-30 |
| EC-24 | procedural:docs/evidence/push-log/PUSH_LOG.md; procedural:docs/evidence/push-log/PUSH_LOG_SCHEMA.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity~events=push+pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#pre-merge-verification~events=pull_request~governs=docs/evidence/push-log/**+scripts/ci/evidence_sequence_check.py | AT-13, AT-20, AT-24, AT-26, AT-27, AT-28, AT-31, AT-32, AT-33, AT-35, AT-36 |
| EC-26 | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory~events=push+pull_request~governs=docs/evidence/reading-note/**+docs/directives/INDEX.md+scripts/ci/reading_note_map_check.py; workflow:proposals/ci/reading-note-map-check.yml#map | AT-06, AT-08, AT-17, AT-23, AT-29, AT-34 |

## 3. What the recovery session found (live evidence; nothing below is inferred)

1. **The Round-23 merge was real and is on `main`.** PR #30 was merged as merge commit
   `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (first parent `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`,
   second parent `d81c874eaedc13c6fdbaed5102ad7e638449036d`), `merged_at`
   2026-10-20T21:55:06Z, `merged_by` `prateekm1007`. The PR carries only a bot `COMMENTED`
   review: **no qualifying approving review from a non-author on the exact head**. The
   post-merge Governance Acceptance run `39504181997` failed on exactly this point
   ("PR 30 has no approving review from a user other than the PR author that applies to the
   exact declared head and predates the merge … review states: chatgpt-codex-connector[bot]:
   latest state 'COMMENTED'"); the Stability Gate run `39504182072` on the same commit
   concluded `success`. The failure is preserved, not reframed.
2. **The gate removal is non-compliant and is repaired forward, not reverted.** On `main`
   the `pre-merge-verification` job exists only as a comment (blob `939261fb…`). This series
   restores the workflow file byte-identically to the verified pre-R23 blob
   `44fe653ec4e526725eaea4266dc1807222d3be06` (commit `443d6063…`), so the job is executable
   again with its `pull-requests: read` permission, its authenticated API canary, its
   fail-closed behaviour and its exact-head review verification. PR #30's fixes are kept; no
   blind revert was performed.
3. **The two governance failures are different failures and are recorded separately.** The
   main-run failure above is the missing independent approval. The later repair-branch
   Governance Acceptance run `39507046623` failed for the push-log reasons: `P-081-INTENT`
   has no CONFIRM, and its declared pre-merge base `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`
   is no longer the live `main` tip (`63da0e3…`). The Stability Gate run `39507046937` on the
   same candidate concluded `success`. This series reconciles the push-log state with the
   schema's own mechanisms (P-085 `push-reconcile`, P-086 `intent-reconcile`); it does not
   manufacture a retroactive CONFIRM and does not treat the historical base as current.
4. **The 25 pruned branch refs are restored at their exact original heads, and the 18
   incidentally closed PRs are open again.** Each ref was created through the git-data API at
   the exact recorded SHA and verified twice (git ref read + independent `ls-remote`);
   GitHub recorded `head_ref_restored` for PR #30's head at 2026-10-20T22:17:45Z. The
   restoration of `activation/ec18-r19` used the live PR #34 head
   `97067ad0320968bfb73f7e045408fcb60ccd5c7d`; a previously transcribed variant of that SHA
   (`…fcb0c8a20e5a`) does not correspond to any live ref or PR head and was **not** used.
   PRs #18–#29 and #31–#36 were reopened with their titles, bodies, heads and review history
   unchanged; PR #30 (merged) was not touched.
5. **The recovery PR is expected to stay red until a genuine approval exists.** The restored
   `pre-merge-verification` job fails while the exact candidate head carries no non-author
   approving review. That failure is the gate functioning, and the PR must remain blocked.

## 4. Commands executed (recovery session; credentials redacted)

- `git ls-remote` / `GET /repos/{repo}/git/ref/heads/{branch}` — pre- and post-restoration
  census of every branch head (55 heads, 26 tag refs, 36 pull-head refs, 14 pull-merge refs).
- `POST /repos/{repo}/git/refs` with `{"ref": "refs/heads/<branch>", "sha": "<recorded sha>"}`
  — 25 creations, HTTP 201 each; verified by a second independent `ls-remote`.
- `GET /repos/{repo}/pulls?state=all` and `GET /repos/{repo}/pulls/{n}` — PR inventory; the
  list endpoint omits the `merged` flag, so every PR's state was re-read from its detail
  endpoint.
- `PATCH /repos/{repo}/pulls/{n}` with `{"state": "open"}` — 18 reopenings, each re-read
  afterwards to confirm `state=open` and unchanged head/title/body.
- `GET /repos/{repo}/branches/main/protection` — effective protection configuration
  (required checks, required review count, strictness, administrator enforcement).
- `git show 443d6063:.github/workflows/governance-acceptance.yml` — the verified pre-R23 gate
  job; `git hash-object` on the restored file returns `44fe653ec4e526725eaea4266dc1807222d3be06`.
- `python3 scripts/ci/generate_preservation_bundle.py` (new; `scripts/ci/generate_preservation_bundle.py`)
  — the preservation bundle, uploaded by the new `.github/workflows/preservation-bundle.yml`
  workflow as a downloadable CI artifact; its URI and digest are recorded in the recovery PR
  body (the artifact, not the repository's `SHA256SUMS` file, is the bundle's checksum
  evidence).

## 5. Open items and stops (this series does not close the recovery)

- The recovery PR needs a **legitimate, non-author approval on its exact head** and
  green required checks before it may merge. Until then it stays blocked; no self-approval,
  retroactive approval or weakened test is acceptable.
- `docs/evidence/push-log/PUSH_LOG.md` on `main` still carries the un-reconciled `P-081`
  state until this series merges; the reconciliation records (P-085, P-086) are authored here
  and become live only through the reviewed merge.
- The OD-16 scientific/production hard stop remains in force: no deployment, no v25 run, no
  held-out evaluation, no threshold or cohort change, no golden-artifact modification is
  authorized by this recovery.
