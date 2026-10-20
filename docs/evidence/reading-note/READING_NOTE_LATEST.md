# Constitution reading note — 2026-10-20, Round-19 mainline repair

- **Directive series:** SERIES-2026-10-20-R19 (auditor directions 1–8; see `docs/directives/INDEX.md`).
- **Discipline:** this read was performed BEFORE any repository edit in this session. The
  prior rounds' reading notes are NOT evidence for this cycle's preflight (auditor direction 1);
  both constitutions were re-read in full from the exact pinned blobs at the SHAs below, with
  `date -u` timestamps captured at read start. No retroactive evidence was manufactured.

## 1. Sources read (verified before reading; Round-19 cycle)

| Document | Pinned blob (full SHA-1) | Verified | Source commit at read time | Read window (UTC, real) |
|---|---|---|---|---|
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | YES — `git ls-tree` + `git cat-file` on the read commit returned the identical blob ID | `main` @ `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` (remote-verified by fresh ls-remote at session start) | start `2026-10-20T12:23:14Z` |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | YES — same method, identical ID; read in full (20572 bytes / 178 lines) | `main` @ `f535a71` | same window, immediately following |
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | YES | `main` @ `f535a71` | read in full (75 lines) |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | YES | `main` @ `f535a71` | read in full |
| `docs/decisions/OWNER_DECISIONS.md` (register rev 7) | `e4f7808d490411e257a523ff8a336515c45028ac` | YES | `main` @ `f535a71` | read in full (Parts 1–10) |
| `docs/P0_STATUS.md` | `40672ace565917ed06b77907ebee5d23fc460c22` | YES | `main` @ `f535a71` | read in full |
| `docs/directives/INDEX.md` | `040db2632ca76e3f3b556e9f18f5553ce4fb9cd7` | YES | `main` @ `f535a71` | read in full (inventory parsed) |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` | schema v2 (two-phase) | YES | `main` @ `f535a71` | read in full |

The first repository write of this cycle (the repair branch's first commit) happened after the
read window above. The repair branch was created from the EXACT current `main` SHA
`f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` (auditor direction 2; single-writer discipline).

## 2. Rule → enforcement → acceptance mapping (EC-26; the inventory in INDEX.md is normative)

Legend: `workflow:FILE#JOB` = a LIVE GitHub Actions job under `.github/workflows/` (Round-19
liveness rule: ACTIVE rules cite unpinned working-tree workflows only — a `proposals/ci/…`
draft is inert and can no longer satisfy an ACTIVE rule); `script:PATH` = checker script;
`procedural:PATH` = a tracked document stating the procedure. `AT-nn` = acceptance test in
`scripts/ci/tests/ACCEPTANCE_TESTS.json`. Policy status and enforcement status are recorded
separately: a rule in force is not proof its checks run.

### 2a. ACTIVE rules (in force; live enforcement points after this repair)

| Rule ID | Anchor | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-§1 | Purpose and authority; fail closed on ambiguity | procedural:AGENTS.md | NONE |
| EC-§2 | Pre-coding protocol (constitutions first) | script:scripts/ci/reading_note_map_check.py; workflow:.github/workflows/governance-acceptance.yml#reading-note-inventory; procedural:AGENTS.md | AT-06, AT-08, AT-17 |
| EC-§3 | Security is a release requirement | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks); workflow:.github/workflows/governance-acceptance.yml#secret-scan-all-refs; pytest:backend/tests/test_api_key_guard.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock |
| EC-§4 | Product and scientific integrity | procedural:AGENTS.md; pytest:backend/tests/test_certification.py | workflow:.github/workflows/toscanini-gate.yml#Unit Tests + Hash Lock |
| EC-§5 | Architecture and maintainability | procedural:AGENTS.md | NONE |
| EC-§6 | Testing and evidence | script:scripts/ci/tests/run_all.py; workflow:.github/workflows/governance-acceptance.yml#acceptance-suite; procedural:AGENTS.md | AT-16 |
| EC-§7 | CI/CD and supply chain (exact-candidate checks; pinned tooling) | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint; procedural:.github/workflows/toscanini-gate.yml | AT-03, AT-14 |
| EC-§8 | Reliability and operations | procedural:AGENTS.md | NONE |
| EC-§9 | Customer experience and cost control | procedural:AGENTS.md | NONE |
| EC-§10 | LLM and external-provider boundaries | procedural:AGENTS.md | NONE |
| EC-§11 | Change, release, incident discipline | procedural:AGENTS.md; procedural:docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md | NONE |
| EC-§12 | Definition of done ("not run" never "passed") | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary; procedural:AGENTS.md | AT-11, AT-19 |
| EC-§13 | Prohibited shortcuts (no gate weakening; no secrets) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks) | workflow:.github/workflows/toscanini-gate.yml#Secret Scan (gitleaks) |
| EC-§14 | Required task artifacts | procedural:AGENTS.md; procedural:SHA256SUMS (package manifest, machine-checked) | NONE |
| EC-§15 | External engineering references | procedural:ENGINEERING_CONSTITUTION.md | NONE |
| OBS-153 | Avidity quantification (contact density) — ACTIVE domain rule | procedural:constitution.md | NONE |
| CRIT-154 | Structural admissibility (warhead-RMSD, steric veto) — ACTIVE domain rule | procedural:constitution.md | NONE |
| EC-17 | Push on completion (durable storage before subsequent work) | procedural:docs/evidence/push-log/PUSH_LOG.md; script:scripts/ci/evidence_sequence_check.py; workflow:.github/workflows/governance-acceptance.yml#push-log-integrity | AT-13, AT-20 |
| EC-22 | Status vocabulary (RUN_PASS needs a real run ID) | script:scripts/ci/status_vocabulary_lint.py; workflow:.github/workflows/governance-acceptance.yml#status-vocabulary | AT-04, AT-11, AT-19 |
| EC-25 | No ungated deploy on push | script:scripts/ci/workflow_deploy_trigger_lint.py; workflow:.github/workflows/governance-acceptance.yml#ec25-trigger-lint; procedural:proposals/remediation/README.md | AT-03, AT-14 |

`NONE` rows are declared absences (honest, not hidden). Before this repair, EC-17/EC-22/EC-25
were in force as POLICY with no live automated enforcement (the drafts under `proposals/ci/`
are inert; that gap is exactly what the Round-19 live workflow closes).

### 2b. PROPOSED rules (not in force; activation pending per direction 6)

**Activation-record update on this candidate (activation/ec21-r19):** EC-21 flip to ACTIVE citing the live governance-acceptance jobs (see the INDEX inventory). This takes effect only when the complete gate is green on this exact candidate AND the auditor records the activation in the register (direction 6); until both happen the rules remain PROPOSED in force.

| Rule ID | Status on this candidate | Enforcement point(s) | Acceptance test(s) |
|---|---|---|---|
| EC-21 | ACTIVE on merge (record pending auditor registration) | ACTIVE | script:scripts/ci/check_declared_base.py, workflow:.github/workflows/governance-acceptance.yml#declared-base |


EC-16, EC-18, EC-19, EC-20, EC-23, EC-24, EC-26 remain PROPOSED. Their checks now run
live in the governance-acceptance workflow (claim-verification, push-log-integrity,
declared-base, secret-scan-all-refs, reading-note-inventory jobs), but RUNNING a check is not
ACTIVATING a rule: each becomes in force only when its activation candidate passes on the
exact final SHA and the auditor records it in the register (direction 6). The activation
states as adjudicated in Round-19: EC-21 accepted in principle (rebase + rerun required);
EC-26 blocked pending this checker correction (now landed; rerun required); EC-16/EC-19/EC-20
blocked on the EC-20 job failures (root cause found and fixed in this repair; rerun required);
EC-18/EC-23/EC-24 candidates to be built and tested. EC-17/EC-22/EC-25 remain in force under
the previous ruling, with live enforcement from this repair forward.

## 3. Round-19 context recorded at read time

- Current `main` = `f535a71` (verified by ls-remote at session start; the adjudicated baseline
  for repair-forward per direction 2). PR #17's merged tree carries the acceptance defects
  this repair addresses: two non-documentation provider-SDK references, a stale remediation
  README, an incomplete 23-entry manifest, and the P-015-CONFIRM collision recorded in
  `docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md`.
- The gate on `main` is green (run 39415713844) but ran no round-specific acceptance checks —
  green base gate is necessary, not sufficient (direction 6).
- Hard stop (OD-16) remains in force for everything outside this repair's authorized scope:
  no deploys, no v25 calibration, no held-out evaluation, no threshold/cohort changes, no
  golden-artifact regeneration. The 40-row cohort metadata table and the independent r2
  archive acceptance remain OPEN (direction 8).
