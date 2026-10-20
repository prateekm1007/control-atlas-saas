# P0_STATUS — Implementation Status, Blockers, Test Evidence, Candidate SHAs

**Status date:** 2026-10-20 (Round 13)
**Authority:** OD-01 target layout; maintained per the Engineering
Constitution (EC §6, §12, §14). Decision status lives in
`docs/decisions/OWNER_DECISIONS.md` (canonical copy on the
`docs/owner-decisions` branch); this file tracks implementation only.

## 1. Round-13 candidate branches (exact SHAs)

Prepared and tested this round per the Round-13 directions. All are
feature branches for review — **none is merged**; the hard stop (OD-16)
holds.

| Branch | Base | Tip SHA | Contents |
|---|---|---|---|
| `hygiene/h1-remove-gemini` | `main` @ `2195a048` | `a761d722…` | H-1: provider credential literal and every provider reference removed from `main`'s tree (21 files); deterministic offline narrative compiler replaces the provider-backed one; acceptance grep clean |
| `governance/canonical-docs` | `main` @ `2195a048` | (this commit) | OD-01: root `AGENTS.md`, root `ENGINEERING_CONSTITUTION.md` (byte-exact blob `b6bb35e8…`), `docs/coder/BRIEF.md`, `docs/P0_STATUS.md` (this file) |
| `docs/owner-decisions` | `main` @ `2195a048` | `d1d766df…` | Register rev 5 at `docs/decisions/OWNER_DECISIONS.md`: OD-01…OD-21 AUDITOR-DECIDED with decision status separate from implementation status; OD-07 per the Round-13 SHA-256 ruling; A-01 accepted risk |
| `feature/od07-audit-id` | `science/p0-correctness-2026-10-19` @ `ebc1fd79` | `73ef84c1…` | OD-07: deterministic full-SHA-256 audit IDs (canonical input + canon hash + matrix hash); `structure_digest` informational; legacy 8-char lookup; OD-08 key scope on lookup; `test_ids.py` asserts the ruling |
| `feature/n8-altlocs-per-residue` | `feature/od07-audit-id` @ `73ef84c1` | `af2c1641…` | N8 / R-4: `per_residue` altloc policy as the ENGINE DEFAULT (per_atom preserved as opt-in); pinned-cohort verdict diff stated (empty) |
| `feature/provenance-v2` | `science/p0-correctness-2026-10-19` @ `ebc1fd79` | `4bddfe7c…` (+ evidence `6ce5ee85…`) | Runner provenance v2: package-ownership audit, installed-artifact route, explicit unavailability, strict fail-closed, `calibration_executed=false` record marker |
| `feature/modrs-census-v2` | `science/p0-correctness-2026-10-19` @ `ebc1fd79` | `dea6826b…` | MODRES census v2: repository-relative input artifact (digest-recorded), canonical-source coordinate fetch with per-target digest binding; Round-11 result reproduced |

Round-12's four local commits (`0c7ea77`, `03d20a1`, `46c79ea`,
`3eb62b5`) were never pushed and were lost in the sandbox reset; the
Round-13 series above is new work, not a restoration.

## 2. Workstream status (auditor's Round-13 checklist: about 54%)

| Workstream | Weight | Done | Contribution |
|---|---|---|---|
| Security and hardening | 20 | 75% | 15.0 |
| Engine correctness | 25 | 65% | 16.25 |
| AlphaFold acquisition | 15 | 70% | 10.5 |
| Calibration validity | 25 | 25% | 6.25 |
| Governance and decisions | 5 | 60% | 3.0 |
| Hygiene | 5 | 60% | 3.0 |
| Deployment readiness | 5 | 0% | 0.0 |

Work counts at full value only where its commit exists on GitHub;
the branches in §1 are candidates pending push/review/merge.

## 3. Test evidence (local, this environment)

Environment: Python 3.12.15, gemmi 0.7.5, pytest 9.x; provider SDK
absent. The 4 `TestUsageWidget` failures are environmental (streamlit
not installed locally); they pass in CI, which installs
`backend/requirements-dev.txt`.

| Tree | Command | Result |
|---|---|---|
| pristine `main` @ `2195a048` | `cd backend && python3 -m pytest tests -q` | 4 failed, 249 passed, 5 skipped, 1 error (pre-existing: TestUsageWidget ×4; `test_b1_tokens` module-level `sys.exit` on import) |
| pristine science @ `ebc1fd79` | same | 4 failed, 502 passed, 5 skipped |
| `feature/od07-audit-id` | same | 4 failed, 513 passed, 5 skipped (−7 old +18 new `test_ids`; zero regressions) |
| `feature/n8-altlocs-per-residue` | same | 4 failed, 526 passed, 5 skipped (+13 per-residue tests; zero regressions) |
| `feature/provenance-v2` | same | 4 failed, 515 passed, 5 skipped (+13 provenance tests; zero regressions) |

Hash locks verified bit-identical on every science-based branch (canon
`6a9cd4b4349b81de`, matrix `ed1731ace09a9a3b…`, law-semantic vs
`governance/law_semantic.sha256`).

Pinned-cohort verdict diff (N8 default flip, direction 6): 40/40
experimental targets audited under both altloc policies from canonical
RCSB bytes — **0 verdict changes, 0 deterministic-score changes, 0
failing-law changes**; 2/40 targets (1EJG, 3NIR) retain 3 fewer atoms
under `per_residue` (chimera prevention; both PASS/83 either way).
Evidence: `docs/audit-evidence/2026-10-20/n8-altlocs/` on the N8 branch.

## 4. CI state (GitHub, pre-Round-13)

- Stability Gate triggers on push to **any** branch (science-branch
  workflow: gitleaks secret scan, lock+dev installs, unit tests, three
  hash pins). `main`'s legacy workflow lacks the secret scan and
  dev-dependency install; `main` CI was already red before Round-13
  (run 22597810035, failed ~22s, logs expired) — an inherited baseline,
  not caused by Round-13 work.
- Science tip `ebc1fd79` CI: **success** (run 39253683033).
- Round-13 branch pushes will trigger CI on their exact SHAs; results
  are to be reported as-is once the pushes execute.

## 5. Blockers and open owner actions

- **Hard stop (OD-16):** no merges to `main`, no deploys, no v25 run, no
  held-out evaluation, no threshold/cohort changes, no golden-artifact
  regeneration.
- **A-02:** provider key rotation at the issuer — tree removal done
  (`hygiene/h1-remove-gemini`); issuer-side revocation unverified.
- **A-03:** EC2 production secrets (`TOSCANINI_API_KEYS`,
  `INTERNAL_CALLBACK_SECRET`, `ALLOWED_ORIGINS`) — owner-executed.
- **A-04:** branch protection required-check configuration on `main` —
  owner-executed (rulesets endpoint returned none to the auditor).
- **A-05:** governance merge — patch prepared (`governance/canonical-docs`);
  merge is owner-gated.
- **Push credential:** the Round-11 PAT was observed rejected
  (`401 Bad credentials`) on 2026-10-20; Round-13 pushes await a valid
  token (A-01 remains an accepted risk per the register; rotation
  recommended, not required).
- Not started (owner-gated by decision or hard stop): OD-03/OD-04/OD-05
  law-semantic changes; OD-18 corrected-cohort encoding; OD-19 strata;
  OD-20 NMR policy; C3–C11.

## 6. Constitution reading

Both constitutions were read in full from their pinned blobs before any
edit this round (domain `211ad2c5…`, engineering `b6bb35e8…`); the
reading note citing every EC/OBS/CRIT ID enforced is committed to every
Round-13 branch at
`docs/audit-evidence/2026-10-20/constitution-reading-note-2026-10-20.md`.
