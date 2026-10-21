# P0_STATUS — Implementation Status, Blockers, Test Evidence, Candidate SHAs

**Status date:** 2026-10-20 (Round 21; earlier-round sections 1–8 below are retained
as historical evidence — the per-round state is cumulative, and each section
carries its own round label; the CURRENT state is sections 9–10 plus the
Round-21 evidence table)
**Authority:** OD-01 target layout; maintained per the Engineering
Constitution (EC §6, §12, §14). Decision status lives in
`docs/decisions/OWNER_DECISIONS.md` (canonical copy on the
`docs/owner-decisions` branch); this file tracks implementation only.

## 1. Round-14 branch state (exact SHAs, all verified on GitHub)

| Branch | Base | Tip SHA | Stability Gate |
|---|---|---|---|
| `hygiene/h1-remove-gemini` | science @ `ebc1fd79` (+ merge of the parallel `66ffcbcb`) | `51494350e43eac85a00aaba6e84c88cc3d500b46` | **SUCCESS** — run 39376518935 |
| `feature/od07-audit-id` | `08a9403` (pushed) + compliance fix | `c02f7742000be14747af0a0015cd18c82625831f` | **SUCCESS** — run 39376279686 |
| `feature/n8-altlocs-per-residue` | `08a9403` (parallel writer's commits `eb7e565`+`b5afd6c`) | `b5afd6ca01a345432852d5865a79999e464636b4` | **SUCCESS** — run 39373864542 |
| `feature/provenance-v2` | science @ `ebc1fd79` | `6ce5ee8503ea40c760239d8b79687966b90afd5c` | **SUCCESS** — run 39376283537 |
| `feature/modrs-census-v2` | science @ `ebc1fd79` | `4328156d443aa1a8bb3804b34c057b1a31b6d0db` | first run failed (gitleaks false-positive on the digest artifact); fixed by a scoped `.gitleaksignore`; latest run: see report |
| `governance/canonical-docs` | main @ `2195a048` | (this commit) | see report |
| `docs/owner-decisions` | main @ `2195a048` (+ merge of the parallel `0dcc5a1`) | (register rev 6 commit) | see report |

Round-13's local-only SHAs (`a761d72`, `73ef84c`, `af2c164`, `a632f59`,
`d1d66d`) were superseded by the Round-14 rebuilds (direction 4); they
remain preserved in the Round-13 durability bundles and patch series
with recorded SHA-256s.

## 2. Incident (open, owner review requested)

A parallel writer using this sandbox's git identity (`Z User
<z@container>`) pushed `feature/od07-audit-id` @ `08a9403` at ~05:07 UTC
and, during this session, `feature/n8-altlocs-per-residue`,
`evidence/round11-fix-2026-10-20`, `hygiene/h1-remove-gemini` @
`66ffcbcb` (main-based, non-compliant with the H-1 ruling), and
`docs/owner-decisions` @ `0dcc5a1` (a divergent "rev 6, 45 rows"
register). Full record: `evidence/incident-push-2/` on the H-1 and
OD-07 branches. No force-push; no deletion; their branches are resolved
by fast-forward merges that keep the compliant trees at the tips. My
Round-13 "PAT dead" claim was my own error (empty credential from a
non-propagated shell variable) and is corrected there. Token-level
attribution of the parallel pushes requires the owner's review (A-09).

## 3. Workstream status (auditor's Round-14 checklist: about 52%)

| Workstream | Weight | Done | Contribution |
|---|---|---|---|
| Security and hardening | 20 | 75% | 15.0 |
| Engine correctness | 25 | 60% | 15.0 |
| AlphaFold acquisition | 15 | 70% | 10.5 |
| Calibration validity | 25 | 25% | 6.25 |
| Governance and decisions | 5 | 50% | 2.5 |
| Hygiene | 5 | 60% | 3.0 |
| Deployment readiness | 5 | 0% | 0.0 |

## 4. Test evidence (local, this environment: Python 3.12.15, gemmi 0.7.5)

| Tree | Result (full backend suite) |
|---|---|
| pristine science @ `ebc1fd79` | 4 env-fail / 502 passed / 5 skipped |
| H-1 P0 stack (`6d6600e` tree) | 4 env-fail / 502 passed / 5 skipped (identical; zero regressions) |
| OD-07 fix (`c02f774`) | 4 env-fail / 523 passed / 5 skipped (their 08a9403 tests + ruling tests; zero regressions) |
| provenance-v2 / census-v2 / N8 | per their evidence files; hash locks bit-identical everywhere |

The 4 failures are the environmental TestUsageWidget cases (streamlit
absent locally; green in CI, which installs `requirements-dev.txt`).

## 5. Blockers and open owner actions

- **Hard stop (OD-16):** no merges to `main`, no deploys, no v25 run, no
  held-out evaluation, no threshold/cohort changes, no golden-artifact
  regeneration.
- **A-09:** owner review of the token history for the parallel pushes;
  the owner provided the repository PAT in-chat (Round 14) and directed
  its use (recorded; no revocation reminders per the owner's
  instruction).
- **A-02:** provider key issuer-side revocation still unverified (tree
  removal is complete on the H-1 branch).
- **A-03:** EC2 secrets (`TOSCANINI_API_KEYS`,
  `INTERNAL_CALLBACK_SECRET`, `ALLOWED_ORIGINS`); the former provider
  key is no longer required (H-1).
- **A-04:** branch protection required-check configuration on `main`.
- **A-05:** governance merge — patch ready on this branch; owner-gated.
- **Register divergence:** the parallel writer's "rev 6, 45 rows"
  register on the remote `docs/owner-decisions` history is superseded by
  this round's canonical rev 6 (Round-14 directions); the divergence is
  recorded in both.
- Not started (owner-gated): OD-03/04/05 law-semantic changes;
  OD-18 corrected-cohort encoding; OD-19 strata; OD-20 NMR policy;
  C3–C11.

## 6. Constitution reading

Both constitutions were re-read from their pinned blobs before any edit
this round (domain `211ad2c5…`, engineering `b6bb35e8…`); the Round-14
reading note citing every EC/OBS/CRIT ID enforced is committed to the
H-1, OD-07, governance, and register branches at
`docs/audit-evidence/2026-10-20/constitution-reading-note-round14.md`.

## 7. Amendment-package correction cycle (added 2026-10-20, SERIES-2026-10-20-EC-FIX)

Branch `proposal/ec16-ec26-enforcement` (base `main` @ `2195a048…`, tip advanced from
`17eb881…` by the 3.1–3.9 correction commits — exact tips in the push log and round report;
all commits appended, no force-push):

| Item | State on the branch tip |
|---|---|
| EC-16…EC-26 + R-AS/R-DL adjudications | Recorded in register **rev 7** (`docs/owner-decisions` @ `4f546160…`); amendments remain PROPOSED / inert |
| 3.1 EC-26 | Rule inventory derived from `docs/directives/INDEX.md`; fail-closed reading-note checker; canonical reading note (real read window recorded); PR-body input removed |
| 3.2 EC-20 | Named acceptance tests EXECUTED via pytest on the exact candidate (evidence JSON: argv, SHA, counts, exit) |
| 3.3 EC-21 | Remote base verification (exact returned-SHA match; fail-closed error taxonomy), line-unique Stack-parent validation, shallow fail-closed |
| 3.4 EC-22 | Machine-checkable `STATUS:` grammar validated against a runs manifest built from Actions records; progress computed, not asserted |
| 3.5 EC-16 | Strict per-entry SHA lists (any malformed entry rejects the list); 404/422 = absent vs 401/403/429/5xx = VERIFICATION_FAILURE; least-privilege workflows |
| 3.6 EC-18/24 | Parent-prefix append-only check; two-phase push log (INTENT/CONFIRM); truthful P-011-CONFIRM with remote-verified SHA + timestamp reconciliation |
| 3.7 EC-25 | Inert dispatch-only deploy remediation + gate-trigger fix under `proposals/remediation/`; lint activation preconditioned on the remediation; `[ain]` claim backed by exact path/blob/line evidence and de-misattributed |
| 3.8 EC-23 | Action pinned to API-verified commit `ff98106e…`; binary checksum-pinned (official checksums.txt); all-refs job (`--log-opts=--all`); seeded positive + clean negative controls |
| 3.9 | Counts corrected (SEVEN workflow drafts / six checkers — receipt R8/R9, append-only); acceptance suite AT-01…AT-16 (91 assertions) all passing locally; CI = NOT RUN (no run IDs exist; no CI-pass claim) |

**Test evidence (local, exact commands in `docs/evidence/ec16-ec26/LOCAL_TEST_LOG.md` on the
proposal branch):** `python3 scripts/ci/tests/run_all.py` → 16 groups / 91 `check()`
assertions / 0 failed (Python 3.12.15, pytest 9.0.2). GitHub Actions for the correction
commits: **NOT RUN** — the amendment workflows are inert drafts under `proposals/ci/`; no
check-run claim is made for any commit in this cycle.

**Unchanged program-level state:** the 40-row cohort metadata table remains OUTSTANDING;
r2 archive acceptance remains BLOCKED pending independent verifier execution against the
delivered bytes; the hard stop (OD-16) is in force; no merge to `main`, no deployment, no
v25 calibration, no held-out evaluation, no threshold/cohort changes, no golden-artifact
regeneration occurred in this cycle.

## 8. Round-19 mainline repair (added 2026-10-20, SERIES-2026-10-20-R19)

Branch `repair/round19-mainline-acceptance`, created from the exact current `main`
@ `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` (auditor direction 2: repair forward through a
reviewed, narrowly scoped PR; no direct push to main; no force-push; no administrative merge).

| Direction | Item | State on the repair branch |
|---|---|---|
| 3A | Two remaining non-documentation provider-SDK references | REMOVED (`proposals/remediation/toscanini-gate.branch-triggers.yml`, `fixtures/toscanini_gate_AIN_DEFECT.yml`); reviewed content carried from the unmerged `c4a0fcc` edits (content-only; no history merge; no push-log copy) |
| 3B | Remediation README | CORRECTED to the current live state (gate `["**"]` triggers, dispatch-only inert deploy); patches labeled SUPERSEDED/HISTORICAL with exact historical blob links; policy vs active enforcement kept distinct |
| 3C | Package manifest | Scope defined and machine-enforced (`scripts/ci/manifest_check.py`); `SHA256SUMS` regenerated over the full 38-file scope; clean-checkout verification recorded in the round report |
| 4 | Push-log collision + force-update deviation | Recorded in `docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md` (dated, evidence-backed; no historical entry rewritten); `P-017-INTENT` appended with a globally unused ID (census across all published branch logs); validator hardened (duplicate IDs, malformed confirmations, missing fields, INTENT/CONFIRM pairing with PENDING disclosure) |
| 5 | EC-20 job failures (PRs #18/#19/#20) | Root cause found and fixed: repo-wide `pytest --collect-only` aborts with INTERNALERROR when importing legacy script-style modules that `sys.exit()` at import time (`gpu_worker/tests/*.py`), so named acceptance tests registered as "not collectable" (runs 39417485854/39417488957/39417492355; reproduced locally on the exact head `748544f` in a pytest-only environment before fixing). Collection is now scoped per-file (AT-18) |
| 6 | EC-26 checker gap | FIXED: ACTIVE rules must map to LIVE workflows (unpinned, under `.github/workflows/`, present in the working tree); negative tests added (AT-17). Rerun on the consolidated activation candidate follows after this repair merges |
| 7 | Live end-to-end enforcement | `.github/workflows/governance-acceptance.yml` added (gate untouched): acceptance suite AT-01..AT-20, reading-note/inventory validation with ACTIVE-liveness, status vocabulary against an in-job runs manifest, push-log schema + parent-prefix, deploy-trigger lint on the actual tree, package manifest + forbidden-reference checks, PR-body claim verification with executed acceptance tests, declared-base verification, full-history all-refs secret scan with controls and a documented baseline classification (known H-1 history fingerprints disclosed every run; A-02 stays open) |

**Status of the activation program after this repair (direction 6 adjudications):** EC-21
accepted in principle (rebase onto the repaired base, records reconciled, all checks rerun on
the exact final candidate — then register update + reviewed-PR merge; NOT yet active);
EC-26/EC-16/EC-19/EC-20 fixes landed here, activation candidates to be built/rerun (NOT yet
active); EC-18/EC-23/EC-24 candidates to be built and tested; EC-17/EC-22/EC-25 remain in
force with live enforcement from this repair forward. PRs #18/#19/#20 (three copies of the
same workflow) are superseded by the consolidated candidate once it passes.

**Unchanged program-level state:** A-01 remains ACCEPTED RISK (no rotation reminder, per
direction 8); the scientific hard stop is preserved (no law semantics, calibration,
thresholds, cohorts, held-out evaluation, golden artifacts, or deployment); the portable
40-row cohort metadata table and the independent r2 archive acceptance remain OPEN and are
not closed by this repair.

## 9. Round-20 acceptance-defect repair (added 2026-10-20, SERIES-2026-10-20-R20)

Auditor ruling: Round-19 is NOT closed. The green checks on PR #30's head
`f44bf4d…` (gate `39443515998`, Governance Acceptance `39443516042`) did not
detect two critical deficiencies in the new controls; a green base gate is
necessary, not sufficient. This section records the repair, implemented on the
same branch `repair/round19-mainline-acceptance` (base `main` @ `f535a71…`,
fix-forward from `f44bf4d…`; no merge approved; no direct push to main; no
force-push).

| Direction | Defect found by the audit | Fix on this candidate |
|---|---|---|
| 3 | EC-23 all-refs scan tolerated every generic-api-key finding under docs/ or evidence/ (a new credential there would pass); positive control accepted any nonzero scanner exit | Broad category exemption REPLACED by the reviewed fingerprint allowlist `docs/evidence/ec16-ec26/secret_scan_baseline.json` (47 documented entries; source/classification/reason each; unknown fingerprints FAIL wherever they live; H-1 separate under A-02) enforced by `scripts/ci/secret_baseline_check.py`; controls verified by `scripts/ci/secret_control_check.py` (expected finding + fired rule `aws-access-token` + key redaction; `--exit-code 101` disambiguates findings from execution errors, which fail closed; negative control verified against its report). AT-21, AT-22 |
| 4 | EC-24 checker verified SHA syntax but never the live remote ref, nor that evidence supports the SHA; the R19 follow-up pushes were unreconciled | `evidence_sequence_check.py --remote-verify`: live-ref comparison per CONFIRM/RECONCILE (CURRENT/HISTORICAL), evidence support (SHA + branch in evidence; contradicting cited run IDs FAIL), branch-under-test tip declaration (unpaired INTENT = PENDING disclosed; undeclared tip FAILS); the eight R19 tips reconciled truthfully by RECONCILE records P-050..P-057 (never backdated; schema v3); the `push-log-integrity` CI job runs the full verification with event-context branch/tip. AT-24 |
| 5 | EC-26 liveness checked path + job but not the `on:` triggers (manual-only would pass) | ACTIVE workflow references must declare `~events=…` and the live workflow must actually trigger on each declared event with branch coverage of `main`; manual-only / branch-disabled / untriggered FAIL. AT-23 (AT-17 retained) |
| 6 | Register rev 8 stranded outside the reviewed path; PR #30 carried rev 7 | Rev-8 content carried into the reviewed path content-only (no history merge of `docs/owner-decisions-r19`) and extended as register **rev 9** (Part 12: the Round-20 records with the decision/policy/implementation/candidate/CI/live-on-main distinctions). This P0_STATUS section and the directives INDEX are updated from verified results; the Round-19 status report remains as historical evidence |
| 7 | Activation adjudications | All eight rules remain PROPOSED per the Round-20 table (register Part 12.6); PR #30 merges only on exact-candidate green + the required independent review; #31-#36 serially revalidated after that merge; #23-#29 close only on proven replacement coverage |

**Status of this repair:** the fixes, tests, and records above are complete
on the candidate; CI runs on the exact final candidate are the completing
evidence (run IDs recorded in the push-log confirmations and the round
report when observed). The two program-level standing items are unchanged
and NOT closed by this work: the independent r2 archive verification and the
portable 40-row cohort metadata table remain OPEN; A-01 remains ACCEPTED
RISK (no rotation reminder, per direction 8); the scientific hard stop is
preserved in full.

## 10. Round-21 acceptance-defect repair (added 2026-10-20, SERIES-2026-10-20-R21)

Auditor ruling: **PR #30 is not approved for merge yet.** Round 20's green CI
was real but necessary-not-sufficient; the independent review found three
critical defects in the new controls plus the unproven protected-main merge
transition. This section records the Round-21 repair, implemented on the same
branch `repair/round19-mainline-acceptance` (fix-forward from `62fab8a…`;
no merge approved; no direct push to main; no force-push).

| Direction | Defect found by the audit | Fix on this candidate |
|---|---|---|
| 3 | The EC-23 control checker could print the report's `Secret` field value in its failure path (disclosing the very secret it detected); a requested capture log was only checked when present; the HEAD-tree clean path trusted scanner exit status | `secret_control_check.py` rewritten: NO path prints `Secret`/`Match` values or raw report/log bytes (safe metadata + booleans only); `--log-file` requested-but-missing/unreadable FAILS closed; the HEAD-tree and clean-repo negative controls now require the report to exist, parse, and contain zero findings (exit status alone is never trusted); the workflow routes the HEAD-tree result through the hardened negative mode with its capture log. AT-25 (10 checks incl. stdout/stderr key-absence proofs on every failure path) |
| 4 | EC-24 accepted any unpaired intent as declaring the branch tip (an old intent could cover later commits); CONFIRMs paired by branch/FIFO rather than identifier; `expected_remote_tip_before_push` was unvalidated | Identifier pairing (`P-058-CONFIRM` consumes exactly `P-058-INTENT`); the single accepted terminal intent must be INTRODUCED BY THE EXACT CANDIDATE (a later commit without a new intent FAILS), declare `expected_remote_tip_before_push` equal to the newest recorded branch state (an intervening ref update or stale declaration FAILS; ABSENT requires no prior record), and the checked remote tip must be the exact expected tip; excess unmatched intents FAIL unless voided by the new `intent-reconcile` verb; terminal state labeled PENDING, never COMPLETE. AT-26 (8 checks) |
| 5 | `HISTORICAL` only proved the SHA was not the live tip — not that the push was real; runs absent from the manifest were tolerated as unverified; failed runs were not distinguished | Every HISTORICAL confirmation must be proven: ancestry from the live tip, or a reviewed `evidence-back` record (new verb) whose SHA is a real remote commit (object store or the commits API) with run citations verified; every run cited anywhere in an entry (field or body) must verify via the manifest, else the Actions API (`--repo` + token) — an unverifiable citation FAILS the gate (a manifest omission is never an approval); failed/cancelled runs are disclosed and cannot prove a push; run `head_branch` connects a SHA to its branch. The four off-lineage legacy confirmations (P-011/P-012/P-013/P-014) are backed by reviewed records P-062..P-065 (commits-API existence + push-event run connections, honestly dated). AT-27 (8 checks) |
| 6 | The post-merge `main` push had no demonstrated declaration path | New `merge-intent` verb (P-061 pre-declares the future reviewed merge of PR #30: base `f535a71…`, head DEFERRED-resolved to the exact final candidate). On the source branch the active merge-intent must bind the exact candidate and the declared base must be the live target tip; on `main`, a push is DECLARED only by a VERIFIED merge-intent: the candidate must be a merge commit whose parents are the declared base + approved head, the intent must ride the PR branch history, and the PR must be merged with an approving non-author review (GitHub API in CI; fixtures in tests). Raced/unrelated/direct/unapproved updates FAIL. Deterministic AT-28 (8 checks: legitimate merge, unapproved, author-approval, unmerged PR, raced main, direct push, source binding, stale intent) |
| 7 | EC-26 checked triggers and branch filters but not path filters or filter types | ACTIVE workflow references must declare `~governs=` (the changes they enforce); `paths`/`paths-ignore` filters that exclude a governed path FAIL; `~governs=**` (tree-wide) cannot be satisfied by a restrictive paths filter; branch/path filters of invalid TYPE (dict/int) FAIL; unknown annotation keys FAIL. All ACTIVE rows in the INDEX inventory now carry `~governs=`. AT-29 (9 checks) |
| 8 | Status/evidence misalignment: P0_STATUS top metadata referred to Round 14; INDEX R19/R20 wording ambiguous; register lacked the Round-21 adjudication; no durable evidence table | This section + the top metadata correction; INDEX R19 row now reads CLOSE-OUT REJECTED / R20 SUPERSEDED-BY-R21 / new SERIES-2026-10-20-R21 row (all explicit that nothing is live on `main`); register **rev 10 Part 13** appended (decision status distinct from implementation, exact candidate, CI status, main-activation status); the Round-21 evidence table is GENERATED FROM LIVE EVIDENCE (ls-remote + the Actions API at generation time; generator committed, output committed as `docs/evidence/ec16-ec26/R21_EVIDENCE_TABLE.md`); schema v4; reading note replaced for this cycle; SHA256SUMS regenerated. The PAT-hygiene sweep is recorded in the coder report as a LOCAL claim only — GitHub CI does not independently verify it, and no token value appears in any file, worklog, or the evidence table |
| 9 | Merge gating | PR #30 remains UNAPPROVED: the full gate must pass on the exact final candidate, the post-merge `main` transition must be demonstrated (AT-28 + the pre-declared P-061), and the required independent review must be recorded; merge only through the protected PR mechanism with a MERGE COMMIT. After that merge: PRs #31–#36 revalidate serially against the resulting `main` with fresh exact-head results, their declared bases updated, and their own merge-intents declared before merging; EC-21 is not waived; #23–#29 close only on proven unique-content coverage |

**Acceptance suite:** AT-01..AT-29 (29 pytest groups; `run_all.py` and
`ACCEPTANCE_TESTS.json` updated). Local result on the exact working tree
before commit: 29/29 passed (~38s, Python 3.12.15, pytest 9.x). No result is
claimed for CI runs that have not yet executed on the final candidate SHA.

**Unchanged program-level state:** the scientific hard stop is preserved in
full; A-01 remains ACCEPTED RISK (no rotation reminder, per the Round-21
PAT instruction); A-02 issuer-side revocation, the independent r2 archive
verification, the portable 40-row cohort metadata table, and the EC2
deployment secrets remain OPEN and are not closed by this work.

## 11. Round-21-followup suite-soundness repair (added 2026-10-20, SERIES-2026-10-20-R21B)

Independent verification of the Round-21 candidate `0fe4840…` (pristine
checkout, pytest run with the Round-19 soundness probe) found that the CI
"Acceptance suite (AT-01..AT-29)" green carried **15 test nodes with hidden
check() failures**: the R20/R21 rewrites of `test_acceptance.py` lost the
conftest soundness fixture, so a `check()` failure recorded inside a test node
never failed that node (the f44bf4d class recurring; EC §6/§13). The hidden
nodes were AT-04/09/10/11/12/14/16/17/18/21/23 (R20 classes) and
AT-24/26/27/28 (R21 classes).

| Defect (hidden failure) | Genuine repair on this candidate |
|---|---|
| Inline JSON passed to `--runs-manifest`/`--api-json` (checkers correctly fail closed on an unusable fixture file) | File-path fixtures (AT-04, AT-10, AT-12) |
| EC-20 evidence bodies missing the fail-closed `Claimed-SHAs: none`/`Directive-ID` fields | Fields added (AT-09, AT-18) |
| AT-11 fixture claimed pct=50 while the records computed 66% | RUN_FAIL record added: 2 accepted / 4 statuses = 50% |
| AT-10 stack-mismatch scenario could not fail (every stack class mapped to the same ref); unrelated-base scenario collided (identical commit timestamps produced the SAME SHA in two synth repos) | Orphan science line; salted `_synth_repo`; proj2 fetches proj3 for the locally-present-but-unrelated case |
| AT-14: the draft's `on:` block still carried the `ain]` defect (blob 3e9e0b4d restored by the R20 rename); the AIN fixture was normalized (did not reproduce the defect); `main]`-substring collision in the body check | Draft repaired at block level (pull_request `branches: [main]`, push `["**"]`); true-defect fixture restored byte-exact; precise-token regex `(?<!m)ain\]` |
| AT-16 expected 6 checkers and a RECEIPT.md that R20/R21 had removed | 11 checkers (R21 adds generate_evidence_table.py); RECEIPT.md completed append-only (R13/R14/R15 correction entries) |
| AT-17/AT-23: ACTIVE citation lacked `~events=`; PROPOSED case declared an event the draft lacks | `~events=push~governs=**` citation; the PROPOSED-without-annotation case plus a NEW explicit negative check (declared-but-absent event on a PROPOSED row fails) |
| AT-21: baseline entry contradicted its own fingerprint (file=docs/x.md vs fingerprint's docs/digests.json) | Entry fields DERIVED from the fingerprint |
| AT-24: single-letter branch names matched inside hex SHAs and the body fallback made the branch check vacuous; citation cases used a stale tip | `topic/b` branch names; branch must be named IN the evidence field (checker hardening); citation cases pinned to the live tip cC |
| AT-26: voided-intent case reused intents introduced by an earlier commit (the R21 binding rule correctly rejects that) | Fresh intent IDs introduced by the same candidate commit |
| AT-27: off-lineage scenarios left the log commit's tip undeclared; the ghost-back case did not bind | Terminal intents (DEFERRED binding) for the log commits; the ghost back now BINDS to its entry so the object-store proof is what rejects it |
| AT-28: the direct-push scenario dropped the log file (checker crashed with an unhandled traceback instead of failing closed) | The pushed log is carried; the checker also fails closed with a clear message on a missing log |
| Checker: the terminal intent's declaration violations were an elif chain (only the first fired) | Independent findings, with a documented-rebuild waiver (a binding evidence-back for the same branch+SHA makes a non-ancestral pre-push tip legitimate) |

**Result:** `29 passed` under pytest WITH the conftest fixture active, three
consecutive runs on the exact candidate, plus `run_all.py` green; every
previously-hidden failure genuinely repaired; no check weakened or deleted
(one explicit negative check ADDED at AT-23). P0_STATUS's own prose carries no
machine-checkable status claims beyond the run records cited above.

**Unchanged program-level state:** the scientific hard stop is preserved in
full; A-01 remains ACCEPTED RISK (no rotation reminder, per direction 8); the
independent r2 archive verification, the portable 40-row cohort metadata table,
A-02 issuer-side revocation, and EC2 deployment secrets remain OPEN. PR #30 is
still NOT approved for merge: the merge gate remains (all checks green on the
exact final candidate + the demonstrated post-merge transition + the required
independent approval through the protected path).

## 12. Round-22 acceptance-defect repair (added 2026-10-20, SERIES-2026-10-20-R22)

**Auditor ruling (binding):** Round-21 is NOT CLOSED and PR #30 is NOT APPROVED
FOR MERGE. Green CI is necessary but not sufficient: the independent audit
found further fail-open cases in EC-23, EC-24 and EC-26. `main` stays
untouched at `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`; work proceeds on
`repair/round19-mainline-acceptance` through ordinary fast-forward pushes and
the reviewed PR path only.

**Preflight (auditor direction 1, this cycle):** both constitutions re-read in
full from the exact pinned blobs (`constitution.md`
`211ad2c5d4b40059d9f3f03de12e0701f9235bc7`, `ENGINEERING_CONSTITUTION.md`
`b6bb35e84956d84a8c9d9e65734c212baf3ed075`) with the real read window
`2026-10-20T19:06:52Z`–`19:07:45Z`; the register (rev 10), P0_STATUS (§1–11),
INDEX, PUSH_LOG_SCHEMA (v4), the full push log through P-073-INTENT, the
reconciliation and incident records, the evidence table and RECEIPT read in
full before the first edit (~19:25Z). A parallel-writer observation is
disclosed: two fast-forward bookkeeping pushes (R21B-3 `47022be`, R21B-4
`908de6d`, by the same sandbox identity as the documented incident class)
landed during the preflight window carrying push-log records only; this cycle
was rebased forward onto `908de6d` with no conflict (they touched no checker,
workflow, or documentation content), and the single-writer pre-check is
re-verified before this cycle's own push.

**Repairs implemented in this cycle (auditor directions 2–7):**

- **EC-23 (direction 2):** `scripts/ci/secret_control_check.py` no longer
  interpolates ANY report-derived value into diagnostics (the Round-22 audit
  found `RuleID`/`File`/`StartLine` leak paths in the failure diagnostics and
  negative-control output). Every report element is type-validated before
  examination (non-object rows, non-string RuleID/File, non-int StartLine,
  non-string Secret/Match fail closed with FIXED messages and row contents
  withheld); unexpected findings report a COUNT only; the raw-key checks cover
  the report bytes even when the schema validation fails; unexpected internal
  errors fail closed with a fixed message. AT-30 proves the assembled
  synthetic key never reaches stdout/stderr on every path (RuleID, File,
  string StartLine, non-object rows, malformed JSON, non-array reports,
  missing reports/logs, negative paths, scanner-error masquerade).
- **EC-24 (direction 3):** `scripts/ci/evidence_sequence_check.py` now
  separates commit existence, branch-push provenance, and CI success as
  INDEPENDENT claims. The runs manifest and the Actions-API verifier retain
  and validate each cited run's EVENT TYPE and HEAD_BRANCH (a push-event run
  on another branch contradicts the entry and FAILs). An off-lineage
  HISTORICAL confirmation must prove the branch push by (a) a push-event run,
  (b) a live Events-API PushEvent, (c) a recorded `push_event=` citation
  adjudicated by a machine-verifiable `review_ref=` into the canonical
  auditor-adjudicated record, or (d) an explicit `push-not-proven`
  reclassification (valid `review_ref=`; disclosed, never counted as
  evidence). A commits-API 200 or a free-text rebuild narrative is NOT push
  provenance. P-063/P-065 prove their pushes by their push-event runs (their
  failure conclusions are disclosed separately); P-062/P-064 prove only
  commit existence and are re-adjudicated in register Part 14 (events-API
  search results recorded there). AT-31 covers: never-proven-pushed commits,
  wrong-branch runs, non-push events, missing API proof, rebuild narratives,
  valid/invalid/duplicate reclassifications, live and recorded PushEvents,
  and the PushEvent before/after binding of the declared pre-push tip (false
  ABSENT, prior unlogged ref updates, raced/absent events).
- **EC-24 (direction 4):** intent reconciliations must target PRECEDING
  eligible intents, match branch scope (push-intent branch; merge-intent
  source_branch — unambiguous semantics, no free-form values), void exactly
  once, and never void consumed intents; a CONFIRM naming a voided intent
  FAILS. The unmatched-intent invariant is GLOBAL: no branch (including
  branches other than the one under test) may hold more than one unvoided
  unmatched intent. AT-32 covers all of these plus the merge-intent
  source-branch scoping.
- **Merge transition (direction 5):** the workflow gains the narrowly scoped
  `pull-requests: read` permission and a live API canary (authenticated
  PR/review/events reads must succeed, else the job fails closed BEFORE the
  checker runs). The post-merge verifier now checks the PR API's
  `merge_commit_sha` against the protected branch's live candidate, the
  source `head.ref`, `merged_at`, and a FRESH non-author approval applied to
  the EXACT head (commit_id binding, submitted before the merge, latest
  state per reviewer — dismissed/superseded approvals rejected), plus the
  merge PushEvent's `before` equaling the declared base. A new
  `pre-merge-verification` job (PR events only) checks the CURRENT
  PR/head/review state BEFORE any merge executes and is EXPECTED TO FAIL
  while no independent approval exists on the exact head — that is the
  functioning pre-merge gate. DEFERRED merge-intent supersession is resolved
  via introducing commits on both the source and target sides. AT-33 covers
  all failure modes plus a second sequential merge (the #30 → #31 pattern).
  A-04 (exact branch-protection required checks) remains an OPEN blocker
  until independently verified.
- **EC-26 (direction 6):** `scripts/ci/reading_note_map_check.py` replaces
  the representative-sample coverage test with a SOUND conservative
  glob-coverage check (a structural witness set generated from each governed
  glob: multiple name shapes per wildcard component, `**` at depths 1–3,
  root/hidden/workflow paths for whole-tree governance); ordered negative
  patterns follow GitHub Actions filter semantics (`!` negations, last match
  wins) for paths AND branches; empty filter lists and malformed event
  configuration types (string/list/int event values) fail closed. AT-34
  covers the auditor's full list (docs/** vs docs/X, `['**','!docs/evidence/**']`,
  `branches: ['**','!main']`, empty lists, malformed configurations, both
  pattern orders, witness-shaped paths-ignore).

**Acceptance suite:** AT-01..AT-34 registered and green on this candidate
(pytest + `run_all.py`, 34/34). The registry, INDEX inventory rows, and the
governance-acceptance job names are updated; no existing check was weakened
or deleted (AT-27/AT-28 fixtures were strengthened to the new required PR
fields: merge_commit_sha, merged_at, commit_id, submitted_at).

**Push-log protocol (direction 8):** this cycle's push carries P-073-CONFIRM
for the `908de6d` push (with its verified run evidence), the Round-22
re-adjudication records for the legacy off-lineage confirmations (per the
events-API search results recorded in register Part 14), a fresh push intent
bound to the exact introducing candidate with the true pre-push tip
`908de6d…`, the P-072 supersession (void + redeclaration bound to the exact
final candidate and the live base), and the regenerated manifest and evidence
table. The terminal intent stays PENDING (disclosed, never COMPLETE); no push
is manufactured to eliminate it.

**Outstanding (not closed by this cycle):** PR #30 still lacks an independent
current-head approval (the pre-merge-verification job will remain red until
it is recorded); the merge awaits all Round-22 blockers fixed on the exact
final head, the full Stability Gate + Governance Acceptance green there, the
passing pre-merge approval verification, and the protected-path MERGE COMMIT;
after any verified merge, PRs #31–#36 revalidate serially against the new
`main` (EC-21 not waived) and #23–#29 close only on proven unique-content
coverage. A-04 branch-protection verification, the independent r2 archive
verification, the portable 40-row cohort metadata table, A-02 issuer-side
revocation, and EC2 deployment secrets remain OPEN. The scientific hard stop
is preserved in full; A-01 remains ACCEPTED RISK (no rotation reminder, no
credential operations, token values never written anywhere).

## 13. Round-24 governance recovery (added 2026-10-21, SERIES-2026-10-21-R24) — OPEN

- **PR #30 was merged into `main` without a qualifying independent approval.** Merge commit
  `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (first parent `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`,
  second parent `d81c874eaedc13c6fdbaed5102ad7e638449036d`), merged 2026-10-20T21:55:06Z by
  `prateekm1007`. The only review recorded on PR #30 is a `COMMENTED` review by
  `chatgpt-codex-connector[bot]`. The merge is recorded as an unapproved merge; it is NOT
  retroactively approved and the Round-23 justification in register Part 14 is withdrawn.
- **The required `pre-merge-verification` job is missing on `main`.** `main` carries workflow
  blob `939261fbbfb200f0f4e83553b7840642f826aaea`, in which the job survives only as a comment
  recording its removal. The Round-23 deletion is **superseded and non-compliant**. The gate is
  restored on `repair/r23-governance-recovery` byte-identically to the verified pre-R23 blob
  `44fe653ec4e526725eaea4266dc1807222d3be06` (commit `443d6063b42c40ae47b6ea0607e45d31c12239d2`),
  keeping its `pull-requests: read` permission, authenticated API canary, exact-head review
  verification and fail-closed behaviour. PR #30's fixes are kept; no blind revert is performed.
- **The governance runs failed and stay visible as failures.** Governance Acceptance run
  `39504181997` on `main` @ `63da0e3…` — FAIL (no qualifying non-author approval on the exact
  PR #30 head). Governance Acceptance run `39507046623` on the repair-branch candidate — FAIL
  (`P-081-INTENT` has no CONFIRM; its declared pre-merge base `f535a71…` is historical while
  `main` is at `63da0e3…`). Stability Gate runs `39504182072` and `39507046937` concluded
  `success` on the same candidates — necessary, not sufficient.
- **The pruned refs and the incidentally closed PRs are restored.** All 25 branch refs exist at
  their exact original heads (created through the git-data API, each verified by a ref read and
  an independent `ls-remote`; GitHub recorded `head_ref_restored` for PR #30's head at
  2026-10-20T22:17:45Z). PRs #18–#29 and #31–#36 are open again with titles, bodies, heads and
  review history unchanged. Inventory: 36 PRs — 4 merged (#1, #15, #17, #30) and 32 open; 55
  branch heads; 26 tag refs.
- **Preservation evidence is produced by automation.** `scripts/ci/generate_preservation_bundle.py`
  builds the machine-readable bundle (refs, pull requests, unique commits, branch diff matrix,
  preservation matrix, runtime dependency map, PR evidence, commands log, SHA256SUMS) and
  `.github/workflows/preservation-bundle.yml` publishes it as a downloadable CI artifact; its URI
  and digest are recorded in the recovery PR body and its checksums are verified from the
  artifact itself. The repository's own `SHA256SUMS` file is not evidence of bundle completeness.
- **Recovery remains OPEN.** The recovery PR is blocked until a legitimate, non-author approval
  exists on its exact head and the required checks are green there; the restored
  `pre-merge-verification` job is expected to fail until then — that is the gate functioning.
  The push-log reconciliation (P-085 `push-reconcile`, P-086 `intent-reconcile`) becomes live
  only when that PR merges through the protected path.
- **OD-16 scientific/production hard stop: IN FORCE.** No deployment, no v25 run, no held-out
  evaluation, no threshold or cohort change, no law-semantic change, and no golden-artifact
  modification is authorized by this recovery. A-01 remains ACCEPTED RISK (no rotation reminder,
  no rotation initiation, no credential operations; token values are never written anywhere).
