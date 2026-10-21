# Receipt — EC-16…EC-26 amendment proposal packaging (append-only, EC-18)

- Series: SERIES-2026-10-20-EC (see `docs/directives/INDEX.md`). Branch: `proposal/ec16-ec26-enforcement`, declared base `main` @ `2195a048e87dcec4aa63e8dc7981b3bb6d237ca0` (EC-21), declared stack parent: main-line.
- Grammar: `## <ISO-8601 UTC> | <ENTRY-ID> | <verb> | <subject> [| RECONSTRUCTED]`. Timestamps are `date -u` at the moment of the action, same host throughout. Entries added after their action are labeled RECONSTRUCTED.

## 2026-10-20T05:52:00Z | R1 | receive | auditor proposal message: EC-16..EC-26 amendments + rule-authoring standard + directive lifecycle + adoption order (EC-17+EC-25 first; EC-22+EC-23 next; rest follow); proposals pending owner approval | RECONSTRUCTED

- RECONSTRUCTED: session-start receipt time approximated from tool-call ordering; exact wall time not stamped at receipt. The message text is reproduced verbatim in the amendment document.

## 2026-10-20T05:53:30Z | R2 | verify | GET-only ls-remote: main=2195a048 (unchanged); three branches advanced beyond worklog Task 35 tips: feature/od07-audit-id 08a9403->c02f774, hygiene/h1-remove-gemini 66ffcbc->5149435, docs/owner-decisions 0dcc5a1->0773f0b; feature/n8-altlocs-per-residue=b5afd6c and evidence/round11-fix-2026-10-20=742e36b unchanged | RECONSTRUCTED

- RECONSTRUCTED: stamped from adjacent observed stamp; the ls-remote output itself is quoted in the worklog. Post-Task-35 commits carry git commit times 05:35:26Z-05:49:22Z (2026-10-20) and no worklog entry — worklog gap recorded in the INDEX.

## 2026-10-20T05:54:22Z | R3 | create | worktree + branch proposal/ec16-ec26-enforcement from declared base main@2195a048

## 2026-10-20T06:03:20Z | R4 | test | acceptance tests AT-01..AT-07 ALL PASS (24 checks, exit 0) after two fixture bugs fixed (AT-01 indented Directive-ID swallowed as continuation; AT-03 wrote path strings instead of YAML contents). First run 06:01:58Z: 21 pass / 3 fail — same root causes, fixed, re-run.

## 2026-10-20T06:03:29Z | R5 | record | EC-25 deploy-trigger lint against main's tree (report-only): .github/workflows/deploy.yml job 'deploy' triggers on push — EC-25 VIOLATION on main confirmed and recorded in docs/evidence/ec16-ec26/EC25_LINT_ON_MAIN.txt; NOT fixed (remediation goes through the approved merge path)

## 2026-10-20T06:04:16Z | R6 | write | package files authored: 3 amendment docs, directive INDEX, proposals/ci (README + 6 workflow drafts + PR template), scripts/ci (6 checkers + tests/run_all.py), EC-25 lint evidence; receipt and push log finalized in commit sequence

## 2026-10-20T06:05:30Z | R7 | verify | H-1 acceptance re-check on current tip 5149435: git grep for the three Gemini tokens (AIzaSy, GEMINI, generativeai) outside docs/audit-evidence/ + docs/security/ = ZERO hits (clean); OD-07 tip c02f774: compute_audit_id triple present in backend/main.py ingestion response | RECONSTRUCTED

- RECONSTRUCTED: action performed between R2 and R3 (before worktree creation, against fetched remote refs); stamped here when recorded, not at action time.

## 2026-10-20T09:52:30Z | R8 | record | correction of R6 inventory counts (auditor item 3.9): this branch carries SEVEN workflow YAML drafts under proposals/ci/ (check-declared-base, evidence-sequence-check, gitleaks-all-branches, no-deploy-on-push-lint, reading-note-map-check, status-vocabulary-lint, verify-claimed-shas) — not six as R6 stated — and SIX checker scripts under scripts/ci/; R6 is preserved unchanged (append-only)

- The R6 text ("6 workflow drafts", "6 checkers + tests/run_all.py") was accurate for checkers and wrong for workflows (the PR template and README were not workflows; verify-claimed-shas/check-declared-base/evidence-sequence-check/reading-note-map-check/status-vocabulary-lint/gitleaks-all-branches/no-deploy-on-push-lint are the seven). The audit checklist verified the same defect. AT-16 now enforces the corrected counts mechanically.

## 2026-10-20T09:52:31Z | R9 | test | acceptance suite EXTENDED and re-run locally: AT-01..AT-16, 91 check() assertions, all passing on the exact candidate (see docs/evidence/ec16-ec26/LOCAL_TEST_LOG.md for exact commands, SHA, environment, counts, exit codes); GitHub Actions CI status for these commits: NOT RUN — the workflows are inert drafts under proposals/ci/ and no run ID exists, so NO CI-pass claim is made or implied

- Supersedes R4's "24 checks" claim: the suite was extended per auditor item 3.9 to cover the 3.1–3.8 defect classes (fail-closed negatives, test execution, remote verification, manifest-backed statuses, strict SHA lists, append-only parent-prefix, remediation fixtures, action pinning, inventory counts). R4 remains preserved above, unchanged.

## 2026-10-20T13:45:00Z | R11 | receive | auditor Round-19 directions 1-8: repair forward from main @ f535a71; remove the two remaining provider-SDK references; correct the remediation README; regenerate the scoped manifest; record the P-015-CONFIRM collision and the a0e321c->c4a0fcc force-with-lease deviation; fix the EC-20 activation failures; fix the EC-26 liveness gap; wire the full acceptance battery into the live CI path; activation adjudications per rule

- Both constitutions re-read FIRST from the pinned blobs on `main` @ `f535a71` (read window and source SHAs in `docs/evidence/reading-note/READING_NOTE_LATEST.md`); earlier rounds' reading notes were not treated as this cycle's preflight evidence (direction 1).

## 2026-10-20T13:46:00Z | R12 | verify | preflight recon on the exact baseline: the two `google-generativeai` lines confirmed at `proposals/remediation/toscanini-gate.branch-triggers.yml:36` and `proposals/remediation/fixtures/toscanini_gate_AIN_DEFECT.yml:23`; `SHA256SUMS` confirmed at 23 entries vs the 35-file intended package; PR #17's merged tree verified as `f535a71`; EC-20 job failures reproduced locally on the exact head `748544f` in a pytest-only environment (collection aborts: INTERNALERROR importing `gpu_worker/tests/test_credits.py`, module-level `sys.exit(1)`)

- The reproduction evidence (checker output, evidence JSON with `tests_missing_from_collection`, collect probe showing the INTERNALERROR) is preserved in the round report; the fix (per-file scoped collection) was verified to EXECUTE the two named tests (counts {tests: 2, failures: 0}) on the same head before being committed.
## 2026-10-20T17:15:00Z | R13 | record | Round-19-fix correction cycle on this repair branch (auditor Round-19-fix directions 1-8): the count correction continues — the package now carries NINE checker scripts under scripts/ci/ (9 checker scripts: the original six plus ec23_scan_check.py plus the Round-19 additions already present) and SEVEN workflow drafts under proposals/ci/ (7 workflow drafts unchanged); the governance-acceptance workflow is the live enforcement path

- Direction 3 (EC-23): the all-refs scan's path-wide `generic-api-key` exemption under docs//evidence/ and the `|| true` error swallowing are REPLACED by an exact-fingerprint allowlist (`scripts/ci/ec23_allowlist.json`, 47 unique entries, census 48 raw findings over all refs on 2026-10-20, baseline report sha256 2d6bab7ebe5c395ac79c6f5e9a4ba22fa78c08ac8232a029ee583be4f923a919) validated by `scripts/ci/ec23_scan_check.py` (classify/positive-control/negative-control modes, all fail-closed); the positive control now asserts the expected rule+file finding AND verifies redaction (exit status alone is never acceptance); the H-1 class is kept separate from the digest false-positive classes; a passing scan is not proof of A-02 revocation. Disposition record: `docs/evidence/ec16-ec26/GITLEAKS_DISPOSITION_R19.md`.
- Direction 4 (EC-18/24): `evidence_sequence_check.py` rewritten to perform the protocol's LIVE remote verification (commits-API existence with exact returned-SHA match; 404/422 = ABSENT; 401/403/429/5xx/network = VERIFICATION_FAILURE fail-closed; branch-tip/ancestry checks), deterministic FIFO pairing with correction-entry semantics, reconciliation integrity (exact entry_id+branch+sha binding, machine-checkable pointers, SUPERSEDED-BY-DOCUMENTED-REBUILD still checks remote existence, SUPERSEDED-BY-CORRECTION-ENTRY names its authoritative correction), and gate mode (`--require-complete`: PENDING and OFFLINE-CLEAN fail; COMPLETE only with live verification). Records: `PUSH_LOG_RECONCILIATION.json` + `.md` (six evidence-verified rebuild entries; every recorded SHA verified this session from fetched objects).
- Direction 5 (EC-26): `reading_note_map_check.py` now validates trigger applicability for ACTIVE workflow points (workflow_dispatch-only, missing/empty triggers, finite branch lists without '**', tags-only push, finite paths filters all fail), retains the R19 unpinned/live-path rules, and computes the Enforcement class (LIVE/PROCEDURAL/INERT/NONE) validated against the declared label in the 5-column INDEX inventory — policy status and automation enforcement status are recorded separately.
- Environment anomaly recorded honestly: reads of `proposals/remediation/toscanini-gate.branch-triggers.yml` (blob 3e9e0b4d) and the old fixture returned content that differed between read paths in this sandbox; the GIT BLOBS are the canonical evidence (codepoint-level verification in-process). Because the file could not be reliably edited here (writes reverted to the committed blob), the superseded draft was REMOVED (append-only: v1 remains in history at blob 3e9e0b4d) and the corrected draft was written to a NEW path `toscanini-gate.branch-triggers-v2.yml` (pull_request.branches repaired to the continuous set ["**"], matching the corrected live state); the defect fixture `toscanini_gate_AIN_DEFECT.yml` was recreated byte-faithfully (defect literal constructed via codepoints) after the tracked copy was found normalized.
- Suite soundness fix: a conftest autouse fixture now fails any pytest node whose check() assertions recorded failures — the prior suite could report green under pytest while check() had recorded failures (the PR #30 tip f44bf4d9 carried 11 such hidden failures: AT-02/04/09/10/11/12/18 classes; all are genuinely green now, with tests corrected to the checkers' real contracts: file-path --api-json/--runs-manifest fixtures, full directive bodies for EC-20 evidence runs, timestamp-collision-free synthetic repos).

## 2026-10-20T18:40:00Z | R14 | record | Round-20-fix soundness cycle on this repair branch (candidate 62fab8a): the count record is brought current — 10 checker scripts under scripts/ci/ at that point (the R20 set with secret_baseline_check.py + secret_control_check.py; the R19 interim ec23_scan_check.py + ec23_allowlist.json superseded by the reviewed baseline) and 7 workflow drafts under proposals/ci/ (unchanged)

- Suite soundness: the R20 rewrite of `test_acceptance.py` lost the Round-19 conftest soundness fixture; on the exact candidate 62fab8a a pristine pytest run reported `24 passed` while a conftest probe surfaced 12 test nodes carrying hidden check() failures. Repaired in the Round-21-followup cycle (R15 below) on top of the R21 candidate.

## 2026-10-20T18:55:00Z | R15 | record | Round-21-followup suite-soundness cycle on this repair branch (candidate 0fe4840 + this push): the count record is brought current — 11 checker scripts under scripts/ci/ (the R21 set adds generate_evidence_table.py) and 7 workflow drafts under proposals/ci/ (unchanged); the governance-acceptance workflow is the live enforcement path

- The R21 commit (0fe4840) repaired the three Round-21 acceptance defects but did NOT restore the conftest soundness fixture; a pristine pytest run on the exact candidate reported `29 passed` while a conftest probe surfaced 15 test nodes carrying hidden check() failures (AT-04/09/10/11/12/14/16/17/18/21/23 R20 classes + AT-24/26/27/28 R21 classes). This cycle restores `scripts/ci/tests/conftest.py` (accumulate-and-report semantics preserved) and repairs every node to its checker's real contract: file-path `--api-json`/`--runs-manifest` fixtures; fail-closed EC-16/EC-19 fields in EC-20 evidence bodies; derived baseline-entry fields; event-annotated ACTIVE citations; a genuine orphan-line STACK_MISMATCH fixture; the true-defect AIN fixture restored; the draft's `on:` block genuinely repaired; the `main]`-substring collision fixed by precise token matching; salted synthetic repos (timestamp-collision class); AT-23 gains an explicit negative check. No check was weakened or deleted (EC §13).
- Checker hardenings (EC-24): the claimed branch must be named IN THE EVIDENCE FIELD (the body fallback made the check vacuous); the terminal intent's declared pre-push tip is validated against the newest recorded branch state INDEPENDENTLY of the ancestry check (both findings can fire); a missing log file fails closed with a clear message instead of an unhandled traceback. The real PUSH_LOG passes the hardened checks with the CI-equivalent invocation.
- Environment anomaly disclosed (recurring from R13/R14): reads of the draft `toscanini-gate.branch-triggers.yml` (blob 3e9e0b4d) again returned differing content across read paths in this sandbox; GIT BLOBS remain canonical and the repaired draft was verified at codepoint level in-process before commit.

## 2026-10-20T20:10:00Z | R16 | record | Round-22 acceptance-defect cycle on this repair branch (candidate built on 908de6d): the count record is unchanged — 11 checker scripts under scripts/ci/ and 7 workflow drafts under proposals/ci/; the governance-acceptance workflow is the live enforcement path and now carries 10 jobs (7 always-run + 3 PR-only, adding pre-merge-verification)

- Auditor Round-22 finding (direction 2, EC-23): the Round-21 fix removed the
  direct `Secret`/`Match` prints but still interpolated report-derived
  `RuleID`/`File`/`StartLine` values into failure diagnostics and
  negative-control output. Repaired: every report element is type-validated
  before examination (non-object rows, non-string RuleID/File, non-int
  StartLine, non-string Secret/Match fail closed with fixed messages and row
  contents withheld); unexpected findings report a COUNT only; the raw-key
  checks cover the report bytes even when the schema fails; unexpected
  internal errors fail closed with a fixed message. AT-30 proves the key
  never reaches stdout/stderr on every path. Local re-run of the pinned
  scanner controls on this cycle's tree: positive control PASS (genuine
  detection, redaction verified), negative control PASS (verified empty
  report), HEAD-tree scan PASS (verified empty report + fail-closed capture
  log), all-refs scan PASS (48 findings / 47 reviewed fingerprints; the two
  A-02 history entries disclosed separately).
- Auditor Round-22 findings (directions 3-4, EC-24): commit existence was
  conflated with branch-push provenance; cited runs' event/head_branch were
  not validated; intent reconciliations lacked branch scope, one-void, and
  consumed-intent rules; the unmatched-intent invariant was branch-local;
  expected_remote_tip_before_push was not bound to push-event before/after
  state. Repaired per schema v5 (see PUSH_LOG_SCHEMA.md): claim separation
  (PUSH-PROVEN-BY-RUN / PUSH-PROVEN-BY-EVENT / NOT-PROVEN reclassification /
  COMMIT-EXISTS-ONLY fail-closed), event+head_branch validation on every
  citation, `review_ref=` adjudication references validated against this
  register, global unmatched-intent invariant, PushEvent before/after binding
  (all-zero before for ABSENT; prior unlogged updates and false ABSENT
  fail). Verified against the LIVE log with a CI-equivalent manifest: P-012/
  P-014-CONFIRM classify PUSH-PROVEN-BY-RUN (their failed conclusions
  disclosed separately); P-011/P-013-CONFIRM fail closed pending this
  cycle's re-adjudication records (register Part 14.4).
- Auditor Round-22 findings (direction 5, merge transition): the workflow
  lacked `pull-requests: read`; the checker trusted merge-commit parents
  without merge_commit_sha/source-ref/merged_at/fresh-approval proof; no
  pre-merge gate existed; DEFERRED supersession was broken. Repaired:
  permission added (least privilege: read-only) with live API canaries in
  both push-log-integrity and the NEW pre-merge-verification job; post-merge
  verification checks merge_commit_sha == live candidate, head.ref ==
  source_branch, merged_at present, latest-state-per-reviewer approvals
  bound to the exact head (commit_id) and submitted before the merge
  (dismissed/superseded rejected), and the merge PushEvent's before ==
  declared base; the pre-merge-verification job verifies the CURRENT PR/
  head/review state on PR events and is EXPECTED RED until the independent
  approval exists (the functioning pre-merge gate; disclosed in P0_STATUS
  §12 and the reading note); DEFERRED supersession resolves via introducing
  commits on both source and target sides. AT-33 covers every failure mode
  plus a sequential second merge.
- Auditor Round-22 findings (direction 6, EC-26): representative-sample
  coverage could be satisfied by a single matching sample; ordered negative
  patterns and empty/malformed filters were mishandled. Repaired:
  conservative witness-based glob coverage (structural witnesses from each
  governed glob: multiple name shapes per wildcard component, `**` at depths
  1-3, root/hidden/workflow paths for whole-tree governance), ordered `!`
  semantics for paths AND branches (last match wins), empty filter lists and
  malformed event configuration types fail closed. AT-34 covers the
  auditor's full list.
- Acceptance suite: AT-01..AT-34 registered and green (34/34 pytest +
  run_all; no check weakened or deleted; AT-27/AT-28 fixtures strengthened
  to the new required PR fields merge_commit_sha/merged_at/commit_id/
  submitted_at — the old fixtures lacked exactly the fields the Round-22
  audit found unverified).

## 2026-10-20T23:20:00Z | R17 | record | Round-24 cycle (branch `repair/round24-claim-head-binding`, based on the merged `main`
`63da0e347002ed9b77071ff1ed0cefda7c0ce88a`): the merged Round-22/Round-23 implementation was audited against the
auditor's defect list (EC-23, EC-24 and EC-26 are CLOSED there and reproduce
against the live remote), and the two defects still open in the merged state were
repaired in ONE ordered commit: (a) the EC-16/19/20 exact-head binding --
`verify_claimed_shas.py --expect-head-sha` (supplied by the workflow from
`github.event.pull_request.head.sha`) requires the PR body to DECLARE the exact
current head and to cite a real, successful run ON that head, and the workflow
reads the LIVE pull-request body; AT-35 is the adversarial test and FAILS on the
prior implementation; (b) the push-log gate's introducing-commit lookup, whose
git pathspec was built from the OS path string and fell outside the repository
where the path spelling differs from git's (Windows short form vs long form) --
the lookup now runs from the log's own directory; (c) AT-13's append-only
comparison copy is now written in binary mode (text mode translated the line
endings twice on Windows and produced a false "SHRUNK" failure). Register Part 15
records the adjudication table and the OPEN item: the merge transition of PR #30
is NOT verified by this repository's own gate (no independent approval; the
declared pre-merge base is no longer the live tip of `main`), so `main`'s
Governance Acceptance run 39504181997 is RED and the repair branch's last push
(d81c874, run 39507046623) is RED for the same reason -- disclosed, never
suppressed, no approval invented, no gate weakened. The count record is brought
current: 11 checker scripts under scripts/ci/ (unchanged this cycle) and 7
workflow drafts under proposals/ci/ (unchanged); the acceptance suite is AT-01..AT-35
with AT-35 registered and mapped in the INDEX rows for EC-16/19/20. No token
value appears in this receipt, in the repository, or in any worklog.

## 2026-10-20T23:45:00Z | R18 | record | Round-24 cycle, second push: the package manifest had been
regenerated over the WORKING-TREE bytes of a CRLF checkout while the repository
stores LF blobs, so every digest mismatched on the Linux runner and the
package-integrity job of the first push (Governance Acceptance run 39515766450,
Stability Gate run 39515766416 conclusion=success) FAILED. Repair: the manifest
is regenerated over the COMMITTED bytes and `scripts/ci/manifest_check.py` now
hashes the committed artifact (the index entry) instead of the materialised
checkout, so the check is platform-independent and pins the shipped artifact.
The failed run is disclosed in P-082-CONFIRM; no digest was fudged, no scope
file was dropped, and the failure is not smoothed over. The push-log protocol
records for this push are P-082-CONFIRM (the first push of the branch) and
P-083-INTENT (the terminal intent bound to the exact candidate with the true
pre-push tip 6ff3402ae7e13b51a131e9c760f9ae0cb1376f1f).

## 2026-10-21T00:10:00Z | R19 | record | Round-25 cycle on `repair/round24-claim-head-binding` (candidate on top of
2ce77e57194702c2b292fc3f3b28cddd2cff15c9): the auditor's Round-25 decision (PR #37 not approved for merge; the
live-body fallback, the claimed-CI-evidence validation, the manifest semantics
and the missing pre-merge gate are blockers). Repaired: (a) the live PR body is
now retrieved by `scripts/ci/fetch_pr_body.py` and FAILS CLOSED (the
event-payload fallback is removed; empty, malformed, null and missing bodies
all fail; the expected head stays an explicit input) - AT-36; (b) the runs
manifest carries the workflow identity and the claim gate requires the
`Toscanini Stability Gate` and `Governance Acceptance` runs explicitly (an
unrelated successful workflow, a failed/cancelled required workflow, a
required workflow on another SHA and an unverifiable run are all rejected) -
AT-37; (c) the package manifest validates the COMMITTED HEAD tree and a
staged index is rejected - AT-38; (d) the executable pre-merge-verification
job is RESTORED from the verified pre-R23 blob 44fe653e (the Round-24
description of the Round-23 removal as sanctioned is superseded, register
16.4); (e) the PR #30 merge incident is recorded without rewriting history
(register 16.5); (f) the preservation bundle is generated by automation,
published as a CI artifact and re-verified from a fresh checkout (register
16.6); (g) A-04 is verified through an authorized route (register 16.7). The
count record is brought current: 13 checker scripts under scripts/ci/ and 7
workflow drafts under proposals/ci/; the acceptance suite is AT-01..AT-38 with
AT-36/AT-37/AT-38 registered and mapped. No token value appears in this
receipt, in the repository, or in any worklog.
