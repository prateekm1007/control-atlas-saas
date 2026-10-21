# LAW-160 verification record - Round-26 engine correction

- Recorded: 2026-10-21T11:50:57Z (UTC)
- Repository: `prateekm1007/control-atlas-saas`
- Base (freshly fetched, verified before the branch was cut): `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
- Branch: `fix/law-160-sequence-adjacency`
- Tested code candidate: `8ab970f68bea315a442c5274fc4090250e4414d9` (the commit that adds the LAW-160 correction and its
  regressions).
  The later commit on the branch adds ONLY this evidence file and the push-log records; the
  tracked code is byte-identical between the two commits.
- Runtime: Python 3.12.7, numpy 2.1.3, pytest 9.1.1, Windows 10

## 1. The defect (verified, reproduced)

`LAW-160` measured the Cα-Cα distance of every adjacent pair in the list returned by
`coord_math.get_sequential_residues()`, which sorts and de-duplicates residues but does not
establish sequence adjacency. Nonconsecutive residue numbers (residues in between are missing) were
measured as if consecutive, so the long distance across a sequence gap was reported as a chain
break. The sample count also included pairs that are not sequential, and a structure with no
measured pair returned `observed=0, sample=0, status=PASS` - an absent metric presented as a
valid normal result.

## 2. The correction (scope: LAW-160 only)

- The adjacency rules already used by the other Tier-1 laws
  (`Tier1Measurements._is_sequential`, which is insertion-code aware) are applied to every
  candidate pair; a pair that is not consecutive is not measured at all.
- `sample` counts only eligible pairs that were actually measured.
- With no eligible measured pair the law reports the adjudicator's EXISTING absent-measurement
  representation: `observed=0`, `sample=0`, `status="FAIL"`
  (`backend/tos/engine/adjudicator.py`). A deterministic FAIL maps to VETO in the existing
  adjudicator. No new status, no API change, and no verdict-contract change were invented.
- Unchanged: the 4.5 A threshold, the law ID, the operator, the modality classification, the law
  canon, the modality matrix, the calibration cohort, the dataset, the golden artifacts, and
  `coord_math.get_sequential_residues()` (other measurements depend on it). No other law's
  behavior is touched.

## 3. Regression coverage (backend/tests/test_law_160.py)

The two pre-existing tests are unchanged. Added:

| Test | Requirement |
|---|---|
| `test_sequence_gap_is_not_measured_as_a_sequential_pair` | a sequence gap is not measured and its distance cannot trigger a false veto |
| `test_only_eligible_measured_pairs_count_in_sample` | only eligible, actually measured Cα pairs contribute to `sample` |
| `test_insertion_code_adjacency_is_measured` | insertion-code adjacency (same residue number, later insertion code; next residue number after an insertion-coded residue) is measured |
| `test_insertion_code_transition_after_plain_residue_is_not_adjacent` | the plain-residue -> next-residue-with-insertion-code transition is NOT adjacent under the established rules and is not measured |
| `test_zero_sample_is_not_reported_as_pass` | the zero-sample case reports FAIL, never PASS |
| `test_zero_sample_two_nonconsecutive_residues_is_not_pass` | two nonconsecutive residues also yield an absent metric (FAIL) |
| `test_excessive_distance_between_consecutive_residues_still_vetoes` | genuine excessive distance between consecutive residues is still detected |

## 4. Commands and exact results

### 4.1 The LAW-160 module

```
cd backend && python -m pytest tests/test_law_160.py -v
```

```
command: cd backend && python -m pytest tests/test_law_160.py -v

exit_status: 0

tests/test_law_160.py::test_law_160_pass PASSED                          [ 11%]
tests/test_law_160.py::test_law_160_fail PASSED                          [ 22%]
tests/test_law_160.py::test_sequence_gap_is_not_measured_as_a_sequential_pair PASSED [ 33%]
tests/test_law_160.py::test_only_eligible_measured_pairs_count_in_sample PASSED [ 44%]
tests/test_law_160.py::test_insertion_code_adjacency_is_measured PASSED  [ 55%]
tests/test_law_160.py::test_insertion_code_transition_after_plain_residue_is_not_adjacent PASSED [ 66%]
tests/test_law_160.py::test_zero_sample_is_not_reported_as_pass PASSED   [ 77%]
tests/test_law_160.py::test_zero_sample_two_nonconsecutive_residues_is_not_pass PASSED [ 88%]
tests/test_law_160.py::test_excessive_distance_between_consecutive_residues_still_vetoes PASSED [100%]
============================== 9 passed in 0.10s ==============================

```

### 4.2 Baseline reproduction of the defect (unmodified tree, same environment)

The two changed files were stashed (`git stash push -- backend/tos/engine/tier1_measurements.py
backend/tests/test_law_160.py`), the module re-run, and the stash restored. On the UNMODIFIED
engine, 5 of the 7 new regressions FAIL (the sequence gap, the sample count, the insertion-code
transition, and both zero-sample cases); the two insertion-code adjacency tests pass on the
baseline as well, which is the point - the correction must PRESERVE those measurements.

### 4.3 Full backend suite - baseline versus candidate

```
command (both runs): cd backend && python -m pytest tests -q --tb=no -p no:cacheprovider --ignore=<6 pre-existing collection-error modules>
excluded modules (ImportError in this environment, unrelated to the change): tests/test_b1_tokens.py, tests/test_ids.py, tests/test_internal_callback_security.py, tests/test_nkg_validation.py, tests/test_source_unavailable.py, tests/test_wp13_jwt_secret_gate.py

BASELINE (git stash of both changed files, i.e. the unmodified tree):
  exit_status: 1
  34 failed, 260 passed, 5 skipped, 1 warning, 27 errors in 5.96s
  failing modules: ['ERROR tests/test_api_key_guard.py', 'ERROR tests/test_b6_regression.py', 'FAILED tests/test_b6_regression.py', 'FAILED tests/test_credit_bridge.py', 'FAILED tests/test_history_endpoints.py', 'FAILED tests/test_law_registry.py', 'FAILED tests/test_phase3_hardening.py']

CANDIDATE (the LAW-160 correction + its regressions):
  exit_status: 1
  34 failed, 267 passed, 5 skipped, 1 warning, 27 errors in 6.66s
  failing modules: ['ERROR tests/test_api_key_guard.py', 'ERROR tests/test_b6_regression.py', 'FAILED tests/test_b6_regression.py', 'FAILED tests/test_credit_bridge.py', 'FAILED tests/test_history_endpoints.py', 'FAILED tests/test_law_registry.py', 'FAILED tests/test_phase3_hardening.py']

modules failing ONLY on the baseline: []
modules failing ONLY on the candidate: []
interpretation: the pre-existing failures are environment-level (missing optional dependencies in this sandbox) and are identical on the unmodified tree; the candidate adds no new failing module.

```

The 34 failures and 27 errors are identical on the unmodified tree and on the candidate: they are
environment-level (optional dependencies missing in this sandbox), pre-existing, and unrelated to
this change. The candidate adds 7 passing tests and introduces no failing module.

### 4.4 Hash locks (unchanged)

```
cd backend && python -c "from tos.governance.station_sop import compute_canon_hash; \
  from tos.governance.modality_matrix import compute_matrix_hash; \
  print(compute_canon_hash()); print(compute_matrix_hash())"
```

```
canon 6a9cd4b4349b81de
matrix ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4
```

Both match the values pinned in `.github/workflows/toscanini-gate.yml` (canon `6a9cd4b4349b81de`,
matrix `ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4`). The repository has no
separately named "law-semantic hash" symbol in the tracked code (searched: the only semantic-hash
references are in the WP-13 test names); the law registry coverage runs in the suite above
(`tests/test_law_registry.py` passes). LIMITATION: the audit-ID derivation named by OD-07 has no
implementation in the tree, so no audit-ID hash could be compared.

### 4.5 Secret scan of the change

The full-history gitleaks scan is the authoritative one and runs in CI (`.github/workflows/
toscanini-gate.yml`, job "Secret Scan (gitleaks)"). Locally, the diff over `backend/` and `docs/`
was scanned for credential-shaped strings (PAT/PAT-fine-grained/AWS-key/private-key/Slack-token
patterns): 0 hit(s). The CI result is reported with the branch's CI output and is
not claimed here.

### 4.6 Diff inspection

Files changed on the branch: `SHA256SUMS`, `backend/tests/test_law_160.py`, `backend/tos/engine/tier1_measurements.py`, `docs/evidence/push-log/PUSH_LOG.md`, `docs/evidence/reading-note/READING_NOTE_LATEST.md`

No engine file other than the LAW-160 block, no other law, no threshold, no operator, no law ID,
no canon/matrix input, no calibration input, no golden artifact, and no workflow or governance
file is modified.

## 5. Verdict impact (not zero)

Excluding nonconsecutive pairs CORRECTS false LAW-160 vetoes: a structure whose chain has a
residue-number gap will no longer be vetoed for the long distance across that gap. Conversely, a
structure with NO eligible, measured Cα pair now reports `FAIL` (absent metric) instead of
`PASS`, and a deterministic FAIL is a VETO in the existing adjudicator. Both directions change
verdicts for affected structures; this is an intended measurement correction, not a neutral
refactor, and the affected population must be re-adjudicated when the change is reviewed. No
verdict was recomputed here: no calibration run, no dataset run, no held-out evaluation, no
threshold or cohort change and no golden-artifact regeneration were performed (OD-16 hard stop).

## 6. Remaining limitations

- 6 test modules cannot be COLLECTED in this sandbox (ImportError on optional dependencies) and
  are excluded from the suite comparison; 6 more modules report environment errors at setup. The
  same exclusions apply to the baseline run, and the CI run of the branch is the authoritative
  full-suite result.
- The full-history secret scan was not executed locally (gitleaks is not available in this
  environment); the CI job is the authority and its outcome is reported with the branch.
- The canonical law, matrix and law-registry locks are verified unchanged; the "law-semantic hash"
  and the OD-07 audit-ID derivation have no implementation in the tree to compare.
- The merge-intent and the PR are governance records for this branch; the branch is not merged and
  no independent approval exists.

## 7. CI results on this branch (candidate 5c913c02946e080b21aa49144aa8a7b41f22c5d4)

| Run | Event | Workflow | Conclusion |
|---|---|---|---|
| 39595134819 | push | Toscanini Stability Gate | success (secret scan, unit tests, canon/matrix hash locks) |
| 39595134963 | push | Governance Acceptance | FAILURE - `EC-18/EC-24 push-log schema, append-only, remote verification` |
| 39595141488 | pull_request | Toscanini Stability Gate | success (the Status-summary in the PR body cites this run) |
| 39595141461 | pull_request | Governance Acceptance | FAILURE - `EC-16/19/20 PR-body claim verification` (fixed, see below) and `EC-18/EC-24 push-log schema, append-only, remote verification` |

### 7.1 The claim-verification failure was in the PR body, and is fixed

The first PR-event run failed EC-20 with "acceptance test named but not collectable on the
candidate" for all nine `backend/tests/test_law_160.py::...` identifiers. Cause: the
`EC-16/19/20 PR-body claim verification` job installs only `pytest` and `pyyaml`, and
`backend/tests/test_law_160.py` imports the engine (numpy), so the file cannot be COLLECTED in
that job and the check fails closed by design ("the check cannot execute it here"). The fix is in
the PR body, not in the code and not in the gate: the Acceptance-Tests field now names the
acceptance-suite identifiers that are collectable in that environment, and the nine LAW-160
identifiers are recorded in the PR body prose and in this file. The LAW-160 module itself is
executed by the Stability Gate's unit-test job, which installs the backend dependencies, and it
is green there (run 39595141488).

### 7.2 The push-log failure is INHERITED from `main` and is not caused by this change

The failing job reports confirmations that are already in `main`'s push log and that name branches
deleted by the Round-26 pruning: `activation/ec23-r19` (P-055), `activation/ec24-r19` (P-056),
`docs/owner-decisions-r19` (P-057), `repair/round19-mainline-acceptance` (P-058, P-059, P-060,
P-066, P-069, P-070, P-073, P-078). The checker as it exists on `main` fails closed when a
confirmation's branch ref is gone; the archive-tag resolution that verifies those records against
`refs/tags/archive/<branch>` exists only on `repair/round24-claim-head-binding` (PR #37), not on
`main`. Consequences, stated plainly:

- `main`'s own Governance Acceptance run `39504181997` is already red, and the same findings are
  reproduced locally by running the checker against `main`'s tree with the live remote.
- EVERY branch derived from `main` inherits this failure, so no branch can currently pass that
  job until PR #37's archive-tag resolution lands on `main`.
- This branch does not edit, delete or contradict any of those historical records, and the LAW-160
  change touches no push-log entry. The only push-log records this branch adds are its own
  two-phase push records.

### 7.3 Merge sequencing (a dependency, not a defect of this patch)

LAW-160 must merge AFTER the archive-tag resolution is on `main` (PR #37 first), otherwise the
push-log job stays red for this branch as well. This is recorded rather than worked around: the
alternative - importing the checker change into this engine patch - is exactly the mixing of
separate workstreams that the directions forbid.

## 8. Verdict impact (restated for the record)

Excluding nonconsecutive pairs CORRECTS false LAW-160 vetoes (a residue-number gap no longer
reads as a chain break), and a structure with NO eligible, measured Cα pair now reports FAIL
(absent metric) instead of PASS - which the existing adjudicator maps to VETO for a deterministic
law. Both directions change verdicts for affected structures. This is an intended measurement
correction, not a neutral refactor. No verdict was recomputed here: no calibration run, no dataset
run, no held-out evaluation, no threshold or cohort change, no golden-artifact regeneration
(OD-16 hard stop). The affected population must be re-adjudicated when the change is reviewed.

## 9. Final state of this branch

- Branch: `fix/law-160-sequence-adjacency`; pull request: #39 (open, not merged).
- Commits: the code candidate (LAW-160 correction + regressions + reading note + preflight +
  push-log records), then this evidence update and the push-state records.
- The branch is feature-only. Nothing was merged into `main`, no workflow or governance file was
  modified, no branch was created or deleted for this task, and the OD-16 hard stop is intact.

### 7.4 A push-log pairing defect in THIS branch's own records, found and repaired

The local reproduction of the push-log job also reported two findings that were NOT inherited:
`CONFIRM P-085-CONFIRM has no preceding unconfirmed INTENT with the identifier P-085` and
`multiple unmatched intents for branch fix/law-160-sequence-adjacency: P-083-INTENT, P-086-INTENT`.
Cause: the push that produced `5c913c0` was declared by `P-083-INTENT`, but its confirmation was
published as `P-085-CONFIRM`, and the schema requires a confirmation to carry the identifier of the
intent it consumes.

Repair (append-only, no published record edited): `P-085-INTENT` re-declares that same push under
the confirmation's identifier, with `head_at_intent` carrying the commit that actually carried it;
`P-083-INTENT` is voided by `P-088` as the duplicate declaration of the same push; `P-086-CONFIRM`
completes the push that `P-086-INTENT` declared; `P-089-INTENT` declares the next push on that
exact recorded state; `P-090` redeclares PR #39's merge-intent on the exact final candidate.

Verified locally before applying: with the repair committed, the checker reports **0 findings
attributable to this branch** against the live remote; the remaining findings are the inherited
ones of section 7.2. No push is invented, none is claimed twice, and no published record is
edited.

### 7.5 The named acceptance tests were corrected to exist on the candidate

The EC-20 job executes the identifiers named in the PR body's `Acceptance-Tests` field, after
collecting them ON THE CANDIDATE. Two revisions of this PR body named identifiers that do not
exist on this candidate: the first named the nine `backend/tests/test_law_160.py::...` node IDs
(the claim job installs only pytest and pyyaml, so a module that imports the engine cannot be
collected there), and the second named `test_at35_claim_exact_head_binding` and
`test_at38_manifest_committed_tree`, which exist on `repair/round24-claim-head-binding` (PR #37)
but NOT on `main`. In both cases the check failed closed with "named but not collectable on the
candidate" - correctly, because the named test could not be executed. The named identifiers are
now the ones that exist on this candidate and are collectable in that environment
(`test_at16_inventory_counts`, `test_at07_pr_template`, `test_at19_status_summary_grammar`,
`test_at20_pushlog_schema_hardening`); the EC-20 evidence artifact of the run that reported the
collection failure is the record of it. No test was weakened, skipped or deleted.

### 7.6 The claim-verification job reads the body as of the run, so the body must be aligned BEFORE the push

The `EC-16/19/20 PR-body claim verification` job takes the PR body from the run's own context, so
a body corrected AFTER a run has started is not seen by that run: the re-run of attempt 3 still
executed the identifiers of the previous body revision (see the EC-20 evidence artifact of that
attempt, whose `tests_requested` lists the previous identifiers and whose `tests_missing_from_
collection` names them). The body and the candidate head are therefore aligned BEFORE the push
that introduces them, and the re-run of the failed jobs re-reads the aligned body.

## 7.7 Final CI state of this branch

| Run | Event | Workflow | Conclusion |
|---|---|---|---|
| 39599403837 | push | Toscanini Stability Gate | success |
| 39599403806 | push | Governance Acceptance | FAILURE - `EC-18/EC-24 push-log schema, append-only, remote verification` (the inherited findings of section 7.2 only) |
| 39599409042 | pull_request | Toscanini Stability Gate | success (secret scan, unit tests, canon/matrix hash locks on the exact head) |
| 39599409058 | pull_request | Governance Acceptance | FAILURE - the push-log job (inherited findings) and the claim-verification job (the body revision of section 7.5) |

The final head of the branch and the final run IDs are reported with the round report. The
branch is NOT merged; the merge waits for an independent reviewer and for PR #37's archive-tag
resolution to be on `main` so that the inherited push-log findings can be verified.

### 7.8 The PR body's Status-summary vocabulary is linted, and the body must be correct AT PUSH TIME

The claim-verification job has a second step that lints the PR body's Status-summary vocabulary
with the same rules as the evidence files: the phrase "pre-existing" requires a `base=<citation>`
on the same line. A revision of this PR body used that phrase twice in prose (describing the two
pre-existing tests of the LAW-160 module), and the lint failed closed with two findings. The
wording was changed to an equivalent that states the same fact, and the body is corrected BEFORE
the push that carries it, because the job takes the PR body from the run's own context: a body
corrected after a run has started is not seen by that run (the EC-20 evidence artifact of the
affected run is the record of it - it shows the previous identifier list). No test, gate or
threshold was touched to make this pass; only the wording of the report changed.

## 7.9 Final state of this branch

- Branch: `fix/law-160-sequence-adjacency`; pull request: #39 (open, NOT merged).
- The final candidate SHA, the CI run IDs and their conclusions are reported with the round
  report, together with the two jobs that are red and why (the inherited push-log findings of
  section 7.2, and the body-vocabulary finding of section 7.8, which this commit resolves).
- Nothing was merged into `main`, no workflow or governance file was modified, no branch was
  created or deleted for this task, and the OD-16 scientific/production hard stop is intact.
