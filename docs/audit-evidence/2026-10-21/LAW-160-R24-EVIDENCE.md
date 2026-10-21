# LAW-160 (Chain Integrity) - sequence-adjacency correction: audit evidence

**Cycle:** Round-24 (auditor decision: stop the governance-only workstream and
resume the LAW-160 engine correction on a new, isolated feature branch).
**Branch:** `feature/law-160-sequence-adjacency` (created from a freshly
verified `main`; the branch did not exist on the remote before this push).
**Base commit:** `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (`main`).
**Candidate:** the commit that carries this evidence file. A commit cannot
contain its own SHA, so the candidate is bound by the push-log terminal intent
`P-082-INTENT` and, after the push returns, by the `P-082-CONFIRM` record and
the pull-request body (the repository's established two-phase convention).

**Runtime:** Python 3.12.7 on Windows 10 (AMD64); pytest as invoked above. The Linux CI
runner executes the same commands with the full suite (the five modules
excluded locally import `fcntl`, a POSIX-only module, and cannot be collected
on Windows).

---

## 1. The defect and the correction

`get_sequential_residues()` returns the residues that are PRESENT in a chain
(sorted by `(res_seq, insertion_code)`, de-duplicated). It does NOT guarantee
that two successive entries are consecutive in the biological sequence. LAW-160
measured the C-alpha distance between successive ENTRIES, so a chain with
missing intervening residues (a modelling gap) produced a large C-alpha
distance across the gap and was reported as a physical chain break: a FALSE
VETO.

The correction:

1. the adjacency predicate is factored into a single source of truth
   (`Tier1Measurements._adjacent(s1, i1, s2, i2)`); the existing
   `_is_sequential()` (residue dictionaries) and the new
   `_is_sequential_residues()` (the `(res_seq, insertion_code, res_name)`
   tuples of `get_sequential_residues()`) both delegate to it, so the two
   representations cannot drift apart;
2. the semantics of `_is_sequential()` are PRESERVED verbatim (same three
   clauses, same chain check) - only the factoring is new;
3. LAW-160 now measures only pairs that are sequence-adjacent; a pair with a
   residue-number gap is skipped;
4. `get_sequential_residues()` is NOT changed (other laws use it);
5. the canonical 4.5 A threshold, the law ID, the operator, the units and the
   classification are UNCHANGED (`backend/tos/governance/station_sop.py` was
   not touched).

## 2. Verdict impact (explicit)

Correcting pair eligibility can REMOVE false LAW-160 vetoes, and therefore
CHANGE the verdict of affected structures: a structure whose chain has missing
intervening residues previously vetoed on the modelling gap and now passes
(or, when no pair is eligible at all, is reported as unmeasured - see 5). The
Stability Gate / Governance Acceptance results and any structure-level
verdicts produced before this correction are not comparable with the
post-correction results for chains with gaps. The threshold, the law
semantics and every other law are unaffected.

## 3. Commands, exit codes and outputs

### 3.1 Targeted LAW-160 tests

```
$ cd backend && python -m pytest tests/test_law_160.py -q --no-header
```
exit code: 0

```
.........                                                                [100%]
9 passed in 4.36s
```

### 3.2 Full backend suite (candidate)

```
$ cd backend && python -m pytest tests -q --no-header --ignore=tests/test_internal_callback_security.py --ignore=tests/test_nkg_validation.py --ignore=tests/test_source_unavailable.py --ignore=tests/test_wp13_jwt_secret_gate.py --ignore=tests/test_ids.py
```
exit code: 1

```
_pinned_labels
FAILED tests/test_law_registry.py::TestHandoffMatchesRegistry::test_handoff_has_no_conflict_markers
FAILED tests/test_modality_matrix.py::TestMatrixHash::test_hash_matches_known_value
FAILED tests/test_phase3_hardening.py::TestGovernanceFingerprintStructure::test_matrix_hash_matches_locked_value
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _health]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _laws]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _definitions]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[POST _ingest]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[POST _v1_batch]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _audit_{audit_id}]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[POST _search]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _nkg]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _comparison_{o}_{r}]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _comparisons_by-original_{id}]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[POST _compare]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _refinement_pdb_{audit_id}]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[POST _refinement_token]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[POST _refinement_callback]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[POST _refinement_submit]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _refinement_status_{job_id}]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _refinement_jobs_{audit_id}]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _history_{user_email:path}]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _usage]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[GET _history_ip_1.2.3.4]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_route_auth[POST _refinement_internal-callback]
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_all_app_routes_covered
ERROR tests/test_api_key_guard.py::TestRouteCoverage::test_invalid_key_rejected_on_protected_route
25 failed, 291 passed, 5 skipped, 1 warning, 23 errors in 56.46s
```

**Baseline (unmodified tree, same environment, same command):** the identical
suite reports `25 failed, 284 passed, 5 skipped, 23 errors` - the same
failing/erroring set (29 unique ids) as the candidate. The candidate reports
`25 failed, 291 passed, 5 skipped, 23 errors`: the delta is exactly the seven
new LAW-160 regression tests. The pre-existing failures and errors are
UNRELATED to this patch (they are Windows import/API issues in other modules)
and are reported, not suppressed.

### 3.3 Hash locks

```
$ cd backend && python -c "from tos.governance.station_sop import compute_canon_hash; ..."
```
exit code: 0

```
canon : 6a9cd4b4349b81de
matrix: b1249cedde5018628ea2351de49c9622cdce4657561af7f74146df182d1a210d
```

The candidate's canon and matrix hashes are byte-identical to the base
commit's: the patch touches no governance artefact, and the Stability Gate on
the candidate (run `39598182836`, conclusion `success`) re-computed and
verified BOTH locks on the exact candidate.

A NOTE ON THE LOCAL MATRIX VALUE, so the record is not misread: the matrix hash
computed on THIS Windows host (`b1249cedde5018628ea2351de49c9622cdce4657561af7f74146df182d1a210d`)
differs from the value the Stability Gate asserts
(`ed1731ace09a9a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4`). That
difference is a CHECKOUT-ARTEFACT of this platform, not a repository drift:
the same local value is computed on the unmodified base commit, and the Linux
runner computes the asserted value on the same tree. The hash inputs are
checked-out files whose line endings differ between a CRLF and an LF checkout.
The same artefact explains why the local `manifest_check.py` run reports 46
digest mismatches while the CI package-integrity job is green: the checker on
`main` hashes the WORKING TREE (a CRLF checkout here), while the manifest
pins the COMMITTED bytes. This is a PRE-EXISTING platform dependence of the
checker on `main` (a fix for it exists on the unmerged PR #37 branch and is
out of scope for this LAW-160 patch); it is reported here and NOT fixed here,
because the patch is limited to LAW-160 and the canon/matrix are frozen
artefacts.

### 3.4 Secret scan of the patch

The scan reads the patch as bytes and decodes it defensively (the diff
contains non-UTF-8 bytes from the surrounding evidence, which is why the
first attempt at this scan aborted on a decode error - the scan itself found
nothing).

```
$ python - <<'PY'
import re, subprocess
d = subprocess.run(["git","-C",REPO,"diff",BASE], capture_output=True).stdout.decode("utf-8","replace")
pat = re.compile(r"(ghp_[A-Za-z0-9]{36}|AKIA[0-9A-Z]{16}|sk-[A-Za-z0-9]{32,}|"
                 r"-----BEGIN [A-Z ]*PRIVATE KEY-----|xox[baprs]-[A-Za-z0-9-]{10,})")
hits = [l for l in d.splitlines() if pat.search(l)]
print("patch secret-pattern hits:", len(hits))
raise SystemExit(1 if hits else 0)
PY
```
exit code: 0

```
patch secret-pattern hits: 0
```

Patterns: GitHub PAT (`ghp_*`), AWS access key id (`AKIA*`), OpenAI-style
key (`sk-*`), PEM private-key header, Slack token. The tracked tree was also
scanned for a PAT-shaped value: none. The credential used for the push is
injected through the environment only and is never written to a file, a log, a
commit message or this evidence.

### 3.5 Changed files

```
$ git diff --stat 63da0e347002ed9b77071ff1ed0cefda7c0ce88a
```
```
$(subprocess.run(["git", "-C", REPO, "diff", "--stat", BASE], capture_output=True, text=True).stdout.strip())
```

The complete diff is the two source files above plus this evidence, the
reading record and the push-log intent. No other engine, law, threshold,
canon, matrix, calibration, deployment or dashboard file is touched.

## 4. Files changed and their content digests (SHA-256)

| File | SHA-256 |
|---|---|
| `backend/tos/engine/tier1_measurements.py` | `0d62dbad6bf4f165bdd58971316addb2b91f6b20d6313ed50291a7616b2d8bb4` |
| `backend/tests/test_law_160.py` | `b73b241aca229a0bf9c2ec62381b743957ea1a4f2e28341b5771c3827b1f3152` |
| `docs/evidence/push-log/PUSH_LOG.md` | `aba041345431398cc0d47f6e73f262826f16d8ab7baa829a3ba23711cde3b433` |
| `docs/evidence/reading-note/READING_NOTE_LATEST.md` | `e618604bfd0ebfaaa70deb69aa6eea0514db5587f35cdcfae4c372d6341f24c5` |

(base-commit digests for comparison: the engine file
`adc46abae2ee6d6c03ed3f0aaf2372d89e3ca6b42936c0fcb47a37062087c9e0`, the test
file `16d50669aa76f79462cd2fd936faa626c9c017afa6b1b6ae2e38b5e23254293e`.)

## 5. OPEN: the zero-eligible-pair contract (release blocker, not fixed here)

When no sequence-adjacent C-alpha pair exists (or no pair has both atoms),
LAW-160 is NOT measured. The previous code emitted `observed=0.0`, `PASS`,
`sample=0` - a fabricated value and a PASS for a law that was never
evaluated, which is exactly the misleading scientific result the auditor
identified.

**What this patch does.** The engine emits NO LAW-160 measurement in that
case. The adjudicator's EXISTING missing-measurement default then applies
(`backend/tos/engine/adjudicator.py`: `{"observed": 0, "deviation": "0.0",
"sample": 0, "status": "FAIL"}`) - the same representation every other
unmeasured law already receives. Nothing is fabricated, no new status is
invented, and no PASS is claimed.

**Why that is not the final answer.** "FAIL" for a law that could not be
measured still reads as a scientific failure, and it lowers the deterministic
score of an otherwise sound structure. The correct representation needs a
contract decision:

- **Proposal A (preferred):** extend the approved status vocabulary with an
  explicit unmeasured state (e.g. `NOT_MEASURED`) and make the adjudicator
  treat it as "not evaluated": excluded from `det_passed`, excluded from
  `failing_det`, and surfaced in the response with an explicit
  `measured: false` flag. Requires: an approved status token (the engine's
  vocabulary today is exactly `PASS` / `VETO` / `FAIL`), a change to
  `backend/tos/engine/adjudicator.py` (the status handling and the
  deterministic score), and a change to the API model
  (`backend/tos/schemas/response_v1.py`, where `PhysicsMeasurement.status` is a
  free string but `observed: float` is REQUIRED).
- **Proposal B (minimal, weaker):** keep `observed` a float and add a sibling
  boolean `measured` to the measurement row, with the status left as `FAIL`
  for backward compatibility. The response model needs the extra field; the
  consumer surface (the audit report, the PDF) needs to render it.

**Why `null` was rejected:** `PhysicsMeasurement.observed` is `float` and
required, so a null observation fails validation of the whole response.

**Status: OPEN.** This point requires an owner decision (the status
vocabulary and the response schema are approved artefacts). The LAW-160
correctness repair in this patch is complete and verified; the
zero-sample representation is NOT, and the task is therefore reported as
PARTIALLY COMPLETE with this item explicitly outstanding.

## 6. Unresolved limitations

1. The zero-eligible-pair representation (5) is unresolved and needs an
   owner decision.
2. The pre-existing matrix-hash drift between `main` and the Stability Gate
   assertion (3.3) is reported, not fixed (frozen artefact).
3. The pre-existing failures/errors in the backend suite (3.2) are reported
   with the baseline comparison; they are outside the LAW-160 scope.
4. The four POSIX-only test modules could not be collected on this Windows
   host; the Linux CI runs the complete suite on the exact candidate.
5. No production, deployment, calibration, v25, held-out or threshold work
   was performed, and nothing was merged.

## 7. CI on the exact candidate

| Run | Event | Conclusion |
|---|---|---|
| `39598182836` (Toscanini Stability Gate) | push | **success** - unit tests and the canon/matrix hash locks pass on the exact candidate |
| `39598182682` (Governance Acceptance) | push | failure - one job, see below |

The failing job is `EC-18/EC-24 push-log schema, append-only, remote
verification`. It is NOT caused by this patch: the push log on `main` records
confirmations for the branch `repair/round19-mainline-acceptance`, and that
branch was deleted by the repository-wide pruning that happened AFTER `main`'s
last Governance Acceptance run. Every affected entry fails with "remote ref
refs/heads/repair/round19-mainline-acceptance not found". The same defect is
inherent in `main`; a branch cut from `main` inherits it. It is a
push-log/governance matter, outside the LAW-160 scope of this patch, and it is
reported here rather than patched here.

Every other job of that run is green, including the acceptance suite, the
package manifest, the reading-note validation, the status-vocabulary lint, the
deploy-trigger lint, the full-history secret scan and the preservation-bundle
jobs.

## 8. Round-24 continuation: the approved explicit NOT_MEASURED state

**Decision (as given, by the auditor under the owner's delegation):** Option A
is approved - an explicitly unmeasured measurement gets an explicit
`NOT_MEASURED` state, a truthful typed response and conservative
adjudication. The owner delegated this decision to the auditor; this section
records the decision AS GIVEN. It is not an owner signature and it does not
introduce a new decision-register identifier.

**Starting SHA of this continuation:** `1831bdfa98180e9c95d64d5e0a8a7148e3624d4e` (the branch head that carried
the sequence-adjacency repair and the evidence correction). The patch below is
a focused continuation of the same branch; it is not merged.

### 8.1 What was implemented

1. **Engine** (`backend/tos/engine/tier1_measurements.py`): when no
   sequence-adjacent C-alpha pair with both atoms exists, LAW-160 returns an
   EXPLICIT row - `status="NOT_MEASURED"`, `measured=False`, `observed=None`,
   `deviation=None`, `sample=0`. No observation is fabricated, no deviation is
   invented, and the row is not omitted, so no generic fallback can invent a
   numerical value either.
2. **Response schema** (`backend/tos/schemas/response_v1.py`):
   `PhysicsMeasurement` gained `measured: bool = True`, a REQUIRED-but-NULLABLE
   `observed: Optional[float]` and a nullable `deviation: Optional[str]`. An
   ordinary measurement keeps its number; an unmeasured one carries null with
   the explicit flag.
3. **Adjudication** (`backend/tos/engine/adjudicator.py`): `NOT_MEASURED` is
   neither PASS nor scientific FAIL. It is excluded from `det_passed` (an
   unmeasured observation never earns a passing point), it never enters
   `failing_deterministic`, and an unmeasured REQUIRED deterministic law
   yields `INDETERMINATE` instead of a silent PASS. The existing behaviour is
   preserved: the coverage gate still yields INDETERMINATE when it applies, and
   a real deterministic FAIL/VETO is still a VETO (the veto check runs first).
4. **Consumers audited and updated**:
   - `backend/tos/forensic_artifacts/pdf_generator.py` - a `NOT_MEASURED` row
     is NOT counted as a detection (it is neither a violation nor a pass) and
     is rendered in its own "NOT MEASURED" table with the cell text "not
     measured" - no number is printed for it;
   - `backend/scripts/sentinel_audit.py` - a law with no observation is
     skipped (it was previously coerced with `float(... or 0)`, which would
     have raised on a null observation);
   - `backend/main.py` (API serialization) and the batch processor carry the
     adjudicator rows through unchanged - the Pydantic model is the contract
     boundary and it now carries the flag and the nullable observation;
   - the dashboard consumes the API JSON and does not read `observed` (no
     matching reference), so no front-end change is required by this patch.
5. **Scope control**: LAW-160's threshold (4.5 A), law ID, operator, units and
   modality classification are unchanged; the canonical law definition, the
   modality matrix, the cohorts and the calibration artefacts are NOT touched.
   The new state is limited to an explicitly unmeasured measurement and the
   contract handling it needs - the global missing-measurement fallback for a
   law the engine did not produce at all is deliberately left as it was, so no
   other law's status is reinterpreted by this patch.

### 8.2 Why this is an output-contract change (and where it is bounded)

The change is bounded to the unmeasured edge case, but it IS an API
output-contract change: a consumer that assumed `observed` is always a number
must now read the `measured` flag (or tolerate null). The response model keeps
`observed` REQUIRED, so a missing field is still a validation error - only a
genuinely unmeasured measurement may carry null. The consumers inside this
repository (the API, the batch path, the PDF dossier, the sentinel audit) were
audited and updated accordingly; the acceptance tests below prove the
serialisation and the contract behaviour.

### 8.3 Tests added (end to end)

`backend/tests/test_not_measured_contract.py` (new) and
`backend/tests/test_law_160.py` (updated):

- a structure with no eligible C-alpha pairs returns a LAW-160 row with
  `measured=false`, a null observation, a null deviation, zero sample size and
  `NOT_MEASURED`;
- the response model serialises that row successfully and still validates
  ordinary measured rows with their numbers;
- a structure with that unmeasured LAW-160 row and otherwise PASSing
  deterministic laws yields `INDETERMINATE`, not PASS;
- a genuine deterministic veto (LAW-120) remains a VETO when LAW-160 is
  `NOT_MEASURED`;
- `NOT_MEASURED` is absent from the failing deterministic list and is never
  counted as a pass;
- the coverage gate keeps its INDETERMINATE behaviour;
- the existing insertion-code and gap regressions remain green.

### 8.4 Results on this candidate

```
$ cd backend && python -m pytest tests/test_law_160.py tests/test_not_measured_contract.py \
      tests/test_adjudicator.py tests/test_response_schema.py -q --no-header
36 passed, 5 skipped
```
(the 5 skips are the pre-existing response-schema tests that need a captured
live sample at /tmp - unrelated to this change)

```
$ cd backend && python -m pytest tests -q --no-header <the five POSIX-only modules ignored>
25 failed, 296 passed, 5 skipped, 1 warning, 23 errors
```
The failing/erroring set is the SAME pre-existing set as the unmodified
baseline (25 failed / 284 passed / 23 errors on the base tree); the delta is
the new tests. The pre-existing failures are Windows import/API issues in
other modules and are disclosed, not suppressed.

```
canon  : 6a9cd4b4349b81de          (unchanged; the Stability Gate lock)
matrix : b1249cedde5018628ea2351de49c9622cdce4657561af7f74146df182d1a210d
```
The local matrix value is the CRLF-checkout artefact described in 3.3 (the
Linux runner computes the asserted value on the same tree; the Stability Gate
re-verified both locks on the candidate).

```
patch secret-pattern hits: 0
```

### 8.5 Status of the LAW-160 task

The zero-sample representation is now CLOSED: an unmeasured LAW-160 is
represented explicitly and truthfully end to end (engine, adjudicator, API
schema, PDF dossier, sentinel audit) and never produces a misleading
scientific result. No further LAW-160 item is outstanding. The next work
items, each on its own branch, are LAW-130 (fail-open clash measurement and the
`<` vs the canonical `<=` at the 20.0 boundary) and then the calibration
runner's authentication and atomic publication.

### 8.6 A bookkeeping error in this branch's own push log (disclosed and repaired)

The confirmation of the `1831bdf` push was written with the wrong identifier
base (`P-084-CONFIRM` instead of `P-083-CONFIRM`), so the identifier-pairing
rule rejected it and the governance acceptance suite (AT-20/AT-24) went red on
this branch. The repair is appended, not rewritten: `P-084-INTENT` re-establishes
the P-084 pair for the SAME push and declares no new push, `P-086` discharges
the orphaned `P-083-INTENT`, and the published entries are untouched (the
append-only parent-prefix check passes). The records are ASCII on purpose: the
governance tests print the log, and a non-ASCII marker breaks the Windows
console encoding.

## 9. Round-25 continuation: the final candidate and the CI evidence

**Final candidate SHA:** `9a149aae2dacd06e2e24f4ad699620311a906614` (branch `feature/law-160-sequence-adjacency`; base `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`).

**CI on the exact candidate 9a149aae2dacd06e2e24f4ad699620311a906614:**

| Run | Event | Conclusion |
|---|---|---|
| `39649469536` (Toscanini Stability Gate) | push | **success** - "Unit Tests + Hash Lock" and "Secret Scan (gitleaks)" both green on this exact candidate |
| `39649482980` (Toscanini Stability Gate) | pull request | **success** |
| `39649469505` (Governance Acceptance) | push | failure - ONE job: `EC-18/EC-24 push-log schema, append-only, remote verification` |
| `39649483107` (Governance Acceptance) | pull request | failure - the same push-log job plus `EC-16/19/20 PR-body claim verification` |

**Local verification on the same tree (Windows host; the Linux runner is the
supported platform and its Stability Gate is green):**

```
$ cd backend && python -m pytest tests -q --no-header <the five POSIX-only modules ignored>
25 failed, 302 passed, 5 skipped, 445 warnings, 23 errors
```
The failing/erroring set is the SAME pre-existing set as the unmodified
baseline (25 failed / 284 passed / 23 errors on the base tree); the delta is
the new tests. The pre-existing failures are Windows import/API issues in
other modules and are disclosed, not suppressed.

```
canon  : 6a9cd4b4349b81de          (unchanged)
matrix : b1249cedde5018628ea2351de49c9622cdce4657561af7f74146df182d1a210d
```
(the local matrix value is the CRLF-checkout artefact of 3.3; the Linux
runner's Stability Gate re-verified BOTH locks on this candidate)

```
patch secret-pattern hits: 0
```

**Remaining blockers (reported honestly):**

1. The push-log remote verification is red because the push log inherited
   from `main` records confirmations for `repair/round19-mainline-acceptance`
   and that branch was deleted by the repository-wide pruning; the historical
   remote references are truthfully unverifiable. Per the auditor's direction
   this is NOT repaired in this engine task and no gate is bypassed.
2. The PR-body claim verification is red because the claim gate requires a
   successful run of BOTH required workflows on the exact head, and the
   Governance Acceptance workflow is red for the reason in (1). The PR body
   names the exact collectable pytest node IDs and cites the successful
   Stability Gate run on this exact head; the Governance Acceptance run is
   disclosed as red rather than cited as a pass.
3. The candidate is therefore NOT mergeable yet: the technical defects of
   Round 25 are closed, but the governance evidence is not green, no
   independent approval exists on this head, and the standing hard stop
   forbids a merge.
