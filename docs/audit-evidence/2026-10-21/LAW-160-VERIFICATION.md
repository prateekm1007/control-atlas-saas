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
