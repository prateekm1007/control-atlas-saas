# Provenance v2 evidence — feature/provenance-v2

Base: science/p0-correctness-2026-10-19 @ ebc1fd79f20b9a6b6a262a70ce6942ba069d19e9
Code commit: 4bddfe7cec42e7032731e992cbdf01d34821f27f (this commit adds
evidence artifacts only; the demo manifest below was generated from that
exact code state, as recorded inside it).

## Demonstration manifest (record-only, explicitly not a calibration run)

    cd backend && python3 -m calibration.runner --emit-provenance \
        ../docs/audit-evidence/2026-10-20/provenance-v2/runner_provenance_demo.json

Output:

    {"written": ".../runner_provenance_demo.json",
     "schema": "toscanini-runner-provenance/2.0.0",
     "code_sha": "4bddfe7cec42e7032731e992cbdf01d34821f27f",
     "working_tree_clean": true}

Contents (abridged):

    schema:               toscanini-runner-provenance/2.0.0
    calibration_executed: False
    record_type:          provenance-record-only
    engine.source_provenance: git-checkout
    engine.code_sha:      4bddfe7cec42e7032731e992cbdf01d34821f27f
    git_context.remotes:  {"origin": "https://github.com/prateekm1007/control-atlas-saas.git"}
    dependencies:         {"numpy": "2.1.3", "gemmi": "0.7.5"}

The runner identifies itself as v25-runner/1.0.0, and the record-only
nature is explicit (`calibration_executed: false`), so this demonstration
cannot be mistaken for a v25 calibration run.

## Round-12 direction 4 compliance matrix

| Required fix | Implementation | Test |
|---|---|---|
| Package-install audit: must not report an unrelated parent repository as the engine source | `_engine_repo_root()` requires work-tree containment AND `git ls-files --error-unmatch` tracking | `test_unrelated_parent_repo_is_rejected` (also asserts the unrelated HEAD never appears in the record) |
| Installed-artifact provenance route, or explicit unavailability, never inference | `_installed_distribution()` (name/version/sanitized direct_url); `source_provenance: "installed-package"` with `code_sha: null`; `"unavailable"` route with reason | `test_installed_package_route_records_distribution`, `test_no_git_metadata_records_unavailable` |
| Fail closed for authoritative measurement | `_engine_provenance(strict=True)` (default; used by `run_calibration`) raises SystemExit off the verified git-checkout route | `test_strict_engine_provenance_fails_closed`, `test_installed_package_route_records_distribution` (strict refusal) |
| Tests: normal checkout; clean/dirty tree; installed package; no Git metadata | all four states covered; dirty trees are recorded, not refused, by the record path | `test_engine_repo_root_in_normal_checkout`, `test_dirty_tree_recorded_not_refused`, installed/no-Git tests above |
| Credential sanitization in HTTPS and SCP-style remotes; no secrets in logs, manifests, or exception text | `_sanitize_git_url` (unchanged behavior, retested); `direct_url.json` URLs sanitized on read; strict-exit text asserted credential-free | `test_sanitize_strips_userinfo`, `test_git_context_sanitizes_poisoned_remotes`, `test_strict_engine_provenance_fails_closed`, `test_emit_provenance_cli_writes_manifest` |
| Demonstration manifest records `calibration_executed=false` | top-level `calibration_executed` + `record_type` fields; run manifests record `True` | `test_manifest_has_required_blocks`, `test_emit_provenance_cli_writes_manifest`, `test_runner_hardening` run-manifest assertions |

## Commands run and output (Python 3.12.15; gemmi 0.7.5)

1. python3 -m pytest tests/test_runner_provenance.py tests/test_runner_hardening.py tests/test_calibration_runner.py -q
   -> 24 passed in 4.98s
2. python3 -m pytest tests -q
   -> 4 failed, 515 passed, 5 skipped in 12.94s
   Baseline (pristine ebc1fd79, same interpreter): 4 failed, 502 passed,
   5 skipped. Delta = +13 (new provenance tests). The 4 failures are the
   pre-existing environmental TestUsageWidget cases (streamlit absent
   locally; green in CI).
3. The demonstration manifest command above, exit code 0.

## Verdict impact

NONE. Calibration-runner provenance only; no law, threshold, operator,
scope, or adjudication semantics touched; no calibration run executed
(hard stop holds; OD-13 unchanged).
