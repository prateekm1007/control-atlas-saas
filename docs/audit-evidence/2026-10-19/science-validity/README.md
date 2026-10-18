# WP-SCI-1 Scientific Validity Evidence (2026-10-19)

Directive: Coder Directive — Scientific Validity Gate and End-to-End
Validation Loop. Branch: `science/law-validity-wp15` (stacked on
`integration/p0-hardening-2026-10-19`, head `0eae6d5`).

| File | Content |
|---|---|
| `preflight.json` | Machine-readable preflight generated BEFORE any edits (policy hashes, baseline SHA, inspected files, call sites, defects SCI-1..8, scope, rollback, acceptance criteria). |
| `changed_files.json` | git status manifest at evidence build time. |
| `patch_hashes.json` | SHA-256 of the tracked diff vs base `0eae6d5` + reproduce commands; new files are committed directly. |
| `independent_recompute_20_of_40.json` | Independent correction of the auditor claim: 0/40 experimental could pass under the recorded v23 rows + current matrix. |
| `law120_independent_verification.json` | Independent numpy construction proving the corrected O(i)-C(i)-N(i+1) measurement; engine agrees to <0.05°. |
| `test_summary.json` | Exact suite result bound to the final head SHA. |
| `../MANIFEST.sha256` | Regenerated to cover this subdirectory (verify with `sha256sum -c`). |

Run artifacts committed under `backend/calibration/runs/v24/` (with their
own `MANIFEST.sha256`). No secrets or credential material appear in this
evidence set; all content is derived from repository sources and public
PDB coordinate digests.
