# Calibration — canonical entry points and history

## Canonical current-version entry point

- **Runner:** `backend/calibration/runner.py`
  (`cd backend && python -m calibration.runner --dataset v23 --out calibration/runs/v24`)
- **Current diagnostic run:** `backend/calibration/runs/v24/`
  (`run_manifest.json` binds code SHA, policy hashes, dataset digest;
  `coordinate_manifest.json` pins every acquired coordinate file by
  SHA-256; `MANIFEST.sha256` covers the run directory).
- **Evaluation protocol:** `backend/calibration/protocol/evaluation_protocol_v1.md`
  with `split_manifest_template.json`.
- **Confusion-matrix tooling:** `backend/calibration/evaluate.py`
  (`python -m calibration.evaluate <labeled_cases.json>`).

The runner measures in-process through the production path components
(IngestionProcessor → Tier1Measurements.run_full_audit → adjudicate_laws)
and performs no authenticated network calls (no API key involved; open
mode is never enabled). It exits non-zero unless every expected target is
acquired and measured.

## Historical files (retained unmodified — do not reinterpret)

| File | What it is |
|---|---|
| `backend/calibration/calibration_results.json` | **v23 run** (2026-02-20, backend 22.5.3): 50 targets, per-law measurement schema. Recorded with the UNCORRECTED LAW-120 metric (thresholds recorded as 10.0°; contains zero-sample LAW-120 "PASS" rows for 1I0T and AF-P62988-F1). Diagnostic only. |
| `calibration_results.json` (repository root) | **Legacy v22.4.3 screen** (2026-02-15): 18 structures (9 crystal + 9 AlphaFold), coarse per-structure schema. Predates the v23 per-law format. |
| `backend/calibration/calibration_report_v23.md` | v23 narrative: raised LAW-120 10.0°→18.0° and reclassified LAW-150 (NMR) / LAW-160 (experimental) from the SAME dataset it now claims to validate — circular by the 2026-10-19 directive's standard; superseded by `protocol/evaluation_protocol_v1.md`. |
| `backend/calibration/dataset_v23.py` | v23 dataset definition (reused by the v24 diagnostic run). |

These files are kept exactly as recorded: historical results stay bound to
the code and policy that produced them and are never silently reinterpreted
under corrected semantics.
