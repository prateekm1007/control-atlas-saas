# C5 — v24 verdict-change audit (OD-03/04/20 scope changes)

**Status:** diagnostic evidence for owner review. Not a threshold validation.
**Base engine:** `c78ef15f696a2b1e2658e7ac7b54461f70c6e040` (`engine/od-scope-v25`).
**New engine:** working tree after OD-03/04/20 changes (before the delivery commit).
**Coordinate source:** v24 `coordinate_manifest.json` (`9722c3e45d316ff58ed326df074eec138b81ceb520191a492088cc0ddef63762`).
**Run manifest:** v24 `run_manifest.json` (`f602351f4844531d409cfc256939508103e3c2d4f57fe21d7f77c34494aa4036`).
**Comparison rule:** identical bytes, identical thresholds, same ingestion/adjudication path; only OD-03/04/20 scope and status semantics differ.

## 1. Acquisition and binding

- v24 manifest: **50 targets**.
- v24 acquired and SHA-256 verified: **40/50**.
- v24 acquisition failures: **10/50** (all AlphaFold `model_v4` URLs returned HTTP 404).
- All 40 acquired files re-fetched from their v24 source URLs and matched the recorded byte count and SHA-256.
- No v24 coordinate bytes were fabricated or substituted for the pinned cohort.

### v24 acquisition failures (not in the pinned table)

| Target | v24 status | v24 detail |
|---|---|---|
| AF-O15111-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-P00533-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-P02768-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-P04637-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-P0DTC2-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-P10636-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-P21802-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-P62988-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-P68871-F1 | failed | HTTPError: HTTP Error 404: Not Found |
| AF-Q9Y6K9-F1 | failed | HTTPError: HTTP Error 404: Not Found |

## 2. Verdict-change table

### 2.1 v24-pinned cohort (40/50)

Old: 34 PASS / 5 VETO / 1 INDETERMINATE.
New: 34 PASS / 5 VETO / 1 INDETERMINATE.

**Verdict changes: 0/40.** All structures are experimental, so OD-03 keeps all atoms in scope and OD-04 atom coverage is 100% in both engines.

| Target | Old verdict | New verdict | Changed law rows |
|---|---|---|---|
| 1AKE | PASS | PASS | none |
| 1AON | PASS | PASS | LAW-182: VETO -> FAIL (Advisory) |
| 1BPI | PASS | PASS | none |
| 1BRS | PASS | PASS | none |
| 1D3Z | PASS | PASS | none |
| 1DDL | VETO | VETO | none |
| 1EJG | PASS | PASS | none |
| 1FIN | PASS | PASS | none |
| 1GB1 | PASS | PASS | none |
| 1GCI | PASS | PASS | none |
| 1GFL | PASS | PASS | none |
| 1HHO | VETO | VETO | none |
| 1I0T | VETO | VETO | none |
| 1L2Y | PASS | PASS | none |
| 1MBO | PASS | PASS | none |
| 1NEK | PASS | PASS | LAW-182: VETO -> FAIL (Advisory) |
| 1PPE | PASS | PASS | none |
| 1RVS | PASS | PASS | none |
| 1TEN | PASS | PASS | none |
| 2ACE | PASS | PASS | LAW-182: VETO -> FAIL (Advisory) |
| 2CBA | PASS | PASS | none |
| 2HBS | PASS | PASS | none |
| 2JOF | PASS | PASS | none |
| 2KOX | PASS | PASS | none |
| 2LYZ | PASS | PASS | none |
| 2MJB | PASS | PASS | none |
| 2N0A | PASS | PASS | LAW-182: VETO -> FAIL (Advisory), LAW-200: VETO -> FAIL (Advisory) |
| 2OV0 | PASS | PASS | none |
| 2VB1 | PASS | PASS | none |
| 3J9I | PASS | PASS | none |
| 3LZM | PASS | PASS | none |
| 3NIR | PASS | PASS | none |
| 3PGK | VETO | VETO | LAW-182: VETO -> FAIL (Advisory) |
| 5A63 | VETO | VETO | LAW-182: VETO -> FAIL (Advisory) |
| 5T4O | INDETERMINATE | INDETERMINATE | LAW-182: VETO -> FAIL (Advisory) |
| 6J6J | PASS | PASS | none |
| 6NB6 | PASS | PASS | LAW-182: VETO -> FAIL (Advisory) |
| 6WLC | PASS | PASS | LAW-182: VETO -> FAIL (Advisory) |
| 7JTL | PASS | PASS | none |
| 7K3G | PASS | PASS | none |

### 2.2 Supplementary predicted cohort (9/10 current AFDB v6; NOT v24-pinned)

The 10 v24 AlphaFold targets had no v24 bytes or SHA-256 (all 404). To test the predicted-model path, current AFDB v6 models were fetched separately via the metadata-first API. These rows are supplementary and must not be reported as v24-pinned evidence. 9/10 were acquired; `AF-P62988-F1` still returns 404.

Old: 5 PASS / 0 VETO / 4 INDETERMINATE.
New: 5 PASS / 0 VETO / 4 INDETERMINATE.

**Verdict changes: 0/9.**

| Target | Old verdict | New verdict | Old coverage | New coverage |
|---|---|---|---|---|
| AF-O15111-F1 | PASS | PASS | 82.68 | 83.92 |
| AF-P00533-F1 | PASS | PASS | 70.66 | 71.10 |
| AF-P02768-F1 | PASS | PASS | 95.73 | 95.49 |
| AF-P04637-F1 | INDETERMINATE | INDETERMINATE | 59.80 | 61.34 |
| AF-P0DTC2-F1 | INDETERMINATE | INDETERMINATE | 55.93 | 57.01 |
| AF-P10636-F1 | INDETERMINATE | INDETERMINATE | 7.52 | 8.01 |
| AF-P21802-F1 | INDETERMINATE | INDETERMINATE | 67.24 | 67.84 |
| AF-P68871-F1 | PASS | PASS | 99.32 | 99.29 |
| AF-Q9Y6K9-F1 | PASS | PASS | 75.89 | 77.32 |

### 2.3 Combined measured set

**0 verdict changes across 49 measured structures** (40 v24-pinned + 9 supplementary).
No deterministic law changed status in either cohort.

## 3. Law-status changes (verdict-neutral)

### 3.1 v24-pinned experimental cohort

| Target | Law | Old status | New status | Observed | Sample |
|---|---|---|---|---|---|
| 1AON | LAW-182 | VETO | FAIL (Advisory) | 0.088 | 3640 |
| 1NEK | LAW-182 | VETO | FAIL (Advisory) | 0.183 | 486 |
| 2ACE | LAW-182 | VETO | FAIL (Advisory) | 0.194 | 227 |
| 2N0A | LAW-182 | VETO | FAIL (Advisory) | 0.0 | 550 |
| 2N0A | LAW-200 | VETO | FAIL (Advisory) | 395.05 | 10130 |
| 3PGK | LAW-182 | VETO | FAIL (Advisory) | 0.155 | 187 |
| 5A63 | LAW-182 | VETO | FAIL (Advisory) | 0.129 | 617 |
| 5T4O | LAW-182 | VETO | FAIL (Advisory) | 0.091 | 2352 |
| 6NB6 | LAW-182 | VETO | FAIL (Advisory) | 0.02 | 1478 |
| 6WLC | LAW-182 | VETO | FAIL (Advisory) | 0.2 | 300 |

- 10 changes: LAW-182 ×9 and LAW-200 ×1, all `VETO -> FAIL (Advisory)`.
- Reason: OD-20/OD-03 make heuristic failures advisory. They were never verdict-driving in `_compute_verdict`; the displayed row status now says so.

### 3.2 Supplementary predicted cohort

| Target | Law | Old status | New status | Observed | Sample |
|---|---|---|---|---|---|
| AF-O15111-F1 | LAW-182 | VETO | FAIL (Advisory) | 0.06 | 300 |
| AF-P00533-F1 | LAW-182 | VETO | FAIL (Advisory) | 0.038 | 471 |
| AF-P02768-F1 | LAW-182 | VETO | FAIL (Advisory) | 0.126 | 247 |
| AF-P04637-F1 | LAW-182 | VETO | FAIL (Advisory) | 0.195 | 154 |
| AF-P04637-F1 | LAW-200 | VETO | PASS | 323.77 | 3060 |
| AF-P0DTC2-F1 | LAW-182 | VETO | FAIL (Advisory) | 0.046 | 521 |
| AF-P10636-F1 | LAW-182 | VETO | FAIL (Advisory) | 0.022 | 275 |
| AF-P10636-F1 | LAW-200 | VETO | FAIL (Advisory) | 771.0 | 5544 |
| AF-P21802-F1 | LAW-182 | VETO | FAIL (Advisory) | 0.038 | 339 |
| AF-Q9Y6K9-F1 | LAW-182 | VETO | FAIL (Advisory) | 0.0 | 158 |
| AF-Q9Y6K9-F1 | LAW-200 | VETO | FAIL (Advisory) | 736.24 | 3368 |

- 11 changes: LAW-182 ×8 and LAW-200 ×3.
- LAW-182: all `VETO -> FAIL (Advisory)`.
- LAW-200: two `VETO -> FAIL (Advisory)`; one `VETO -> PASS` (`AF-P04637-F1`, observed 323.77 -> 111.79 after core-only scope). LAW-200 is heuristic, so the verdict remains INDETERMINATE from the coverage gate.

## 4. Predicted-cohort measurement changes (core scope)

No measurement changed in the 40 experimental structures (all atoms remain in scope). In the 9 supplementary predicted structures, 36 law measurements changed:

| Target | Coverage old -> new | LAW-110 n old -> new | LAW-130 n old -> new | LAW-155 observed / n old -> new | LAW-200 observed / n old -> new |
|---|---|---|---|---|---|
| AF-O15111-F1 | 82.68 -> 83.92 | 744 -> 607 | 5932 -> 4978 | 6.12 / 5932 -> 6.25 / 4978 | 157.40 / 5932 -> 134.87 / 4978 |
| AF-P00533-F1 | 70.66 -> 71.10 | 1209 -> 842 | 9392 -> 6678 | 5.97 / 9392 -> 6.04 / 6678 | 222.06 / 9392 -> 144.94 / 6678 |
| AF-P02768-F1 | 95.73 -> 95.49 | 608 -> 582 | 4858 -> 4639 | 6.16 / 4858 -> 6.14 / 4639 | 83.92 / 4858 -> 80.80 / 4639 |
| AF-P04637-F1 | 59.80 -> 61.34 | 392 -> 230 | 3060 -> 1877 | 5.91 / 3060 -> 6.17 / 1877 | 323.77 / 3060 -> 111.79 / 1877 |
| AF-P0DTC2-F1 | 55.93 -> 57.01 | 1272 -> 670 | 9940 -> 5667 | 5.98 / 9940 -> 6.15 / 5667 | 234.47 / 9940 -> 297.52 / 5667 |
| AF-P10636-F1 | 7.52 -> 8.01 | 757 -> 46 | 5544 -> 444 | 5.58 / 5544 -> 6.14 / 444 | 771.00 / 5544 -> 1599.85 / 444 |
| AF-P21802-F1 | 67.24 -> 67.84 | 820 -> 538 | 6456 -> 4380 | 6.03 / 6456 -> 6.12 / 4380 | 272.26 / 6456 -> 183.02 / 4380 |
| AF-P68871-F1 | 99.32 -> 99.29 | 146 -> 145 | 1131 -> 1123 | 5.77 / 1131 -> 5.77 / 1123 | 59.55 / 1131 -> 59.98 / 1123 |
| AF-Q9Y6K9-F1 | 75.89 -> 77.32 | 418 -> 315 | 3368 -> 2604 | 6.18 / 3368 -> 6.31 / 2604 | 736.24 / 3368 -> 569.43 / 2604 |

- LAW-110/130/155/200 now use core atoms for predicted models.
- Atom coverage changed slightly because the denominator includes all selected atoms and the numerator only core atoms.
- The four predicted structures below 70% remain INDETERMINATE; the five above remain PASS. No verdict crossed the gate.

## 5. Baseline note: c78 vs v24 stored results

For the 40 v24-pinned structures, c78 reproduces 39/40 v24 verdicts. The exception is `1PPE`: v24 stored `VETO` (LAW-110/LAW-160), while c78 and the new engine both return `PASS`. This is the already-merged N8 per-residue altloc correction, not an OD-03/04/20 effect, and is recorded here to avoid attributing it to this change set.

## 6. Result

- **No verdict changed in the v24-pinned cohort.**
- **No deterministic law status changed in either cohort.**
- The observable effect is advisory status semantics (LAW-182/200) and core-scope measurement denominators/observations for predicted models.
- The 10 missing v24 AlphaFold targets remain a coverage limitation; the supplementary current-AFDB run is explicitly separated.
