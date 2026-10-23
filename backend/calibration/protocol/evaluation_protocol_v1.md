# Independent Evaluation Protocol v1 (WP-SCI-1)

**Status:** ACTIVE protocol definition — Predetermined before any held-out
results are inspected.
**Effective:** 2026-10-19
**Scope:** Independent validation of the corrected LAW-120 (Bond Angle,
definition v2) and of the overall deterministic verdict path, per the
2026-10-19 scientific-validity directive §5 and the SaaS Engineering
Constitution §4 (Proposed, PR #11).

---

## 0. What this protocol is and is not

- This protocol defines how Toscanini's deterministic verdicts are to be
  **independently validated** before any performance claim is made.
- The v23 cohort (50 structures) is **diagnostic data only**. It was used
  to select the 18.0° LAW-120 threshold under the *uncorrected* metric;
  it cannot establish validity for the corrected metric and is not a
  held-out set.
- No threshold may be tuned, adjusted, or "confirmed" against the held-out
  set. The held-out set is inspected exactly once, after acceptance
  criteria are frozen.
- Until this protocol completes with its acceptance criteria met, LAW-120
  remains **contained** (`station_sop.LAW_CONTAINMENT`): a measured LAW-120
  violation cannot independently issue an authoritative VETO.

## 1. Correction record: the "at most 20 of 40" statement

An earlier auditor statement asserted that *at most 20 of 40 experimental
structures could pass*. **Independent recalculation from the recorded v23
rows** (backend/calibration/calibration_results.json, 50 targets) with the
current modality matrix applied to the recorded per-law statuses gives:

- **0 / 40** experimental structures could pass — 38 blocked by LAW-120 at
  its recorded 10.0° threshold, 2 by LAW-170 (1I0T, 1DDL).
- The recorded verdict column itself shows 44 VETO / 5 INDETERMINATE /
  1 PASS across all 50 targets.

The "at most 20 of 40" number is **not reproducible** from the saved rows
and current matrix rules and is hereby corrected to 0/40 for the recorded
run. Any what-if simulation re-using the v23 cohort must be labelled as
*simulation on diagnostic data*, never as validated performance.

## 2. Cohorts (no overlap by structure identity)

### Cohort A — acceptable structures (independently qualified)

- Source: PDB entries with a wwPDB validation report (same coordinate
  files, downloaded and digest-pinned like calibration inputs).
- Inclusion gates, defined from wwPDB validation data **only** (never from
  Toscanini output):
  - zero Ramachandran outliers (MolProbity), clashscore ≤ 5,
    zero covalent-geometry serious outliers at entry level;
  - resolution known; deposition released ≥ 2016 (validation-report era).
- Strata: X-ray (<1.5 Å, 1.5–2.5 Å, >2.5 Å), cryo-EM, NMR (first model),
  at least 15 structures per stratum where availability allows.
- Predicted structures (AlphaFold DB, pinned version) form a separate
  predicted stratum with label "acceptable-by-source-quality" only for
  high-pLDDT regions; they are reported separately and never pooled into
  experimental rates.

### Cohort B — defective structures (labels independent of Toscanini)

- **B1 — perturbed-by-construction:** deterministic perturbation operators
  applied to Cohort A seed coordinates, with the perturbation recipe and
  magnitude recorded per file. Label "defective" is true **by
  construction**, independent of any engine verdict. Operators (v1):
  - `bond-stretch`: N/CA/C bonds ±0.15–0.30 Å on a documented residue set
  - `angle-distort`: backbone angles distorted by +25° … +60°
    (LAW-120-relevant)
  - `chain-break`: Cα–Cα of adjacent-numbered residues displaced to 8–12 Å
  - `clash-insert`: inter-residue heavy-atom overlaps < 1.5 Å
  - `peptide-flip`: omega driven to ~90° (planarity violation)
- **B2 — known-defective deposits:** entries flagged by wwPDB with serious
  geometry outliers or retracted/corrected status, if available.
- Every B file carries a `defect_manifest.json` (operator, parameters,
  seed residues, expected affected laws) — created by the perturbation
  tooling, not by the engine.

### Cohort C — held-out set

- Drawn from A/B1/B2 candidate pools **before any engine run**, frozen by
  committing `split_manifest.json` (structure IDs + digests + strata).
- Size: ≥ 60 structures (≥ 30 acceptable, ≥ 30 defective) or 20% of the
  pooled candidates, whichever is larger.
- Used exactly once, after calibration thresholds are frozen and accepted.

## 3. Split rules

- No structure identity appears in more than one of
  {v23-cohort, calibration-v24, validation, held-out} (by PDB/AF ID).
- **Near-identity guard:** structures with ≥ 90% sequence identity
  (UniProt-clustered or pairwise alignment over the same chain) or
  identical sequence + same method must land in the same split or be
  dropped. The split builder MUST run this check and record it.
- Strata proportions are mirrored across splits (documented in the split
  manifest).

## 4. Comparison metrics (same coordinate inputs)

For every structure in A/B/C, pull the wwPDB validation XML (or
MolProbity-equivalent metrics) for the SAME coordinate digest:

- MolProbity clashscore, Ramachandran outlier %, rotamer outlier %,
  RSRZ / R-free (X-ray), EM fit metrics (cryo-EM).
- Report per-law correlation of Toscanini's observed values with the
  corresponding wwPDB metrics (LAW-125 vs Ramachandran %, LAW-130 vs
  clashscore, LAW-150 vs rotamer %). Disagreements are analysed, not hidden.

## 5. Metrics computed (fixed before held-out inspection)

From `backend/calibration/evaluate.py` (labels vs verdicts; defective =
positive class):

- sensitivity (VETO on defective), specificity (PASS on acceptable),
  false-veto rate (VETO on acceptable), false-accept rate (PASS on
  defective) — each with 95% Wilson intervals;
- INDETERMINATE rate by cohort (never counted toward sensitivity or
  specificity);
- per-law blocking contribution, per-law NOT_MEASURED / NOT_APPLICABLE
  rates;
- performance by modality/resolution stratum.

## 6. Acceptance criteria (predetermined)

LAW-120 containment may be lifted, and LAW-120 may be described as
validated, only if ALL of the following hold on the **validation split**,
and are then **confirmed on the held-out set**:

1. False-veto rate on Cohort A ≤ 10% with 95% Wilson upper bound ≤ 20%.
2. Sensitivity on B1 `angle-distort` ≥ 80% (detection of the injected
   angle defects); sensitivity on other B1 operators reported (not gated).
3. False-accept rate on B1 `angle-distort` ≤ 10%.
4. LAW-120 NOT_MEASURED rate on Cohort A ≤ 5%.
5. No modality stratum inverts the direction of the effect (high-res
   X-ray false-veto ≤ low-res X-ray false-veto + 5 points).
6. Threshold selection is documented as derived ONLY from calibration/
   validation data, with the selection run recorded and digest-bound.

## 7. Kill criteria (any one halts the claim)

1. False-veto upper bound > 20% on Cohort A.
2. Sensitivity < 60% on B1 `angle-distort`.
3. NOT_MEASURED rate > 10% on any qualified cohort stratum.
4. Evidence that labels leaked across splits (near-identity guard failure
   discovered post hoc).
5. Threshold tuning against the held-out set (automatic kill — the split
   is voided and a new one must be drawn).

## 8. Containment lifting procedure

1. Calibration/validation run artifacts committed under
   `backend/calibration/runs/` with manifests.
2. A reviewed policy PR updates `station_sop.LAW_CONTAINMENT` (removal),
   `governance/law_registry.yaml` (containment block removal, validated
   threshold citation), and re-pins `governance/law_semantic.sha256` —
   CI enforces the explicit re-pin.
3. README/SAAS_HANDOFF claims updated in the same PR.

## 9. Operator prerequisites (recorded, not assumed)

- Network egress to files.rcsb.org (verified in the v24 run) and to the
  wwPDB validation service for labels.
- AlphaFold DB bulk file acquisition was **blocked** from the WP-SCI-1
  sandbox (HTTP 404/403 on all known bucket paths); predicted-stratum
  cohort construction requires an approved egress path or a pinned local
  mirror. This is recorded in the v24 run (10/50 acquisition failures) and
  is a prerequisite, not a result.
- Design-partner interviews / pilot recruitment are out of scope for this
  protocol: no interview, usage, or acceptance claim may be made without
  the actual recorded evidence.
