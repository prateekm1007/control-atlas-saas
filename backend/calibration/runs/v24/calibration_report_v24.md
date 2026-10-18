# TOSCANINI v24 Diagnostic Calibration Report (WP-SCI-1)

**Run type: DIAGNOSTIC — NOT a threshold validation.**
**Date:** 2026-10-19 · **Runner:** `v24-runner/1.0.0` · **Completeness: 40/50 measured, 10/50 acquisition failures (AlphaFold DB unreachable from this environment) → run status INCOMPLETE (runner exits non-zero by design).**

Evidence binding: `run_manifest.json` (code SHA, canon hash
`6a9cd4b4349b81de`, law-semantic hash, matrix hash, registry hash,
dataset digest), `coordinate_manifest.json` (40 SHA-256-pinned coordinate
files), `MANIFEST.sha256`.

---

## 1. Why this run exists

The v23 recorded run is stale and scientifically unusable as validation:

- it recorded LAW-120 thresholds of **10.0°** while the live canon said 18.0°;
- it contains **impossible rows** (1I0T and AF-P62988-F1 with
  `sample_size = 0`, `observed = 0.0`, `status = PASS` for LAW-120);
- its LAW-120 decision (10°→18°) was selected from the same cohort it
  claims to validate (circular), using a **measurement-definition bug**
  (below).

The LAW-120 measurement definition was repaired in WP-SCI-1: the peptide
carbonyl angle is now **O(i)–C(i)–N(i+1)** over valid sequential pairs
(previously O(i)–C(i)–**N(i)** — the same residue's amide nitrogen, which
is not bonded to C(i), compared against the Engh-Huber peptide ideal
123.0°). LAW-160 now requires true sequence adjacency, and every law now
distinguishes NOT_MEASURED / NOT_APPLICABLE / MEASUREMENT_FAILED from
measured PASS / VETO / FAIL.

This run documents the **corrected metrics** on the v23 cohort (diagnostic
reuse only). Per `protocol/evaluation_protocol_v1.md`, **no threshold is
selected, confirmed, or re-validated from this run**, and LAW-120 remains
**contained** (a measured LAW-120 violation cannot independently issue an
authoritative VETO).

## 2. Completeness and acquisition failures (honest representation)

- 40/50 targets acquired from files.rcsb.org and measured; every
  coordinate file pinned by SHA-256 in `coordinate_manifest.json`.
- **10/10 AlphaFold predicted targets FAILED acquisition**: the AlphaFold
  DB file bucket returns HTTP 404 (NoSuchKey) for all known version
  paths, and 403 without a UA; EBI FTP-mirror and GCS public-dataset
  paths also 404. These are recorded as explicit acquisition failures —
  **not** counted as calibration cases, not substituted with synthetic
  coordinates, and the runner exits non-zero.
- Operator prerequisite recorded in the protocol (§9): predicted-stratum
  work requires an approved egress path or a pinned local mirror.

## 3. Corrected-metric results (40 measured experimental structures)

Top-level verdicts: **33 PASS / 6 VETO / 1 INDETERMINATE** (v23 recorded:
44 VETO / 5 INDETERMINATE / 1 PASS under the uncorrected metric).

| Law | Result (40 measured) | Reading |
|---|---|---|
| LAW-120 (corrected) | 39 measured PASS, 1 NOT_MEASURED, **0 VETO**; observed 0.51–7.04°, median 2.03°, max sample 24 003 angles | Real bond-angle deviations of well-validated structures are ~0.5–7° RMSD. The uncorrected metric's 13–18° "RMSD" was an artifact of the pseudo-angle. No structure approaches the pinned 18.0° threshold, which therefore remains **unvalidated and non-binding (contained)** — plausible, but its acceptance awaits the protocol. |
| LAW-160 (adjacency) | 38 measured PASS, 1 NOT_MEASURED, 1 advisory FAIL; **no deterministic breaches** | All 7 v23 "breaches" (5A63 29.68 Å, 6NB6 24.91, 5T4O 14.08, 2ACE 9.82, 1PPE 7.12, 7JTL 6.39, 1BRS 4.95) were **numbering-gap artifacts**. With true adjacency, only 1PPE shows a genuine 7.12 Å Cα–Cα discontinuity between adjacent-numbered residues — a real coordinate break, honestly measured and (per the existing matrix) advisory for experimental structures. |
| LAW-195 | 9 PASS (structures with disulfides), **30 NOT_APPLICABLE**, 1 measured VETO (5A63) | Vacuous PASSes are gone; the no-pair case is recorded as not applicable. |
| All others | NOT_MEASURED appears on 9 law rows total | Honest empty samples (e.g. single-strand or terminal cases), never zeros. |

### VETO causes (measured, per current policy)

LAW-170: 1I0T, 1DDL · LAW-150 (deterministic for X-ray/cryo-EM): 1HHO,
3PGK, 5A63 · LAW-125: 3PGK · LAW-110: 1PPE · LAW-195: 5A63.
5T4O is INDETERMINATE (LAW-150 missing measurement). LAW-182 (heuristic)
flags 3PGK and 5A63 but cannot veto by design.

### Dataset curation finding: 1I0T is not a protein

1I0T (v23: "Ferredoxin", X-ray 0.89 Å) is a **Z-DNA duplex**
(DC/DG nucleotides + spermine + water; 240 atoms, no amino-acid
backbone). The corrected engine records NOT_MEASURED for every backbone
law — which is exactly why v23's zero-sample `PASS` row was an
impossible value rather than a real measurement. **1I0T must be replaced
in the next cohort revision** (protocol §2); the v24 run keeps it with
honest statuses so the curation error stays auditable.

## 4. The "20 of 40" auditor claim — independent correction

Applying the current modality matrix to the recorded v23 per-law rows:
**0/40** experimental structures could pass (38 blocked by LAW-120@10.0°,
2 by LAW-170). The "at most 20 of 40" statement is not reproducible from
the saved rows and is corrected to **0/40** for the recorded run. See
protocol §1.

## 5. What this run does NOT establish

- No LAW-120 threshold is validated. The 18.0° canon value is retained
  for byte-stability of the protected canon hash and is contained.
- No performance claim (sensitivity/specificity/false-veto) is made —
  those require the labeled cohorts and held-out set of
  `protocol/evaluation_protocol_v1.md`, which have not been run.
- Predicted-structure behavior of the corrected metric is unmeasured
  (acquisition blocked).

## 6. Reproduction

```bash
cd backend
python -m calibration.runner --dataset v23 --out calibration/runs/v24 \
    --verify-digests-from calibration/runs/v24/coordinate_manifest.json
sha256sum -c calibration/runs/v24/MANIFEST.sha256
```

The run exits non-zero while any target is unacquired; that is the
designed fail-closed behavior, not a defect.
