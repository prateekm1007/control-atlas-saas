# MODRES Census — v23 40-Target Experimental Cohort (U9, Round-11)

**Recorded as documentation.** No cohort, code, threshold, or calibration
change follows from this census (U9).

## Question

Do any of the 40 experimental coordinate files of the recorded v23
cohort carry modified residues?

## Method and provenance binding

- Cohort and per-file digest binding: the Round-9/10 item-level table
  (`n7_od16_40target_table.json`, which itself records
  `dataset_v23_blob_sha1 = 3a5358374b2037b9475bd2d494754944485da5ac`,
  verified to MATCH the pinned `backend/calibration/dataset_v23.py`
  blob on `science/p0-correctness-2026-10-19` @ `ebc1fd79…`).
- Each local coordinate file's SHA-256 was verified against the table's
  `c5_coordinate_sha256_fresh` **before** counting. 40/40 bound; 0 digest
  failures.
- Counted: PDB `MODRES ` records; and, for context, HETATM atoms whose
  residue name is in the engine's `MODIFIED_RESIDUE_PARENTS` mapping
  (MSE, SEC, PYL, HYP, CME, CSO, CSD, OCS, CSX, CAS, CSM, MLY, M3L,
  KCX, SEP, TPO, PTR, PCA, NLE).
- Script: `evidence/modres-census/` (copy of
  `scripts/modres_census.py`), run log `census_run.log`, machine-readable
  results `modres_census_results.json` (per-target, with digests).

## Results

| Measure | Value |
|---|---|
| Targets censused | 40 / 40, all digest-bound to the recorded table |
| Files with any `MODRES` record | **1 of 40** (`5A63`) |
| Total `MODRES` records | **11** (all in `5A63`) |
| HETATM atoms of modified-residue names (engine mapping) | **0 across all 40 files** |

**The one non-zero file, characterized.** All 11 `MODRES` records in
`5A63` read `ASN … GLYCOSYLATION SITE` (chain A, residues 45, 55, 187,
264, 387, 435, 464, 506, 530, 562, 573). These are *site annotations*:
each named residue is deposited as a **standard ASN in ATOM records**
(5A63 has 9,646 ATOM records; ASN appears as a normal residue), and the
actual carbohydrate moieties (NAG, BMA) are separate HETATM records
(357 HETATM records total), which the engine already excludes as
non-polymer heteroatoms with no amino-acid parent.

## Conclusion — both readings recorded, neither smoothed

1. **Literal PDB-record reading:** the cohort is NOT annotation-free —
   `5A63` carries 11 MODRES annotation records. The statement "zero
   MODRES records in the v23 cohort" would be false as a bare claim.
2. **Operative coordinate reading (the engine's sense):** zero modified
   residues are present in the coordinates of the recorded cohort — no
   MSE/SEP/TPO-class residue replaces a standard one anywhere; the
   engine's `MODIFIED_RESIDUE_PARENTS` retention path is not exercised
   by this cohort. The auditor's "zero modified residues in the v23
   cohort" is correct in this sense.

Both statements are recorded side by side because they answer slightly
different questions; the operative one is what matters for calibration
semantics. Note also that `5A63` is among the identity-error rows of the
quarantined v23 metadata (recorded as TRPV1; actually gamma-secretase —
see the item-level table), which is consistent with a cryo-EM entry
carrying glycosylation-site annotations.
