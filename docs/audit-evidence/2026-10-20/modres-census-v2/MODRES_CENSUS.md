# MODRES Census v2 — v23 40-Target Experimental Cohort

**Recorded as documentation.** No cohort, code, threshold, or calibration
change follows from this census. This is the v2 (Round-13) portability
fix of the Round-11 census: the input artifact is now repository-relative
and digest-recorded, and the coordinate bytes are fetched from the
canonical public source with per-target digest verification — the bundle
no longer depends on any coder's local filesystem.

## Question

Do any of the 40 experimental coordinate files of the recorded v23
cohort carry modified residues?

## Method and provenance binding (v2)

- **Input artifact (repo-relative):**
  `docs/audit-evidence/2026-10-20/modres-census-v2/n7_od16_40target_table.json`
  — the Round-9/10 item-level table, extracted byte-for-byte from the
  sealed r2 evidence archive
  (`outputs/control-atlas-evidence-2026-10-20-round9-r2.tar.gz`, member
  `./02-n7-od16-cohort-audit/n7_od16_40target_table.json`, archive
  SHA-256 `5ac52b9332c2fe2c2510b443c6aeb26ce785eccf61ad6d44cebce6151a3257cd`).
  Input artifact SHA-256 (recorded in the results JSON):
  `4190a774fb46e8546c4b82d4fda700677991cdf9f864aaab165c8ea796994196`.
  The table records `dataset_v23_blob_sha1 =
  3a5358374b2037b9475bd2d494754944485da5ac`, which MATCHES the pinned
  `backend/calibration/dataset_v23.py` blob on this branch.
- **Coordinate bytes:** fetched from the canonical source
  `https://files.rcsb.org/download/<pdb_id>.pdb`; each file's SHA-256 is
  verified against the table's `c5_coordinate_sha256_fresh` BEFORE
  counting. 40/40 bound; 0 digest failures.
- **Counted:** PDB `MODRES ` records; and, for context, HETATM atoms
  whose residue name is in the engine's `MODIFIED_RESIDUE_PARENTS`
  mapping (MSE, SEC, PYL, HYP, CME, CSO, CSD, OCS, CSX, CAS, CSM, MLY,
  M3L, KCX, SEP, TPO, PTR, PCA, NLE).
- **Script:** `modres_census.py` (this directory); run log
  `census_run.log`; machine-readable results
  `modres_census_results.json` (per-target, with per-target source URL,
  SHA-256, digest binding, and counts). No absolute local paths appear
  in any output.

## Results (identical to the Round-11 census — reproduced, not assumed)

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
   by this cohort.

Both statements are recorded side by side because they answer slightly
different questions; the operative one is what matters for calibration
semantics. This v2 re-run on canonical RCSB bytes with recorded input
digests reproduces the Round-11 result exactly (11 / 0 / 40/40 bound),
which is the portability demonstration the Round-12 direction required.
