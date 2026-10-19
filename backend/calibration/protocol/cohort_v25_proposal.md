# Cohort v25 rebuild proposal (P0 item 16 — numeric criteria are OWNER decisions)

**Status:** PROPOSAL. Nothing in this document is executed or frozen. Per
the audit direction: "The owner sets all numeric criteria." The held-out
split is NOT written here — it may only be written after the owner
approves the cohort definition and BEFORE any threshold work
(protocol §2, §7).

## 1. Removals from the v23 cohort (verified, recorded)

### 1I0T — REMOVE (two independent grounds)

1. **It is not a protein.** 1I0T is a Z-DNA duplex (DC/DG nucleotides +
   spermine + water; 240 atoms, no amino-acid backbone). The v23 dataset
   labels it "Ferredoxin". The corrected engine records NOT_MEASURED for
   every backbone law and VETO on LAW-170 (nucleotides are non-standard
   residues) — verified in the v24 run record.
2. **Resolution provenance is inconsistent.** The committed
   `dataset_v23.py` records 0.89 Å; the v24 run recorded 0.6 Å for the
   same target. One of the two is wrong, and the v24 run's provenance is
   already superseded (item 3). The v25 cohort must re-derive every
   annotation from wwPDB/RCSB metadata rather than copying either value.

### AF-P0DTC2-F1 — REPLACE OR FOLLOW MIGRATED ENTITY (owner decision)

Live-verified 2026-10-19 (evidence/item-11/): the AFDB API no longer
serves `AF-P0DTC2-F1`; the entry now answers as `AF-0000000365840314`
(provider VR3D, ColabFold v1.5.2, globalMetric 67.14). Options:
(a) keep the biological target and follow the migrated entity ID (then
the cohort records the mapping), or (b) replace the target with another
membrane/protein target. The UniProt entry itself is still active and
reviewed.

### AF-P62988-F1 — REMOVE (ground recorded)

Live-verified 2026-10-19: the AFDB API returns HTTP 404 for P62988 and
UniProt records `entryType: Inactive` (merged entry — UBIQ_HUMAN was
merged into another ubiquitin entry). No prediction can be acquired;
removal from the predicted stratum is justified by the upstream merge.

## 2. Proposed quality filter (PISCES-style; ALL numbers are [owner])

Source pool: wwPDB/RCSB X-ray entries, then:

| Criterion | Proposed default | [owner] |
|---|---|---|
| Resolution | ≤ 2.0 Å | to set |
| R-free | ≤ 0.25 | to set |
| Completeness | ≥ 90 % | to set |
| wwPDB validation-report anomalies | no serious entry-level covalent-geometry outliers; clashscore ≤ 5 (protocol §2 Cohort A gates) | to set |
| Sequence identity (PISCES-style culling) | ≤ 30 % pairwise identity | to set |
| Chain count | single-chain or per-chain evaluation | to set |
| Size window | [owner] | to set |
| Target count | 200–400 structures total (per prior analysis; protocol suggests ≥ 15 per stratum minimum) | to set |

Strata (kept separate, never pooled): X-ray (<1.5 / 1.5–2.5 / >2.5 Å),
cryo-EM, NMR (first model), predicted (AFDB, pinned versions, reported
separately per protocol).

## 3. Negative controls (Cohort B, per protocol)

Structural defect injection — Cα displacement, induced clashes, chain
breaks, glycine substitution in helices — as fresh PDB bytes with
recorded digests. None exist today (`control/negative_controls.py` is an
unrelated MI scrambled-null utility). Implementation starts only after
the owner approves the cohort plan (it is part of the same decision).

## 4. Order of operations (binding)

1. Owner sets every numeric criterion above and approves the target
   counts per stratum.
2. Cohort built and digest-pinned (RCSB/AFDB downloads, wwPDB
   validation-report gates applied).
3. **Held-out split written and committed BEFORE any measurement or
   threshold work** (protocol §2 Cohort C: 20 % held out, disjoint by
   structure identity from {v23-cohort, calibration-v24, validation}).
4. v25 diagnostic run under the hardened runner (50 of 50 or the
   approved count; executing-commit provenance; append-only run dir).
5. Threshold derivation per protocol §5-§6 on the training split only;
   held-out inspected exactly once after criteria are frozen.

## 5. What the coder needs from the owner

- Approve or amend every [owner] number in §2.
- Decide the P0DTC2 disposition (§1).
- Approve the negative-control defect set (§3).
- Confirm the split fraction (protocol default 20 %) and the exclusion
  sets for held-out disjointness.
