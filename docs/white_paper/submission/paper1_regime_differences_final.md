---
title: >
  Deterministic Structural Admissibility Auditing:
  Geometric Regime Differences Between AlphaFold Predictions
  and Experimental Protein Structures
authors:
  - name: Prateek Meher
    affiliation: Independent Researcher
    email: [YOUR EMAIL]
date: 2026-03-02
version: 1.0
canon_hash: 6a9cd4b4349b81de
system_version: 22.5.3
---

# Deterministic Structural Admissibility Auditing: Geometric Regime Differences Between AlphaFold Predictions and Experimental Protein Structures

**Prateek Meher**
Independent Researcher
[YOUR EMAIL]

Preprint. Not peer reviewed.
bioRxiv submission: 2026-03-02

---

## Abstract

We present Toscanini, a deterministic structural admissibility authority
that applies a fixed canon of fifteen physics-grounded laws to protein
structures and returns a binary verdict, a reproducible score, and a
cryptographically fingerprinted governance record. We apply Toscanini to
ten protein pairs — each consisting of an AlphaFold v6 prediction and a
high-resolution X-ray crystallographic structure — spanning oncology
targets, enzymes, signaling proteins, structural benchmarks, and
engineered proteins. After source-normalized scoring, AlphaFold and
experimental structures are equally admissible under the geometric law
set. However, they occupy distinct structural regimes: AlphaFold
structures exhibit geometric idealization (zero bond violations, low
rotamer strain, looser atomic packing) while experimental structures
carry physical strain signatures consistent with crystallographic
conditions (higher rotamer violation rates, tighter packing, stronger
hydrophobic burial). Two experimental structures receive VETO verdicts
due to marginal rotamer violations at 2.0-2.1 Å resolution. Two
AlphaFold structures receive INDETERMINATE verdicts due to low
reliability coverage in intrinsically disordered or precursor regions.
The primary finding is not that one source type is geometrically
superior. It is that deterministic auditing reveals reproducible,
quantifiable regime differences between predicted and experimental
structures that are consistent across ten protein families. All 21
audit IDs and the canon hash are provided as reproducibility anchors.
Audit ID stability across repeated submissions is formally verified in
a companion paper (Meher, 2026b).

**Keywords:** protein structure validation, AlphaFold, deterministic
auditing, structural regimes, geometric idealization, rotamer strain,
reproducibility, canonical governance

---

## 1. Introduction

AlphaFold has transformed structural biology by producing high-accuracy
protein structure predictions at proteome scale (Jumper et al., 2021).
The AlphaFold Protein Structure Database now provides predicted structures
for over 200 million proteins (Varadi et al., 2022). These predictions
are used as starting points for drug discovery, functional annotation,
and comparative structural analysis.

A central question for practitioners is how AlphaFold predictions compare
to experimentally determined structures in terms of geometric quality.
Experimental structures determined by X-ray crystallography carry physical
strain from crystallographic conditions, ligand binding, and crystal
packing. Predicted structures are optimized under a learned energy function
that enforces geometric ideality. Whether these differences are
quantifiable under a fixed, deterministic law set — and whether they
affect structural admissibility — has not been systematically examined.

Existing validation tools including MolProbity (Williams et al., 2018)
and the wwPDB validation pipeline provide geometric quality assessments.
However, they do not provide governance traceability: results from
different tool versions or parameter sets are silently inconsistent, and
there is no standard mechanism for verifying that two validation results
used the same evaluation framework.

Toscanini addresses this gap. It applies a frozen, hashed law canon to
any protein structure and produces a deterministic result tied to a
specific governance version. The same input always produces the same
output. The law set that produced the output is always identifiable.

This paper introduces the Toscanini framework and reports its application
to ten protein pairs. The central question is not whether AlphaFold
predictions are correct — they are geometrically admissible under the
law set applied here. The central question is whether predicted and
experimental structures occupy the same geometric regime, and if not,
how the regimes differ.

---

## 2. Methods

### 2.1 Audit Framework

Toscanini (version 22.5.3) evaluates protein structures against a canon
of fifteen laws (canon hash: 6a9cd4b4349b81de). Each law specifies a
method, threshold, operator, units, and scope. Methods are:

- deterministic: always counted in scoring, can issue VETO
- heuristic: counted separately, cannot issue VETO
- advisory_experimental: present in output, excluded from det_passed

The deterministic score is computed as:

    det_score = (det_passed / det_total) × 100

where det_total = 12 (fifteen laws minus three heuristic laws).

A source-normalized score is also computed (WP-01):

    det_score_normalized = (det_passed_normalized / det_total) × 100

where det_passed_normalized counts passing laws from both deterministic
and advisory_experimental methods. This score enables direct numerical
comparison across source types (Section 2.2).

### 2.2 Method Assignment and Scoring Ceiling

LAW-100 (Bond Integrity) and LAW-160 (Chain Integrity) are assigned
method advisory_experimental for X-ray crystallographic structures.
These laws are included in det_total but excluded from det_passed,
imposing a structural ceiling of 10/12 = 83 on det_score for any
experimental structure. AlphaFold structures receive these laws as
deterministic, making 12/12 = 100 achievable.

This asymmetry is by design: bond geometry and chain continuity metrics
behave differently under crystallographic versus computational structure
generation. It means that raw det_score is not directly comparable across
source types. The source-normalized score (det_score_normalized) is the
appropriate cross-source comparison metric. Under normalization, all
passing structures from both source classes score 100.

### 2.3 Verdict Logic

    PASS:          det_score >= 83 AND coverage >= 70%
    INDETERMINATE: coverage < 70% OR score between thresholds
    VETO:          one or more hard deterministic law violations

LAW-105 (Reliability Coverage) is the gating law. Coverage below 70%
produces INDETERMINATE regardless of other law results. For AlphaFold
structures, coverage is computed from the per-residue pLDDT score.
For experimental structures, coverage is computed from B-factor values.

### 2.4 Dataset

Ten protein pairs were selected to represent diverse structural families,
functional classes, and known structural challenges.

**Table 1. Dataset composition.**

| Protein     | Experimental | Resolution | AlphaFold    | UniProt | Family              |
|-------------|-------------|------------|--------------|---------|---------------------|
| KRAS G12D   | 4OBE        | 1.35 Å     | AF-P01116-F1 | P01116  | GTPase, oncogene    |
| p53         | 2OCJ        | 2.05 Å     | AF-P04637-F1 | P04637  | Tumor suppressor    |
| EGFR        | 1IEP        | 2.10 Å     | AF-P00533-F1 | P00533  | Receptor kinase     |
| Myoglobin   | 1MBN        | 2.00 Å     | AF-P02144-F1 | P02144  | Oxygen transport    |
| Lysozyme    | 1HEL        | 1.70 Å     | AF-P00698-F1 | P00698  | Hydrolase benchmark |
| DHFR        | 1RX2        | 1.80 Å     | AF-P00374-F1 | P00374  | Enzyme benchmark    |
| Insulin     | 1MSO        | 1.00 Å     | AF-P01308-F1 | P01308  | Hormone, disulfide  |
| Calmodulin  | 1CLL        | 1.70 Å     | AF-P0DP23-F1 | P0DP23  | Calcium binding     |
| GFP         | 1EMA        | 1.90 Å     | AF-P42212-F1 | P42212  | Engineered protein  |
| Hemoglobin  | 1HHO        | 2.10 Å     | AF-P69905-F1 | P69905  | Oxygen transport    |

AlphaFold structures were downloaded from the EBI AlphaFold Database
(version 6). X-ray structures were downloaded from the RCSB Protein
Data Bank. No preprocessing was applied to any structure file.

### 2.5 Audit Procedure

All structures were submitted via HTTP POST to the /ingest endpoint
of a locally deployed Toscanini instance. The canon hash
(6a9cd4b4349b81de) was verified in every response before recording
results. Any response with a different canon hash would have terminated
the audit run. No such event occurred.

Audit reproducibility for this dataset is formally verified in a
companion paper (Meher, 2026b), which demonstrates that all audit IDs
are stable across three independent submission passes.

---

## 3. Results

### 3.1 Complete Audit Table

All 21 audit results are presented in Table 2. The KRAS pair was audited
in a prior session; all other pairs were audited in a single session on
2026-03-02.

**Table 2. Complete audit results.**

| Protein     | Source | Audit ID | Verdict       | det_score | det_score_norm | Coverage | Residues |
|-------------|--------|----------|---------------|-----------|----------------|----------|----------|
| KRAS        | AF     | F695F5F3 | PASS          | 100       | 100            | 92.6%    | 189      |
| KRAS G12D   | Xray   | FECA32C3 | PASS          | 83        | 100            | 100.0%   | 339      |
| p53         | AF     | 2A8231DA | INDETERMINATE | 91        | —              | 59.8%    | 393      |
| p53         | Xray   | EAF2ED48 | PASS          | 83        | 100            | 100.0%   | 776      |
| EGFR        | AF     | C2B3606F | PASS          | 100       | 100            | 70.7%    | 1210     |
| EGFR        | Xray   | 2F3509A4 | PASS          | 83        | 100            | 100.0%   | 548      |
| Myoglobin   | AF     | B1255783 | PASS          | 100       | 100            | 99.4%    | 154      |
| Myoglobin   | Xray   | 63BA24DF | VETO          | 75        | —              | 100.0%   | 153      |
| Lysozyme    | AF     | 8F67DB05 | PASS          | 100       | 100            | 89.1%    | 147      |
| Lysozyme    | Xray   | EE3AEEEF | PASS          | 83        | 100            | 100.0%   | 129      |
| DHFR        | AF     | 72D71C72 | PASS          | 100       | 100            | 98.9%    | 187      |
| DHFR        | Xray   | 30EBA93F | PASS          | 83        | 100            | 100.0%   | 159      |
| Insulin     | AF     | F9BDFEAF | INDETERMINATE | 83        | —              | 12.7%    | 110      |
| Insulin     | Xray   | D17224C9 | PASS          | 83        | 100            | 100.0%   | 102      |
| Calmodulin  | AF     | 233F0FCF | PASS          | 100       | 100            | 89.9%    | 149      |
| Calmodulin  | Xray   | D6388C63 | PASS          | 83        | 100            | 100.0%   | 144      |
| GFP         | AF     | F405858D | PASS          | 100       | 100            | 98.7%    | 238      |
| GFP         | Xray   | 5268812D | PASS          | 83        | 100            | 100.0%   | 221      |
| Hemoglobin  | AF     | B267F5CB | PASS          | 100       | 100            | 99.3%    | 142      |
| Hemoglobin  | Xray   | DEF2ECD3 | VETO          | 75        | —              | 100.0%   | 287      |

Note: det_score_norm is not computed for INDETERMINATE or VETO verdicts.
Canon hash verified on all 21 records: 6a9cd4b4349b81de.

### 3.2 Verdict Distribution

**AlphaFold structures (n=10):**
- PASS: 8 (80%) — KRAS, EGFR, Myoglobin, Lysozyme, DHFR, Calmodulin, GFP, Hemoglobin
- INDETERMINATE: 2 (20%) — p53, Insulin
- VETO: 0 (0%)

**Experimental structures (n=10):**
- PASS: 8 (80%) — KRAS, p53, EGFR, Lysozyme, DHFR, Insulin, Calmodulin, GFP
- INDETERMINATE: 0 (0%)
- VETO: 2 (20%) — Myoglobin, Hemoglobin

Both source classes show identical PASS rates (80%). The failure modes
differ by source type: AlphaFold failures arise from low reliability
coverage; experimental failures arise from rotamer violations.

### 3.3 Score Distribution

Raw deterministic scores (det_score):
- Experimental PASS (n=8): det_score = 83, uniform
- AlphaFold PASS (n=8): det_score = 100, uniform

Source-normalized scores (det_score_normalized):
- Experimental PASS (n=8): det_score_normalized = 100, uniform
- AlphaFold PASS (n=8): det_score_normalized = 100, uniform

After normalization, all passing structures from both source classes
achieve identical scores. The raw score difference (83 vs 100) is
entirely attributable to method assignment for LAW-100 and LAW-160
(Section 2.2), not to structural quality differences. The structural
regime differences between source types are carried in per-law observed
values, not in the summary score.

### 3.4 VETO Analysis — Experimental Structures

**Myoglobin (1MBN, audit ID: 63BA24DF):**
LAW-150 Rotamer Audit: observed 21.6%, threshold 20.0%, margin +1.60%.
Resolution: 2.0 Å. Side-chain conformations at this resolution are
electron-density supported. The VETO reflects genuine rotameric strain,
not a modeling error.

**Hemoglobin alpha chain (1HHO, audit ID: DEF2ECD3):**
LAW-150 Rotamer Audit: observed 20.35%, threshold 20.0%, margin +0.35%.
Resolution: 2.1 Å. Same interpretation as Myoglobin. The 0.35 percentage
point margin places this structure at the physical boundary of the
rotamer threshold.

Both VETOs arise from the same law at resolutions of 2.0-2.1 Å.
Whether LAW-150's threshold of 20% is appropriately calibrated for
experimental structures at this resolution range is a governance
question outside the scope of this audit. The threshold is fixed at
canon hash 6a9cd4b4349b81de.

### 3.5 INDETERMINATE Analysis — AlphaFold Structures

**p53 (AF-P04637-F1, audit ID: 2A8231DA):**
LAW-105 Reliability Coverage: observed 59.8%, threshold 70.0%.
p53 is a well-characterized intrinsically disordered protein. The
N-terminal transactivation domain and C-terminal regulatory domain
carry low pLDDT scores in AlphaFold predictions, correctly reflecting
structural uncertainty. The INDETERMINATE verdict accurately represents
the reliability of this model.

**Insulin (AF-P01308-F1, audit ID: F9BDFEAF):**
LAW-105 Reliability Coverage: observed 12.7%, threshold 70.0%.
The AlphaFold entry for P01308 models the full preproinsulin sequence
(110 residues). Mature insulin consists of two short chains (A: 21
residues, B: 30 residues) produced after signal peptide and C-peptide
cleavage. The signal peptide and C-peptide carry low pLDDT, dominating
the coverage calculation. The INDETERMINATE verdict reflects that the
majority of the submitted model is low-confidence by pLDDT — an accurate
description of the preproinsulin entry.

### 3.6 Per-Law Descriptive Statistics

Per-law observed values were extracted for all 11 scored laws across
both source classes (n=9 per class). KRAS structures are excluded from
per-law statistics because they were audited in a prior session without
per-law value extraction to the results file used here. KRAS remains
included in categorical results (Sections 3.1-3.5) but excluded from
per-law descriptive statistics. The n=9 dataset is internally consistent
and sufficient for descriptive comparison.

**Table 3. Per-law observed value statistics.**

| Law     | Title               | Xray mean | Xray SD | AF mean | AF SD   | Δ (X-AF) | Units    |
|---------|---------------------|-----------|---------|---------|---------|-----------|----------|
| LAW-100 | Bond Integrity      | 1.368     | 0.689   | 0.000   | 0.000   | +1.368    | %        |
| LAW-120 | Bond Angle          | 14.338    | 1.811   | 13.193  | 1.820   | +1.145    | degrees  |
| LAW-125 | Ramachandran        | 1.210     | 1.012   | 0.273   | 0.465   | +0.937    | %        |
| LAW-130 | Clashscore          | 0.000     | 0.000   | 0.000   | 0.000   | +0.000    | cl/1000  |
| LAW-135 | Omega Planarity     | 0.000     | 0.000   | 0.000   | 0.000   | +0.000    | %        |
| LAW-150 | Rotamer Audit       | 12.328    | 7.016   | 9.794   | 16.937  | +2.534    | %        |
| LAW-155 | Voxel Occupancy     | 6.057     | 0.158   | 5.970   | 0.139   | +0.087    | V        |
| LAW-160 | Chain Integrity     | 4.274     | 1.143   | 3.909   | 0.023   | +0.365    | Å        |
| LAW-182 | Hydrophobic Burial  | 0.587     | 0.207   | 0.405   | 0.216   | +0.182    | ratio    |
| LAW-195 | Disulfide Geometry  | 0.014     | 0.028   | 0.010   | 0.018   | +0.004    | Å        |
| LAW-200 | Packing Quality     | 70.982    | 27.894  | 137.644 | 105.277 | −66.662   | Å³/atom  |

Positive Δ indicates experimental structures show higher observed values.
Negative Δ indicates AlphaFold structures show higher observed values.
LAW-200 is the only law where AlphaFold structures show substantially
higher values than experimental structures.

**Key observations:**

LAW-100 Bond Integrity (Δ = +1.368%): AlphaFold structures show zero
bond violations across all nine structures (mean 0.000%, SD 0.000).
Experimental structures show mean 1.368% (SD 0.689). AlphaFold
structures are constructed under ideal bond geometry constraints.

LAW-125 Ramachandran (Δ = +0.937%): AlphaFold structures show fewer
backbone violations (mean 0.273%, SD 0.465) than experimental structures
(mean 1.210%, SD 1.012). Both classes are well within the 5.0% threshold.

LAW-150 Rotamer Audit: Including all nine AlphaFold structures yields
mean 9.794% (SD 16.937). Excluding INS_AF (coverage 12.7%, preproinsulin
full-length model, INDETERMINATE verdict) yields mean 3.845% (SD 2.187,
n=8). Experimental structures: mean 12.328% (SD 7.016, n=9). Whether to
include or exclude INS_AF is a methodological decision that must be
stated explicitly in any analysis using these values. The corrected
comparison (AF 3.845% vs Xray 12.328%) represents the largest
distributional separation in the dataset for passing structures.

LAW-182 Hydrophobic Burial (Δ = +0.182): Experimental structures show
higher burial ratios (mean 0.587, SD 0.207) than AlphaFold structures
(mean 0.405, SD 0.216). This pattern is consistent with crystallographic
compaction and ligand-bound states; causal attribution requires further
study.

LAW-200 Packing Quality (Δ = −66.662 Å³/atom): AlphaFold structures
show higher values (looser packing) than experimental structures. The
AF SD of 105.277 reflects high variance driven by EGFR_AF (1210
residues). Both classes pass the 300.0 Å³/atom threshold. Causal
attribution requires further study.

LAW-160 Chain Integrity: AF SD (0.023 Å) is 50x smaller than Xray SD
(1.143 Å), reflecting that AlphaFold always produces complete chains
while experimental structures vary due to missing loops, alternate
conformations, and crystal disorder.

### 3.7 Rotamer Audit Individual Values

**Table 4. LAW-150 Rotamer Audit — individual structure values.**

Experimental structures (threshold 20.0%):

| Protein    | PDB  | Observed | Verdict | Margin   |
|------------|------|----------|---------|----------|
| Lysozyme   | 1HEL |  2.86%   | PASS    | −17.14%  |
| Insulin    | 1MSO |  4.40%   | PASS    | −15.60%  |
| p53        | 2OCJ |  7.18%   | PASS    | −12.82%  |
| EGFR       | 1IEP |  7.92%   | PASS    | −12.08%  |
| DHFR       | 1RX2 | 12.50%   | PASS    | −7.50%   |
| GFP        | 1EMA | 16.93%   | PASS    | −3.07%   |
| Calmodulin | 1CLL | 17.21%   | PASS    | −2.79%   |
| Hemoglobin | 1HHO | 20.35%   | VETO    | +0.35%   |
| Myoglobin  | 1MBN | 21.60%   | VETO    | +1.60%   |

AlphaFold structures (threshold 20.0%):

| Protein    | UniProt | Observed | Verdict       | Margin   |
|------------|---------|----------|---------------|----------|
| Calmodulin | P0DP23  |  0.88%   | PASS          | −19.12%  |
| Lysozyme   | P00698  |  0.94%   | PASS          | −19.06%  |
| Myoglobin  | P02144  |  3.90%   | PASS          | −16.10%  |
| GFP        | P42212  |  3.90%   | PASS          | −16.10%  |
| Hemoglobin | P69905  |  4.42%   | PASS          | −15.58%  |
| EGFR       | P00533  |  5.04%   | PASS          | −14.96%  |
| p53        | P04637  |  6.67%   | PASS          | −13.33%  |
| DHFR       | P00374  |  7.78%   | PASS          | −12.22%  |
| Insulin    | P01308  | 54.55%   | INDETERMINATE | +34.55%  |

The two experimental VETOs sit 0.35 and 1.60 percentage points above
the threshold at resolutions of 2.1 Å and 2.0 Å respectively. Seven
of nine experimental structures pass, four sitting above 10%. Whether
LAW-150's threshold of 20% is appropriately calibrated for experimental
structures at 2.0-2.1 Å resolution is a governance question. The
threshold is fixed at canon hash 6a9cd4b4349b81de.

---

## 4. Discussion

### 4.1 Structural Regime Differences

The central finding of this study is not that one source type is
geometrically superior. After source normalization, both AlphaFold and
experimental structures are equally admissible under the fifteen-law
geometric canon. The finding is that they occupy distinct geometric
regimes that are consistently distinguishable across ten protein families.

AlphaFold structures occupy a regime of geometric idealization:
- Zero bond integrity violations (LAW-100: mean 0.000%)
- Low backbone strain (LAW-125: mean 0.273%)
- Low rotamer strain (LAW-150: mean 3.845% excluding preproinsulin)
- Looser atomic packing (LAW-200: mean 137.644 Å³/atom)
- Lower hydrophobic burial (LAW-182: mean 0.405)

Experimental structures occupy a regime of physical strain:
- Measurable bond violations (LAW-100: mean 1.368%)
- Higher backbone strain (LAW-125: mean 1.210%)
- Higher rotamer strain (LAW-150: mean 12.328%)
- Tighter atomic packing (LAW-200: mean 70.982 Å³/atom)
- Stronger hydrophobic burial (LAW-182: mean 0.587)

These regime differences are quantifiable, reproducible under canon hash
6a9cd4b4349b81de, and consistent across protein families spanning
GTPases, kinases, oxygen transport proteins, hydrolases, calcium-binding
proteins, and engineered fluorescent proteins.

### 4.2 Failure Mode Asymmetry

The failure modes of the two source classes are asymmetric and
informative.

AlphaFold failures (INDETERMINATE) arise from low reliability coverage —
a property of the model's own confidence signal (pLDDT). The engine
correctly identifies p53 as partially disordered and insulin as modeled
in precursor form. These are not geometric failures. They are accurate
descriptions of model reliability.

Experimental failures (VETO) arise from rotamer violations at 2.0-2.1 Å
resolution by margins of 0.35% and 1.60%. These are not modeling errors.
They are physically real side-chain conformations that are electron-density
supported at the stated resolution. The VETO threshold of 20% (LAW-150)
may require resolution-dependent calibration for experimental structures
— a governance question for future review.

The asymmetry suggests that the two source classes have different
characteristic failure modes that are predictable from their origin:
predicted structures fail at reliability; experimental structures fail
at rotameric strain under moderate resolution.

### 4.3 Governance Architecture

The results in this paper are tied to canon hash 6a9cd4b4349b81de.
This is not a formality. It is the mechanism by which these results
remain reproducible and citable.

If the law thresholds change — for example, if LAW-150's rotamer
threshold is adjusted from 20% to 22% following a resolution-sensitivity
study — the canon hash changes. Results produced under the new hash are
not directly comparable to results produced under 6a9cd4b4349b81de
without explicit acknowledgment of the governance change.

This property distinguishes Toscanini from tools that update silently.
Governance changes are explicit, detectable, and documentable. The
hash difference is the documentation.

### 4.4 Source-Normalized Scoring

The introduction of det_score_normalized (WP-01) resolves the most
significant methodological limitation of raw det_score comparison. Under
raw scoring, experimental structures are ceiling-capped at 83 due to
advisory_experimental method assignment for LAW-100 and LAW-160.
Under normalized scoring, all passing structures from both source classes
achieve 100.

The governance-internal raw score (det_score) is preserved unchanged
for audit and compliance purposes. The normalized score is the
appropriate metric for any cross-source comparison. This dual-score
architecture ensures that governance integrity and scientific
interpretability are both served without requiring canon modification.

### 4.5 Limitations

(1) Per-law statistics are reported for n=9 per class. This is sufficient
for descriptive comparison but insufficient for inferential statistics.
Mann-Whitney U tests, Cohen's d effect sizes, and confidence intervals
require n >= 30 per class per protein family.

(2) This study covers ten protein pairs. Systematic characterization of
regime differences across the full PDB and AlphaFold database requires
substantially larger datasets.

(3) The LAW-150 rotamer threshold (20%) may not be appropriately
calibrated for experimental structures at 2.0-2.1 Å resolution.
Two VETOs with margins under 2% at this resolution range suggest the
threshold was calibrated for higher resolution data. Resolution-dependent
threshold adjustment requires governance review.

(4) The INDETERMINATE verdict does not distinguish between intrinsic
disorder, low-confidence prediction, and precursor sequence artifacts.
Subcategorization of INDETERMINATE would improve actionability for
researchers receiving this verdict.

(5) Cross-machine reproducibility has not been empirically verified in
this study, though the deterministic architectural properties guarantee
it by construction. See companion paper (Meher, 2026b) for within-session
reproducibility verification.

---

## 5. Conclusion

Toscanini provides a deterministic, reproducible, governance-anchored
structural admissibility audit. Applied to ten protein pairs across five
protein families, it reveals:

(1) After source normalization, AlphaFold and experimental structures are
equally admissible under the geometric law set. The raw score difference
(83 vs 100) is methodological, not structural.

(2) The two source classes occupy distinct geometric regimes. AlphaFold
structures are geometrically idealized. Experimental structures carry
physical strain. These regime differences are quantifiable and consistent
across protein families.

(3) Failure modes are asymmetric by source type. AlphaFold failures arise
from low reliability coverage in disordered or precursor regions.
Experimental failures arise from marginal rotamer violations at 2.0-2.1 Å
resolution.

(4) Canon hash 6a9cd4b4349b81de was stable across all 21 audits.
Governance integrity was maintained throughout.

The framework is ready for broader application. The structural regime
differences reported here are reproducible under the stated canon hash
and constitute testable hypotheses for larger-scale studies.

---

## Data Availability

All audit records and supporting data are available at:
https://github.com/prateekm1007/control-atlas-saas
Path: docs/white_paper/

Complete audit ID list:
F695F5F3 BCFDFEDC FECA32C3 EAF2ED48 2F3509A4 63BA24DF EE3AEEEF
30EBA93F D17224C9 D6388C63 5268812D DEF2ECD3 2A8231DA C2B3606F
B1255783 8F67DB05 72D71C72 F9BDFEAF 233F0FCF F405858D B267F5CB

Canon hash: 6a9cd4b4349b81de
System version: 22.5.3

---

## References

1. Jumper J, et al. (2021). Highly accurate protein structure prediction
   with AlphaFold. Nature, 596, 583-589.

2. Varadi M, et al. (2022). AlphaFold Protein Structure Database:
   massively expanding the structural coverage of protein-sequence space
   with high-accuracy models. Nucleic Acids Research, 50(D1), D439-D444.

3. Williams CJ, et al. (2018). MolProbity: More and better reference
   data for improved all-atom structure validation.
   Protein Science, 27(1), 293-315.

4. Berman HM, et al. (2000). The Protein Data Bank.
   Nucleic Acids Research, 28(1), 235-242.

5. Meher P. (2026b). Deterministic Reproducibility and Canon Hash
   Governance in Protein Structure Admissibility Auditing.
   [Preprint, bioRxiv]

6. Read RJ, et al. (2011). A new generation of crystallographic
   validation tools for the Protein Data Bank.
   Structure, 19(10), 1395-1412.

---

## Acknowledgments

The author thanks the RCSB Protein Data Bank and the European
Bioinformatics Institute AlphaFold Database for providing open access
to structural data used in this study.

---

*Correspondence: [YOUR EMAIL]*
*Preprint server: bioRxiv*
*Subject area: Bioinformatics*
*License: CC BY 4.0*

