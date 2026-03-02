---
title: >
  Deterministic Reproducibility and Canon Hash Governance
  in Protein Structure Admissibility Auditing
authors:
  - name: Prateek Meher
    affiliation: Independent Researcher
    email: [YOUR EMAIL]
date: 2026-03-02
version: 1.0
canon_hash: 6a9cd4b4349b81de
system_version: 22.5.3
---

# Deterministic Reproducibility and Canon Hash Governance in Protein Structure Admissibility Auditing

**Prateek Meher**
Independent Researcher
[YOUR EMAIL]

Preprint. Not peer reviewed.
bioRxiv submission: 2026-03-02

---

## Abstract

Reproducibility is a foundational requirement for computational structural
biology tools, yet most validation frameworks provide no formal mechanism
for verifying that results are stable across runs, machines, or time.
We demonstrate that Toscanini, a deterministic structural admissibility
authority, produces identical audit identifiers, scores, and verdicts
across repeated independent submissions of the same protein structure files.
Ten protein structures — spanning five pairs of AlphaFold v6 predictions
and high-resolution X-ray crystallographic structures — were submitted three
times each to a locally deployed Toscanini instance (version 22.5.3).
All thirty audit records are identical across passes. The governance canon
hash (6a9cd4b4349b81de) remained stable throughout all submissions. We
describe the three architectural mechanisms that guarantee this property:
input content hashing, canon freezing, and stateless evaluation. We further
propose audit ID verification and canon hash comparison as standard
reproducibility checks for deterministic structural validation tools. These
checks are binary, fast, and unambiguous — properties that narrative
reproducibility statements cannot provide.

**Keywords:** protein structure validation, reproducibility, deterministic
auditing, structural bioinformatics, AlphaFold, canonical governance

---

## 1. Introduction

Protein structure validation tools are central to structural biology
workflows. MolProbity, PROCHECK, and the wwPDB validation pipeline are
widely used for assessing geometric quality of experimentally determined
and computationally predicted structures. These tools produce scores and
reports that are cited in publications and required for PDB deposition.

However, the reproducibility of these tools is rarely formally verified.
Reproducibility in structural validation has two distinct requirements that
are often conflated:

(1) Numerical reproducibility: the same structure produces the same scores.
(2) Governance reproducibility: the law set used to produce those scores is
fixed, versioned, and auditable across time and across users.

Most tools satisfy the first requirement implicitly through deterministic
algorithms. The second requirement — governance reproducibility — is almost
never formally addressed. When a validation tool updates its scoring
parameters or thresholds, results from before and after the update are
silently inconsistent. There is no mechanism for a reader to verify that
a validation result cited in a 2023 paper used the same law set as one
cited in a 2026 paper.

Toscanini addresses both requirements through a deterministic audit engine
and a cryptographically fingerprinted law canon. Every audit record embeds
the canon hash — a fixed identifier for the exact law set used. If the law
set changes, the hash changes. Results from different hashes are explicitly
incomparable, not silently inconsistent.

This paper provides formal empirical verification of both reproducibility
properties across ten protein structures and three independent submission
passes, and proposes audit ID stability as a standard reproducibility check
for deterministic structural validation tools.

---

## 2. Methods

### 2.1 System Description

Toscanini (version 22.5.3) is a deterministic structural admissibility
authority deployed as a containerized FastAPI service. It evaluates protein
structures against a fixed canon of fifteen physics-grounded laws covering
bond geometry, backbone continuity, stereochemistry, packing quality, and
reliability coverage. Each law specifies a method (deterministic, heuristic,
or advisory_experimental), a threshold, units, and scope.

The system returns a structured JSON response containing:
- A binary verdict (PASS, VETO, or INDETERMINATE)
- A deterministic score (0-100)
- An 8-character audit ID
- A governance fingerprint containing the canon hash
- Per-law observed values and pass/fail status

### 2.2 Audit ID Construction

The audit ID is an 8-character uppercase hexadecimal identifier derived
deterministically from the byte content of the input PDB file in
combination with the law canon. Identical files submitted to the same canon
always produce the same audit ID. Any byte-level difference in the input
file produces a different audit ID. This property makes input identity
verifiable without transmitting the file itself.

### 2.3 Canon Hash

The law canon is loaded once at system initialization and hashed. The
resulting hash is embedded in every audit record under the field:

    governance.governance_fingerprint.canon_hash

The canon hash for all experiments in this study is: 6a9cd4b4349b81de

Any modification to any law definition, threshold, operator, or metadata
changes this hash. The canon hash therefore serves as a version identifier
for the governance framework — enabling readers to verify that two audit
results were produced under identical evaluation conditions.

### 2.4 Dataset

Ten protein structures were selected from an ongoing comparative study of
AlphaFold predictions and X-ray crystallographic structures (Meher, 2026a).
The dataset was designed to include diverse protein families, both source
types, and structures spanning the full range of verdicts (PASS, VETO,
INDETERMINATE).

| Label      | Source      | PDB/UniProt  | Resolution | Verdict       |
|------------|-------------|--------------|------------|---------------|
| KRAS_4OBE  | X-ray       | 4OBE         | 1.35 Å     | PASS          |
| KRAS_iso1  | AlphaFold   | P01116       | N/A        | PASS          |
| p53_2OCJ   | X-ray       | 2OCJ         | 2.05 Å     | PASS          |
| p53_AF     | AlphaFold   | P04637       | N/A        | INDETERMINATE |
| EGFR_1IEP  | X-ray       | 1IEP         | 2.10 Å     | PASS          |
| EGFR_AF    | AlphaFold   | P00533       | N/A        | PASS          |
| MYO_1MBN   | X-ray       | 1MBN         | 2.00 Å     | VETO          |
| MYO_AF     | AlphaFold   | P02144       | N/A        | PASS          |
| LYS_1HEL   | X-ray       | 1HEL         | 1.70 Å     | PASS          |
| LYS_AF     | AlphaFold   | P00698       | N/A        | PASS          |

AlphaFold structures were downloaded from the European Bioinformatics
Institute AlphaFold Protein Structure Database (version 6). X-ray
structures were downloaded from the RCSB Protein Data Bank.

### 2.5 Submission Procedure

Each structure was submitted three times via HTTP POST to the /ingest
endpoint with identical parameters (mode, candidate_id, file). Submissions
were sequential within a single session on a single machine running
Ubuntu 22.04 via Windows Subsystem for Linux 2. The audit ID, deterministic
score, and canon hash were extracted from each response and recorded to a
plain-text log file. The three pass files were compared using Unix diff.

No preprocessing was applied to any structure file. Files were submitted
as downloaded from their respective sources.

---

## 3. Results

### 3.1 Audit ID Stability

All ten structures produced identical audit IDs across all three
independent submission passes. The results are presented in Table 1.

**Table 1. Audit ID stability across three submission passes.**

| Structure  | Pass 1   | Pass 2   | Pass 3   | Stable |
|------------|----------|----------|----------|--------|
| KRAS_4OBE  | FECA32C3 | FECA32C3 | FECA32C3 | YES    |
| KRAS_iso1  | F695F5F3 | F695F5F3 | F695F5F3 | YES    |
| p53_2OCJ   | EAF2ED48 | EAF2ED48 | EAF2ED48 | YES    |
| p53_AF     | 2A8231DA | 2A8231DA | 2A8231DA | YES    |
| EGFR_1IEP  | 2F3509A4 | 2F3509A4 | 2F3509A4 | YES    |
| EGFR_AF    | C2B3606F | C2B3606F | C2B3606F | YES    |
| MYO_1MBN   | 63BA24DF | 63BA24DF | 63BA24DF | YES    |
| MYO_AF     | B1255783 | B1255783 | B1255783 | YES    |
| LYS_1HEL   | EE3AEEEF | EE3AEEEF | EE3AEEEF | YES    |
| LYS_AF     | 8F67DB05 | 8F67DB05 | 8F67DB05 | YES    |

Stability rate: 10/10 (100%).
Unix diff between pass_1 and pass_2 log files: empty.
Unix diff between pass_2 and pass_3 log files: empty.

### 3.2 Canon Hash Stability

The canon hash 6a9cd4b4349b81de was present and identical in all 30
audit records across all three passes. No governance drift was detected
at any point during the experiment.

### 3.3 Score and Verdict Stability

Deterministic scores and verdicts were identical across all three passes
for all ten structures (Table 2). No score variance or verdict variance
was observed.

**Table 2. Score and verdict stability across three submission passes.**

| Structure  | Score | Verdict       | Consistent |
|------------|-------|---------------|------------|
| KRAS_4OBE  | 83    | PASS          | YES        |
| KRAS_iso1  | 100   | PASS          | YES        |
| p53_2OCJ   | 83    | PASS          | YES        |
| p53_AF     | 91    | INDETERMINATE | YES        |
| EGFR_1IEP  | 83    | PASS          | YES        |
| EGFR_AF    | 100   | PASS          | YES        |
| MYO_1MBN   | 75    | VETO          | YES        |
| MYO_AF     | 100   | PASS          | YES        |
| LYS_1HEL   | 83    | PASS          | YES        |
| LYS_AF     | 100   | PASS          | YES        |

The dataset intentionally includes all three verdict types (PASS, VETO,
INDETERMINATE) to demonstrate that reproducibility holds regardless of
outcome. The VETO verdict for MYO_1MBN (Myoglobin, 1MBN) and the
INDETERMINATE verdict for p53_AF (p53, AF-P04637-F1) are reproduced
exactly across all three passes.

---

## 4. Discussion

### 4.1 Mechanisms of Determinism

Toscanini achieves deterministic reproducibility through three
architectural properties that operate independently and in combination.

**Input content hashing.** The audit ID is derived from a hash of the
input file byte content. Identical files always produce identical hashes.
Any byte-level difference in the input — including whitespace, line
endings, or metadata differences — produces a different audit ID. This
makes input identity verifiable without file transmission: two researchers
can confirm they audited the same structure by comparing audit IDs.

**Canon freezing.** The law set is loaded once at system initialization
and hashed. The canon hash is embedded in every audit record. Researchers
can verify that two audit records were produced under identical evaluation
conditions by comparing canon hashes. Records produced under different
hashes used different law sets and should not be aggregated or compared
without explicit acknowledgment.

**Stateless evaluation.** Each audit is independent. No state from prior
audits influences subsequent ones. The audit pipeline is a pure function
from (input file, canon) to (audit result). There are no session-dependent
variables, no cache effects, and no order dependencies between submissions.

### 4.2 Audit ID as a Reproducibility Primitive

We propose that the audit ID functions as a structural reproducibility
primitive — a compact, verifiable identifier that encodes both the input
and the evaluation context.

Current practice in structural biology papers is to describe validation
results narratively: "MolProbity score of X" or "Ramachandran outliers
Y%." These descriptions cannot be independently verified without access
to the same tool version, the same parameter set, and the same input file.

An audit ID citation — "structure validated under audit ID FECA32C3,
canon hash 6a9cd4b4349b81de" — is verifiable. Any reader with access to
Toscanini and the original PDB file can reproduce the exact result in
under five minutes.

This property is particularly valuable for:

(1) Journal reviewers verifying validation claims without re-running
    complete analyses.

(2) Authors defending validation results across revision cycles where
    the structure may be refined between submissions.

(3) Database curators maintaining long-term validation records where
    tool version consistency is critical.

### 4.3 Canon Hash as a Governance Primitive

The canon hash addresses a problem that is currently invisible in
structural biology: silent inconsistency between validation runs
produced under different tool versions.

When MolProbity updates its Ramachandran reference data, results from
before and after the update are inconsistent — but there is no standard
mechanism for detecting or flagging this inconsistency in published
results. Papers from different years citing MolProbity scores are not
guaranteed to be comparable.

Toscanini's canon hash makes this explicit. Results produced under
different hashes used different law sets. The hash difference is
detectable, documentable, and can be cited. Governance drift is not
silent — it changes a 16-character string that appears in every audit
record.

### 4.4 Limitations

This study demonstrates reproducibility within a single session on a
single machine. Cross-machine reproducibility — submitting the same
files to independent Toscanini deployments on different hardware — has
not been empirically verified in this study. The deterministic
architectural properties described in Section 4.1 are designed to
guarantee cross-machine reproducibility, but empirical verification
across independent deployments is needed.

The study covers ten structures across three passes. Larger-scale
verification — hundreds of structures across multiple machines and
multiple time points — would strengthen the empirical foundation.

The audit ID is derived from file byte content, not from parsed
coordinate data. Two PDB files containing identical atomic coordinates
but different header metadata will produce different audit IDs. Users
should ensure file consistency before comparing audit IDs across
sources.

---

## 5. Conclusion

Toscanini produces identical audit IDs, scores, and verdicts across
repeated independent submissions of the same protein structure files.
In this study, ten structures were submitted three times each, producing
thirty audit records. All audit IDs are identical across passes. Canon
hash 6a9cd4b4349b81de remained stable throughout. Unix diff between
independent pass log files is empty.

We propose audit ID stability and canon hash verification as standard
reproducibility checks for deterministic structural validation tools.
These checks are binary — either the audit IDs match or they do not —
and require no statistical interpretation.

The architectural mechanisms described here — input hashing, canon
freezing, and stateless evaluation — are transferable design principles
for any deterministic computational biology tool where reproducibility
and governance traceability are required.

---

## Data Availability

All pass log files are available at:
https://github.com/prateekm1007/control-atlas-saas
Path: docs/white_paper/reproducibility/

Raw audit records:
Pass 1 audit IDs: FECA32C3 F695F5F3 EAF2ED48 2A8231DA 2F3509A4
                  C2B3606F 63BA24DF B1255783 EE3AEEEF 8F67DB05
Pass 2: identical to Pass 1.
Pass 3: identical to Pass 1.

Canon hash: 6a9cd4b4349b81de
System version: 22.5.3

---

## References

1. Meher P. (2026a). Deterministic Structural Admissibility Auditing:
   A Reproducible Framework for Evaluating AlphaFold Predictions Against
   Experimental Ground Truth. [Preprint, bioRxiv]

2. Williams CJ, et al. (2018). MolProbity: More and better reference
   data for improved all-atom structure validation.
   Protein Science, 27(1), 293-315.

3. Jumper J, et al. (2021). Highly accurate protein structure prediction
   with AlphaFold. Nature, 596, 583-589.

4. Berman HM, et al. (2000). The Protein Data Bank.
   Nucleic Acids Research, 28(1), 235-242.

5. Varadi M, et al. (2022). AlphaFold Protein Structure Database:
   massively expanding the structural coverage of protein-sequence space
   with high-accuracy models. Nucleic Acids Research, 50(D1), D439-D444.

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

