# Deterministic Reproducibility and Canon Hash Governance
# in Protein Structure Admissibility Auditing

## Abstract

Reproducibility is a foundational requirement for computational structural
biology tools. We demonstrate that Toscanini, a deterministic structural
admissibility authority, produces identical audit identifiers, scores, and
verdicts across repeated submissions of the same protein structure files.
Ten protein structures spanning five pairs of AlphaFold predictions and
X-ray crystallographic structures were submitted three times each to a
locally deployed Toscanini instance. All thirty audit records are identical
across passes. The governance canon hash (6a9cd4b4349b81de) remained stable
throughout. We describe the architectural mechanisms that guarantee this
property and propose audit ID verification as a standard reproducibility
check for deterministic structural validation tools.

---

## 1. Introduction

Computational tools for protein structure validation are increasingly central
to structural biology workflows. AlphaFold has produced predictions for
virtually the entire known proteome. Cryo-EM structures are deposited at
unprecedented rates. Validation tools must keep pace — and must themselves
be reproducible.

Reproducibility in structural validation has two distinct requirements:

1. Numerical reproducibility: the same structure produces the same scores.
2. Governance reproducibility: the law set used to produce those scores is
   fixed, versioned, and auditable.

Most existing tools satisfy the first requirement implicitly. None, to our
knowledge, formally address the second. Toscanini addresses both through
a deterministic audit engine and a cryptographically fingerprinted law canon.

This paper provides formal verification of both properties across ten protein
structures and three independent submission passes.

---

## 2. Methods

### 2.1 System

Toscanini version 22.5.3, deployed as a containerized FastAPI service.
Endpoint: POST /ingest. Authentication: API key header.

### 2.2 Audit ID Construction

The audit ID is an 8-character uppercase hexadecimal identifier derived
deterministically from the input PDB file content and the law canon.
The same file submitted to the same canon always produces the same audit ID.
Audit IDs serve as cryptographic anchors — a change in audit ID indicates
either a change in input or a change in canon.

### 2.3 Canon Hash

The law canon is a fixed set of fifteen physics-grounded laws with defined
thresholds, operators, units, and scope. The canon is hashed at system
initialization. The hash is embedded in every audit record under:

    governance.governance_fingerprint.canon_hash

Canon hash for this study: 6a9cd4b4349b81de

Any modification to any law definition, threshold, or metadata changes
this hash. Results produced under different canon hashes are not comparable.

### 2.4 Dataset

Ten protein structures were selected from the white paper dataset
(Toscanini White Paper, 2026-03-02):

| Structure  | Source      | PDB/AF ID    | Mode         |
|------------|-------------|--------------|--------------|
| KRAS G12D  | X-ray       | 4OBE         | experimental |
| KRAS       | AlphaFold   | AF-P01116-F1 | predicted    |
| p53        | X-ray       | 2OCJ         | experimental |
| p53        | AlphaFold   | AF-P04637-F1 | predicted    |
| EGFR       | X-ray       | 1IEP         | experimental |
| EGFR       | AlphaFold   | AF-P00533-F1 | predicted    |
| Myoglobin  | X-ray       | 1MBN         | experimental |
| Myoglobin  | AlphaFold   | AF-P02144-F1 | predicted    |
| Lysozyme   | X-ray       | 1HEL         | experimental |
| Lysozyme   | AlphaFold   | AF-P00698-F1 | predicted    |

### 2.5 Procedure

Each structure was submitted three times via POST /ingest with identical
parameters. Submissions were sequential within a single session on a single
machine. The audit ID, deterministic score, and canon hash were recorded
for each submission. Pass files were compared using Unix diff.

---

## 3. Results

### 3.1 Audit ID Stability

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

10/10 structures stable across all three passes.
Unix diff between all pass files: empty (no differences).

### 3.2 Canon Hash Stability

Canon hash 6a9cd4b4349b81de was observed in all 30 audit records.
No governance drift detected across any submission.

### 3.3 Score Stability

Deterministic scores were identical across all three passes for all
ten structures:

| Structure  | Score | Verdict       |
|------------|-------|---------------|
| KRAS_4OBE  | 83    | PASS          |
| KRAS_iso1  | 100   | PASS          |
| p53_2OCJ   | 83    | PASS          |
| p53_AF     | 91    | INDETERMINATE |
| EGFR_1IEP  | 83    | PASS          |
| EGFR_AF    | 100   | PASS          |
| MYO_1MBN   | 75    | VETO          |
| MYO_AF     | 100   | PASS          |
| LYS_1HEL   | 83    | PASS          |
| LYS_AF     | 100   | PASS          |

No score variance observed. No verdict variance observed.

---

## 4. Discussion

### 4.1 Mechanisms of Determinism

Toscanini achieves deterministic reproducibility through three architectural
properties:

**Input hashing:** The audit ID is derived from a hash of the input file
content. Identical files produce identical hashes. Any byte-level difference
in input produces a different audit ID — making input identity verifiable.

**Canon freezing:** The law set is loaded once at system initialization and
hashed. The canon hash is embedded in every audit record. Operators can
verify that two audit records were produced under the same law set by
comparing canon hashes. Records produced under different hashes are not
comparable and should not be aggregated.

**Stateless evaluation:** Each audit is independent. No state from prior
audits influences subsequent ones. The audit pipeline is a pure function
from input to output.

### 4.2 Implications for Structural Biology

The audit ID functions as a structural fingerprint. Two researchers
submitting the same PDB file to Toscanini will receive the same audit ID,
the same score, and the same verdict — regardless of when or where the
submission occurs, provided the canon hash is identical.

This property enables:

1. **Citation by audit ID:** A paper can cite a specific structural audit
   by its ID. Any reader can reproduce the result by submitting the same
   file to the same canon version.

2. **Governance traceability:** If the law set changes, the canon hash
   changes. Results from different canon versions are explicitly
   incomparable — not silently inconsistent.

3. **Collaborative auditing:** Multiple groups can audit the same structure
   independently and verify agreement by comparing audit IDs, without
   sharing raw computation.

### 4.3 Limitations

This study demonstrates reproducibility within a single session on a single
machine. Cross-machine and cross-session reproducibility has not been
formally verified here. The deterministic properties described suggest
cross-machine reproducibility is guaranteed by construction, but empirical
verification across independent deployments is needed.

The study covers ten structures. Larger-scale verification across hundreds
of structures is the appropriate next step.

---

## 5. Conclusion

Toscanini produces identical audit IDs, scores, and verdicts across repeated
submissions of the same protein structure files. Canon hash
6a9cd4b4349b81de remained stable across all 30 audit records in this study.
Unix diff between independent pass files is empty.

We propose that audit ID stability and canon hash verification be adopted
as standard reproducibility checks for deterministic structural validation
tools. These checks are binary, fast, and unambiguous — properties that
narrative reproducibility statements cannot match.

---

## Appendix: Raw Reproducibility Data

Pass 1:
KRAS_4OBE  FECA32C3  83  6a9cd4b4349b81de
KRAS_iso1  F695F5F3  100 6a9cd4b4349b81de
p53_2OCJ   EAF2ED48  83  6a9cd4b4349b81de
p53_AF     2A8231DA  91  6a9cd4b4349b81de
EGFR_1IEP  2F3509A4  83  6a9cd4b4349b81de
EGFR_AF    C2B3606F  100 6a9cd4b4349b81de
MYO_1MBN   63BA24DF  75  6a9cd4b4349b81de
MYO_AF     B1255783  100 6a9cd4b4349b81de
LYS_1HEL   EE3AEEEF  83  6a9cd4b4349b81de
LYS_AF     8F67DB05  100 6a9cd4b4349b81de

Pass 2: identical to Pass 1.
Pass 3: identical to Pass 1.

diff pass_1.jsonl pass_2.jsonl: (empty)
diff pass_2.jsonl pass_3.jsonl: (empty)

System: Toscanini v22.5.3
Canon:  6a9cd4b4349b81de
Date:   2026-03-02
