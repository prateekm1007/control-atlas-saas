# Asset Dossier: CHAMP-005

## Identity
- **Sequence:** KAWAKKAAAKAEAAKAEAAKYWPTG
- **Length:** 25 aa
- **Architecture:** Helical peptide (AEAAK repeats + YWPTG warhead)

## Tier 1 Geometric Audit

### Pre-Minimization (Raw Chai-1)
| Metric | Value |
|--------|-------|
| Min Distance | 1.63 Å |
| Clashes (< 2.5Å) | 7 |
| Clashes (< 1.9Å) | 2 |
| Contact Density | 199 |
| Status | CLASH_VETO |

### Post-Minimization (100 steps Amber14/OBC2)
| Metric | Value |
|--------|-------|
| Min Distance | **2.86 Å** |
| Clashes (< 2.5Å) | **0** |
| Clashes (< 1.9Å) | **0** |
| Contact Density | 147 |
| Energy | -21,181 kJ/mol |
| Status | **✅ SOVEREIGN_PASS** |

## Important Note
Severe pre-minimization clashes resolved by peptide repositioning,
resulting in 26% loss of interface contacts (199 → 147).

This suggests the binding mode may differ from the raw Chai-1 prediction.

## Verdict
**TIER 1: PASSED** after energy minimization, with caveat of reduced contacts.

## Artifact Provenance
- **Raw Structure:** structures/champ005/structure.cif
- **Minimized:** structures/champ005/minimized.pdb
- **Rescue Results:** structures/champ005/rescue_results.json
- **SHA256 (raw):** 8c0e283004eda5a18f711fc0c3e5e59ea26d49571af1a8516e4136386002e131
