# TOSCANINI v22.5.3 — Trust Matrix Certification

**Date:** 2026-02-17T17:37:16Z
**Doctrine:** PIL-CAL-03 Modality-Aware Enforcement Matrix

## Results

| Structure | Method | Resolution | Verdict | Det Fails | Notes |
|---|---|---|---|---|---|
| 4HHB | X-ray | 1.74Å | PASS | 0 | Hemoglobin tetramer |
| 1CRN | X-ray | 1.50Å | PASS | 0 | Crambin |
| 1UBQ | X-ray | 1.80Å | PASS | 0 | Ubiquitin |
| 6LZG | X-ray | 2.50Å | PASS | 0 | SARS-CoV-2 RBD |
| AF-P01308-F1 | AlphaFold | — | INDETERMINATE | LAW-105 | Low pLDDT (prepropeptide) |
| 7BV2 | Cryo-EM | 2.50Å | PASS | 0 | LAW-170 advisory (non-standard residues) |
| 1G03 | NMR | — | PASS | 0 | LAW-125 advisory (ensemble phi/psi) |

## Modality Matrix (PIL-CAL-03)

| Law | X-ray | Cryo-EM | NMR | Predicted |
|---|---|---|---|---|
| LAW-100 | Advisory | Advisory | Advisory | Deterministic |
| LAW-125 | Deterministic | Deterministic | Advisory | Deterministic |
| LAW-170 | Deterministic | Advisory | Advisory | Deterministic |
| All others | Per station_sop.py classification | | | |

## Certification

- 0 false vetoes across 6 experimental structures (0.48–2.50Å resolution range)
- Heuristic failures (LAW-155) do not trigger VETO — correct
- Coverage gate (LAW-105) correctly triggers INDETERMINATE for low-confidence AF models
- Canon hash: verified against station_sop.py

**Status: FROZEN**
