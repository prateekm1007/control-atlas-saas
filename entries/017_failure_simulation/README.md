# Entry #017 — Failure Simulation (Resistance Prediction)

## Prerequisites
- Entry 015 LOCKED (state manifold validated)
- Entry 016 LOCKED (physics gate validated)

## Objective
Predict control failure modes BEFORE chemistry by simulating pocket perturbations.

## Core Question
> "Which mutations collapse the druggable pocket?"

## Method
1. Load drug-bound structure (6OIM - highest exposure/volume)
2. Identify pocket-lining residues (Switch II contact shell)
3. Simulate perturbations:
   - Point mutations (bulk change: small→large, polar→hydrophobic)
   - Side chain deletion (Ala scanning proxy)
4. Recalculate exposure/volume after each perturbation
5. Flag mutations that collapse pocket metrics below threshold

## Perturbation Types
| Type | Simulation | Expected Effect |
|------|------------|-----------------|
| Steric block | Gly→Trp | Volume collapse |
| Charge flip | Asp→Lys | Exposure change |
| Anchor loss | Tyr→Ala | Binding site disruption |

## Known Resistance Mutations (Validation Set)
| Mutation | Clinical | Expected Metric Change |
|----------|----------|------------------------|
| Y96D | Sotorasib resistance | Exposure drop |
| R68S | Reported | Volume change |
| Q99L | Emerging | Hydrophobic shift |

## Success Criteria
- [ ] Pocket-lining residues identified
- [ ] Perturbation simulation implemented
- [ ] Known resistance mutations show metric collapse
- [ ] Novel resistance candidates predicted

## Abort Condition
If known resistance mutations show NO metric change → model insufficient

## Output
- 
