# Entry #016 — Pocket Physics Discrimination

## Objective
Prove that concavity + polarity distinguish true control sites from geometric decoys.

## Core Claim
> Geometry finds candidates; physics decides control.

## Method
1. Load matched structures (H-Ras, RhoA, Ubiquitin)
2. Render surfaces on Switch II equivalent regions only
3. Apply Molecular Lipophilicity Potential (MLP)
4. Visual comparison: concave/polar pocket vs convex/uniform surface

## Success Criteria
| Protein | Expected Surface | Expected Polarity | Result |
|---------|------------------|-------------------|--------|
| H-Ras (5P21) | Concave pocket | Mixed hydrophobic-polar | PASS |
| Ubiquitin (1UBQ) | Convex/flat | Uniform | REJECT |

## Lock Condition
Visual contrast obvious without explanation.
