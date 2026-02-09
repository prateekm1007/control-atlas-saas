# Entry #018 — Chemistry Grammar

## Prerequisites
- Entry 016 LOCKED (pocket physics)
- Entry 017 LOCKED (resistance prediction)

## Objective
Define rules mapping pocket geometry to compatible chemical scaffolds.

## Core Question
> What chemical features can address this control pocket?

## Method
1. Analyze pocket properties (volume, polarity, exposure)
2. Define pharmacophore requirements
3. Map to chemical scaffold classes
4. Output actionable chemistry rules

## Pocket Properties (from 015/016)
| Property | Drug-Open State |
|----------|-----------------|
| Volume | 2231 A³ |
| Exposure | 36.63 |
| Hydrophobic | 18.8% |

## Pharmacophore Rules
| Pocket Feature | Required Chemistry |
|----------------|-------------------|
| Concave cavity | Rigid aromatic core |
| Partial burial | Hydrophobic anchor |
| Polar edge | H-bond donor/acceptor |
| Cys12 (G12C) | Electrophile warhead |

## Output
- chemistry_rules.json
- scaffold_compatibility.csv

## Status
[ ] Pocket features extracted
[ ] Pharmacophore rules defined
[ ] Scaffold classes mapped
[ ] Rules validated against known drugs
