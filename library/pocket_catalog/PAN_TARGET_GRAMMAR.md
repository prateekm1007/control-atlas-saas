# Pan-Target Grammar v1.0

## Overview
Universal druggability rules derived from 15 validated oncogenic pockets across 7 pocket classes.

## Pocket Classes

### SwitchII (RAS family)
- Targets: KRAS G12C/G12D, HRAS G12V, NRAS Q61R
- Large volume, high exposure
- Covalent warhead compatible

### KinaseHinge (RTK family)
- Targets: BRAF V600E, EGFR L858R, ALK F1174L, RET M918T
- Medium-large volume, moderate exposure
- ATP-competitive scaffold compatible

### ActivationLoop
- Targets: KIT D816V, PDGFRA D842V
- Medium volume, balanced properties

### DNABinding (Tumor suppressors)
- Targets: TP53 R175H/R248Q
- Large volume, lower exposure

### Other Classes
- PHDomain: AKT1 E17K
- Kinase: PIK3CA H1047R
- PhosphoSite: CTNNB1 S45F

## Universal Chemistry Rules

| Property | Range | Rationale |
|----------|-------|-----------|
| MW | 250-650 Da | Pocket volume correlation |
| Rings | 2-5 | Core rigidity requirement |
| HBA | 2-10 | Exposure-dependent |
| HBD | 0-4 | Permeability balance |
| cLogP | 1-6 | Hydrophobicity match |
| RotBonds | ≤10 | Entropy penalty |

## Class-Specific Modifiers

| Class | Core Size | Polar Tolerance | Anchor |
|-------|-----------|-----------------|--------|
| SwitchII | Large | High | Weak-Medium |
| KinaseHinge | Medium-Large | Medium | Medium |
| ActivationLoop | Small-Medium | Medium | Strong |
| DNABinding | Large | Low-Medium | Weak |

## Application

1. **New Target Assessment**: Compute physics → classify → apply grammar
2. **Scaffold Design**: Use class rules to constrain generation
3. **Library Filtering**: Apply universal + class rules to prioritize compounds
