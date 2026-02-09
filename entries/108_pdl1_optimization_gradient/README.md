# Entry 108 — PD-L1 Optimization Gradient

## Objective
Map PD-L1 response to local chemical perturbations.

## Tested Variants
- GGYWPG (C-terminal Glycine)
- GGYWP  (Baseline)
- GGYFP  (Phe swap)
- GGYWPR (C-terminal Arg)

## Findings

### Flexibility Reward
- GGYWPG achieved 3.01 Å clearance
- Interpretation: Entropic tail relaxation improves docking

### Charge Veto
- GGYWPR collapsed to ~2.2–2.4 Å
- Interpretation: PD-L1 pocket intolerant of C-terminal bulk/charge

### Hydroxyl Mandate
- GGYFP collapsed (~2.3 Å)
- Interpretation: Tyrosine hydroxyl required for physical anchoring

## Laws Extracted
1. Tail flexibility is rewarded
2. C-terminal charge is vetoed
3. Tyrosine hydroxyl is mandatory

## Champion
- GGYWPG promoted for Tier 3 (10ns) verification
