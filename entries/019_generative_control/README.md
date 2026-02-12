# Entry #019 — Scaffold Space Definition (Generative Control)

## Prerequisites
- Entry 018 LOCKED (Chemistry Grammar v2)

## Objective
Define the "legal" chemical scaffold space permitted by the Pocket Physics + Chemistry Grammar rules.

## Core Question
> "What abstract chemotypes are chemically legal for this control pocket?"

## Method (Combinatorial Assembly)
1. **Define Fragment Classes** based on Entry 018 regions:
   - Core (Buried)
   - Anchor (Concave)
   - Warhead (Covalent)
   - Tail (Solvent)
2. **Combinatorial Generation**:
   - Assemble `Core-Anchor-Warhead-Tail` combinations
3. **Grammar Filtering**:
   - Apply Entry 018b Rules (Quantitative, Resistance, Forbidden)
4. **Identify Valid Scaffold Classes**

## Fragment Source (Virtual Library)
- Cores: Quinazoline, Pyridopyrimidine, Indole, Purine, Pyridine
- Anchors: Isopropyl, Methyl, Fluorine, Phenyl (resistance test)
- Warheads: Acrylamide, Vinyl Sulfonamide, Propyl (negative control)
- Tails: Piperazine, Morpholine, t-Butyl (hydrophobic control)

## Success Criteria
- [ ] Valid scaffolds identified (e.g., Pyridopyrimidine-Isopropyl-Acrylamide)
- [ ] Invalid scaffolds rejected (e.g., Phenyl anchor at Gly60 clash)
- [ ] Sotorasib scaffold recovered
- [ ] Novel valid scaffolds discovered

## Output
- `legal_scaffolds.csv`
- `scaffold_space_stats.txt`

## Status
[ ] Fragment library defined
[ ] Assembler implemented
[ ] Space generation complete
[ ] Valid scaffolds verified
