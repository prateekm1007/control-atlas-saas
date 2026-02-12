# Entry #021 — Library Filtering

## Prerequisites
- Entry 020 LOCKED (Gatekeeper)

## Objective
Apply the Gatekeeper to filter large compound libraries before docking or synthesis.

## Core Question
> "Which compounds in this library are worth testing against KRAS G12C?"

## Input
- SMILES file (one compound per line)
- Format: `SMILES,ID` or just `SMILES`

## Process
1. Parse each SMILES
2. Apply Entry 020 validation (Grammar + Warhead + Quantitative)
3. Classify as PASS or REJECT
4. Log rejection reasons

## Output
- `passed_candidates.csv` — Compounds ready for screening
- `rejected_candidates.csv` — Failed compounds with reasons
- `filtering_stats.txt` — Summary statistics

## Success Criteria
- [ ] Processes 1000+ compounds without error
- [ ] Known actives (Sotorasib-like) pass
- [ ] Decoys rejected with correct reasons
- [ ] Rejection rate matches expected (~90%+ for random libraries)

## Status
[ ] Filter script implemented
[ ] Test library processed
[ ] Stats generated
[ ] Committed
