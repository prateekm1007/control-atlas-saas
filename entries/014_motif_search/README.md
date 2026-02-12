# Entry #014 — Motif Search (RMSD-Based Discovery)

## Objective
Find other proteins containing the same control geometry as M001_Switch_GTPase.

## Core Question
> "Where else in biology does this exact switch geometry exist?"

## Method
- Reference: `library/motifs/M001_Switch_GTPase.pdb`
- Search space: GTPase superfamily (Rab, Ran, Arf, Rho)
- Metric: Backbone RMSD ≤ 1.5 Å = strong match

## Outputs
- `library/matches/M001_hits.csv`
- `Atlas_Entry_014_Motif_Search.png` (overlay of top hits)

## Status
[ ] Candidate structures pulled
[ ] RMSD calculations complete
[ ] Hits catalogued
[ ] Proof image rendered
