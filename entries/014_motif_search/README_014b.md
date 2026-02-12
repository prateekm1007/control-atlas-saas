# Entry #014b: Geometry-Only Motif Search

## Intent
Implement a blind, sequence-agnostic geometry engine that scans proteins
for structural matches to a control motif using RMSD after rigid-body alignment.

## Method
- Sliding window over Cα coordinates
- Kabsch superposition (ProDy `superpose`)
- RMSD computed on aligned coordinates only
- No sequence alignment, no annotations

## Role in Atlas
- Validates the geometric search engine
- Precursor to composite motif logic (Entry 014c)

## Output
- library/matches/M001_geometry_hits_014b.csv
