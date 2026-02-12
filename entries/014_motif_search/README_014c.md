# Entry #014c: Composite Motif Search — Switch II

## Intent
Identify proteins sharing the Switch II control geometry independent of sequence.

## Probe
- M-001 Sub-motif B
- Switch II (17 residues, continuous)

## Method
- Sliding window over Cα backbone
- Rigid-body superposition (Kabsch)
- RMSD threshold: < 2.5 Å

## Expected Validation
- 5P21 (Ras): MATCH
- 2RGN (RhoA): MATCH
- 1UBQ (Ubiquitin): REJECT

## Significance
Demonstrates geometry-first discovery of druggable control surfaces.

