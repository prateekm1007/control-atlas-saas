# Per-law inventory (P0 item 8, completing ticket v2 PR-2 item 2)

**Generated from:** `backend/tos/governance/station_sop.py` (LAW_CANON,
LAW_METHOD_CLASSIFICATIONS), `backend/tos/engine/tier1_measurements.py`,
`backend/tos/engine/adjudicator.py`, `governance/law_registry.yaml` at
commit `dfb9735`. This is the engine-side inventory as implemented, not
as documented.

Canonical method is shown; the modality matrix
(`governance/modality_matrix.json`) may resolve a law to
`advisory_experimental` per structure source (e.g. LAW-100 for X-ray),
which downgrades VETO to advisory FAIL.

| Law | Title | Method | Measures | Sample unit | Empty-sample status | Atom scope | Inter-residue pairing |
|---|---|---|---|---|---|---|---|
| LAW-100 | Bond Integrity | deterministic | RMSZ of covalent bond lengths vs Engh–Huber (exp) / % outliers z>4σ (pred) | bonds measured | NOT_MEASURED | core | intra-residue |
| LAW-105 | Reliability Coverage | deterministic | % residues in core (pLDDT ≥ 70) over **source** residues | residues (sample = parsed; denominator = source) | FAIL at coverage < 70; denominator ≥ 1 | all (denominator counts source) | n/a |
| LAW-110 | Backbone Gap | deterministic | count of C(i)–N(i+1) > 2.0 Å | **observed peptide bonds** (file-order, no missing span) | NOT_MEASURED | core | file order + no-missing-span |
| LAW-120 | Bond Angle | deterministic | RMSD of N-CA-C, O-C-N, CA-C-N vs ideals | angles (intra + bonded peptide pairs) | NOT_MEASURED | core | file order + `peptide_bonded` (C–N ≤ 2.0 Å) |
| LAW-125 | Ramachandran | deterministic | % phi/psi outliers | residue triples | NOT_MEASURED | core | file order + bonded both sides |
| LAW-130 | Clashscore | deterministic | non-bonded heavy-atom pairs < 1.5 Å per 1000 atoms | atoms (all) | NOT_MEASURED / MEASUREMENT_FAILED | **all atoms** (see item 10) | neighbor exclusion by (chain, \|Δseq\| ≤ 1) |
| LAW-135 | Omega Planarity | deterministic | % omega > 30° from planar ideal | bonded peptide pairs | NOT_MEASURED | core | file order + bonded |
| LAW-145 | Chirality | deterministic | count of D-amino improper dihedrals (0–90°) | residues with N/CA/C/CB | NOT_MEASURED | core | intra-residue |
| LAW-150 | Rotamer Audit | deterministic | % chi1 outside rotamer wells | side chains with chi1 | NOT_MEASURED | core | intra-residue |
| LAW-155 | Voxel Occupancy | heuristic | occupied 2 Å voxels per residue | atoms (all) | NOT_MEASURED | **all atoms** (see item 10) | n/a |
| LAW-160 | Chain Integrity | deterministic | max Cα–Cα over sequential pairs | adjacent CA pairs (file order, no missing span) | NOT_MEASURED | core via get_atom_pos | file order + no-missing-span |
| LAW-170 | Residue Identity | deterministic | count of non-standard residues | all parsed residues | n/a (always measurable) | all residues | n/a |
| LAW-182 | Hydrophobic Burial | heuristic | min per-chain fraction of hydrophobic CAs within 15 Å of centroid | hydrophobic CAs | NOT_APPLICABLE (none present) / NOT_MEASURED (present, unmeasurable) | all residues | n/a |
| LAW-195 | Disulfide Geometry | deterministic | max deviation of SG–SG from 2.033 Å | CYS SG pairs < 3.0 Å | NOT_APPLICABLE | all residues | n/a |
| LAW-200 | Packing Quality | heuristic | bounding-box volume per atom | atoms (all) | NOT_MEASURED | **all atoms** (see item 10) | n/a |

Contained law: **LAW-120** only (`LAW_CONTAINMENT` status contained) — a
measured LAW-120 violation cannot independently issue an authoritative
VETO.

Retired laws (registry `status: retired`, never measured by pinned
code): LAW-190 Ring Planarity, LAW-210 Helix Geometry, LAW-220 Strand
Geometry, LAW-230 Valence.

## Status semantics (adjudicator, post item 8)

- Absent law (engine emitted no row) → `NOT_MEASURED` with reason
  "no measurement row emitted by the engine" (WP-SCI-1, retained).
- Present-but-statusless row (malformed) → `NOT_MEASURED` (item 8;
  previously defaulted to `FAIL`).
- Both paths remain fail-closed: NOT_MEASURED on a deterministic law
  enters `indeterminate_causes` and blocks an authoritative PASS.

## LAW-110 counting (item 8)

LAW-110's sample is the number of peptide bonds actually observed
(`peptide_pairs_checked`), never `len(residues) - 1`. Asserted by:
`tests/test_residue_ordering.py` (samples 3 and 2 on fixtures with 4
residues each) and the degenerate shapes in
`tests/test_honest_statuses_coverage.py`.
