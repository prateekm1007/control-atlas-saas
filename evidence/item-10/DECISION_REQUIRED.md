# DECISION REQUIRED (owner): pLDDT scope of LAW-130 / LAW-155 / LAW-200

**Status:** OPEN — owner gate. No code change made to the laws' scope.

## Verified current behavior (commit `d5fa87a`)

For **predicted** structures, `Tier1Measurements.run_full_audit` computes:

- `core` = residues with mean pLDDT ≥ 70 (LAW-105 coverage basis)
- LAW-100/110/120/125/135/145/150/160 measure over `core`
- **LAW-130 (Clashscore), LAW-155 (Voxel Occupancy), LAW-200 (Packing
  Quality) measure over ALL atoms of the structure**, regardless of
  per-residue confidence.

Code references (`backend/tos/engine/tier1_measurements.py`):
`all_coords = np.array([a.pos for a in structure.atoms])` feeds all
three laws; the `core` list is never applied to them.

## Why the choice matters (demonstrated)

`evidence/item-10/scope_demo.py` (output in `scope_demo_output.txt`)
builds an AlphaFold-like structure: a compact 20-residue core at
pLDDT 90 plus a 10-residue extended tail at pLDDT 40 — the shape
low-confidence AF regions genuinely take.

| Law | All atoms (current) | Core only (alternative) |
|---|---|---|
| LAW-155 Voxel Occupancy | 2.57 → PASS | 3.35 → PASS |
| LAW-200 Packing Quality | **477.01 → VETO** | **68.09 → PASS** |

LAW-200's verdict flips. This is the same mechanism as the recorded
AF-tau case (771 Å³/atom against the 300 threshold): an
intrinsically-disordered low-confidence tail vetoes the whole
prediction under the current scope.

## The two defensible positions

**A. Keep all-atoms scope (current behavior).**
Rationale: clashscore, voxel sparsity, and packing describe the
physical artifact as delivered — a deployed structure with an exploded
tail IS an exploded structure regardless of model confidence; a
customer receives those atoms. Low-confidence regions are still
shipped, still clash, still occupy volume.

**B. Apply the core filter consistently.**
Rationale: every other law already restricts to the confident core;
scoring uncertain coordinates against thresholds derived from
well-validated structures punishes honest disorder (IDRs) rather than
detecting error. Under B, LAW-130/155/200 become statements about the
confident domain only, and verdicts on multi-domain or disordered
proteins change materially (see table above).

## What the coder needs from the owner

Pick A or B (or a per-law mix). Any change from A is
verdict-changing and will be implemented as its own PR with a
verdict-change table per structure class, under the constitution's
scientific-semantics change discipline (EC-4.7, EC-13). Until then the
current scope is pinned by test
(`backend/tests/test_plddt_scope_pin.py`) so it cannot drift silently.

Note: under either option, the thresholds themselves (2.0 / 300 / 20)
remain unvalidated-by-protocol — threshold derivation is gated on the
cohort rebuild (item 16) regardless of this decision.
