# Item 4 evidence — file-order residue ordering

**Base SHA:** `296f10c` (science head) · **Program branch:** `science/p0-correctness-2026-10-19`
**Executed:** 2026-10-19, venv `/home/z/my-project/.venv-ci`
(install: `pip install -r backend/requirements-lock.txt -r backend/requirements-dev.txt`)

## What changed

| File | Change |
|---|---|
| `backend/tos/engine/coord_math.py` | `get_sequential_residues` no longer sorts by `(res_seq, icode)` — file order (first appearance) per chain; mmCIF `label_seq_id` takes precedence when atoms carry it. `residues_are_adjacent` rewritten as a no-missing-span predicate: chain equality + `seq_b - seq_a <= 1`; **alphabetical icode comparison removed**. New `peptide_bonded()` + `PEPTIDE_BOND_MAX_ANGSTROM = 2.0` (C–N chemical validation). |
| `backend/tos/engine/tier1_measurements.py` | `_extract` returns residues in file order grouped by chain first-appearance (sort removed). LAW-120 peptide angles, LAW-125 phi/psi, LAW-135 omega now require `peptide_bonded` — angles are never measured across a discontinuity. LAW-110/LAW-160 semantics unchanged (break detectors over no-missing-span pairs). |
| `backend/tos/ingestion/processor.py` | `Atom.label_seq_id: Optional[int] = None` (mmCIF ordering hook). |
| `backend/tests/test_residue_ordering.py` | 17 new tests (see below). |

## Acceptance criterion (audit direction item 4)

> a synthetic fixture with 183, 184A, 184, 185 gives zero breaks

`test_acceptance_insertion_before_base_zero_breaks` — **PASS**:
LAW-110 observed 0 / sample 3, LAW-160 observed 3.808 ≤ 4.5 / sample 3,
LAW-120 PASS with sample > 0. The fixture is built from Engh–Huber ideal
internal coordinates (planar extended peptide), so every C–N is exactly
1.329 Å in construction order.

## Test commands and results

```
$ cd backend && python -m pytest tests/test_residue_ordering.py -q --tb=short
17 passed in 0.32s

$ cd backend && python -m pytest tests -q --tb=short
452 passed, 5 skipped, 1209 warnings in 8.29s
```

(baseline before this item: 435 passed, 5 skipped — no regressions, 17 new)

## Raw 1PPE evidence (real deposited bytes)

- Source: `https://files.rcsb.org/download/1PPE.pdb`
- sha256 `fe51b4280055199c644214049f59a752e6124d858eaea29c1a2f721fb03ed0bb`, 203,958 bytes — **matches the digest stated in the audit direction**.
- Script: `evidence/item-4/1ppe_evidence.py` (verifies the digest before printing anything)
- Output: `evidence/item-4/1ppe_output.txt`
- ATOM lines for the insertion regions: `evidence/item-4/1ppe_chainE_insertion_ATOM_lines.txt` (189 lines)

Key facts printed from the real bytes:

1. Chain E file order around insertions is `183, 184A, 184, 185`,
   `187, 188A, 188, 189`, `220, 221A, 221, 222` — insertion-coded
   residues listed BEFORE their base number.
2. In file order every window pair is a genuine peptide bond
   (C–N 1.263–1.351 Å); Cα–Cα across insertion pairs 3.69–3.86 Å.
3. Under the removed `(res_seq, icode)` sort the SAME coordinates give
   9 fabricated gaps (C–N 2.894–5.808 Å) — the artifact behind v24's
   "1PPE genuine 7.12 Å coordinate break" conclusion (superseded, item 3).
4. `217 -> 219` (numbering gap, residue 218 unmodeled) is geometrically
   bonded (C–N 1.331 Å) but is still skipped by the missing-span rule —
   absent density is not measured across, by design.

## Verdict change on real 1PPE (engine-level, this commit)

| Law | Before (sorted order) | After (file order) |
|---|---|---|
| LAW-110 | gaps present around 184A/188A/221A → VETO cause | observed 0, sample 244, PASS |
| LAW-160 | 7.12 Å "genuine break" (v24 report) | observed 3.94 Å, sample 244, PASS |
| LAW-120 | pseudo-angles across mis-paired bonds | observed 2.63°, sample 740, PASS |
| Top-level verdict | VETO (LAW-110 cause, v24) | **PASS** |

Reproduction (engine on the pinned bytes):

```
$ cd backend && python - <<'EOF'
from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements
from tos.engine.adjudicator import adjudicate_laws, AdjudicationInput
data = open('/path/to/1PPE.pdb','rb').read()   # sha256 fe51b428...
s = IngestionProcessor.run(data, '1PPE.pdb', 'Audit', 'experimental')
raw, cov, _ = Tier1Measurements.run_full_audit(s)
res = adjudicate_laws(AdjudicationInput(raw_measurements=raw, coverage=cov,
    is_experimental=True, method_type='xray', architecture='NONE'))
print(res.verdict)   # PASS
EOF
```

## Notes and limits

- The two insertion-code test files the auditor had locally do not
  contain this pattern and therefore could not test it; 1PPE (above) and
  the synthetic fixture do.
- mmCIF ordering by `label_seq_id` is implemented and unit-tested with
  synthetic `Atom` records; a real mmCIF path through the parser arrives
  with item 13 (maintained mmCIF parser).
- PDB files whose chains interleave are grouped by chain first
  appearance; within a chain, file order always wins.
