# ⚠️ v24 RUN SUPERSEDED — DO NOT CITE

**Status:** SUPERSEDED (2026-10-19). This entire run directory is retained
for audit history only. None of its numbers may be cited as evidence of
engine behavior, and none of its conclusions may be relied on.

## Why

The run manifest (`run_manifest.json`, unmodified beside this file) records:

- `"code_sha": "0eae6d575cd91376588fe6b36d1b9212352df7fc"` — the head of
  `integration/p0-hardening-2026-10-19`, which is the PARENT of the
  WP-SCI-1 science commits (`863e6dd`, `296f10c`).
- `"backend_version": "22.5.3"`.

Commit `0eae6d5` **contains the pre-repair LAW-120 measurement** — the
intra-residue pseudo-angle `angle_deg(atoms["O"], atoms["C"], atoms["N"])`
at `backend/tos/engine/tier1_measurements.py` in that revision — while this
run's report describes "corrected metrics" that only exist from `863e6dd`
onward. Verified by:

```bash
git show 0eae6d5:backend/tos/engine/tier1_measurements.py | grep -n "angle_deg(atoms\[\"O\"\], atoms\[\"C\"\], atoms\[\"N\"\])"
# present (pre-repair term)
git show 863e6dd:backend/tos/engine/tier1_measurements.py | grep -n "r2\[\"_atoms\"\]\[\"N\"\]"
# the repaired O(i)-C(i)-N(i+1) term first appears here
```

The run therefore either executed code that was not what `code_sha`
records (uncommitted working-tree state), or executed the recorded commit
(pre-repair). Either way the provenance binding required of an evidence run
is broken: the report's numbers cannot be tied to a specific, reproducible
code state.

Additionally:

- `run_status: incomplete` — 40 of 50 targets measured; all 10 AlphaFold
  acquisitions failed (HTTP 404 on the hardcoded `model_v4` file URLs).
- The report's top-level verdict counts (33 PASS / 6 VETO / 1
  INDETERMINATE) and its per-law readings (including the "1PPE genuine
  7.12 Å coordinate break" conclusion) are **not citable**: the 1PPE
  discontinuity was later shown to be an artifact of the engine's
  `(res_seq, icode)` sort reordering insertion-coded residues that the
  file lists before their base number, and the measurement provenance
  itself is the open defect documented above.

## What supersedes it

A v25 run must be executed under `backend/calibration/protocol/
evaluation_protocol_v1.md` after the P0 program closes items 3–15
(`docs/P0_STATUS.md`), with:

- the executing backend's commit recorded from the code actually imported
  at measurement time (plus working-tree cleanliness), and
- acquisition through the metadata-first AlphaFold path (no hardcoded
  `model_v4` URLs), and
- 50 of 50 targets measured before `run_status` may be `complete`.

Until that run exists, no corrected-metric distribution claim from v24
should be repeated in reports, papers, or READMEs.
