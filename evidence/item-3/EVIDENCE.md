# Item 3 evidence — v24 superseded

Executed: 2026-10-19T11:10:05Z

## Commands and outputs
```
$ git show 0eae6d5:backend/tos/engine/tier1_measurements.py | sed -n '150,170p'
            # Intra-residue: O-C-N (requires O in same residue, N is same residue's N)
            if all(a in atoms for a in ["O", "C", "N"]):
                measured = angle_deg(atoms["O"], atoms["C"], atoms["N"])
                ideal = IDEAL_TABLE["O-C-N"]
-> commit recorded in v24 run_manifest.json (code_sha 0eae6d575cd9...) contains the
   pre-repair intra-residue pseudo-angle O-C-N(i).

$ git show 863e6dd:backend/tos/engine/tier1_measurements.py | grep -c "r2\[\"_atoms\"\]\[\"N\"\]"
4
-> repaired O(i)-C(i)-N(i+1) term first exists at 863e6dd, AFTER the recorded code_sha.

$ python3 -c "import json;m=json.load(open('backend/calibration/runs/v24/run_manifest.json'));print(m[('code_sha')],m[('backend_version')],m[('expected_target_count')])"
0eae6d575cd91376588fe6b36d1b9212352df7fc 22.5.3 50

$ sha256sum backend/calibration/runs/v24/run_manifest.json backend/calibration/runs/v24/calibration_results.json
f602351f4844531d409cfc256939508103e3c2d4f57fe21d7f77c34494aa4036  backend/calibration/runs/v24/run_manifest.json
fd1cc8a62baeee03094b9f716ebcbc04d33a3ccfe15abaa21ba8275b7ea21bb8  backend/calibration/runs/v24/calibration_results.json
-> JSON artifacts NOT modified (marker file added beside them only).
```

## Files changed
- added: backend/calibration/runs/v24/calibration_v24_SUPERSEDED.md
- added: docs/P0_STATUS.md
- added: evidence/reading-note/READING_NOTE.md (item 1)
- added: evidence/item-3/EVIDENCE.md (this file)

## Verdict changes
- None (documentation only; no code path touched).
