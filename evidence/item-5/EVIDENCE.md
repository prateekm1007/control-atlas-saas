# Item 5 evidence — HETATM retention incl. MSE mapped to parent

**Base SHA:** 044de69 (item 4 commit) · **Executed:** 2026-10-19T11:20:55Z

## Acceptance (audit direction item 5)
> 1A8O chain A yields 70 residues. The ATOM-only parser yields 66.

```
$ sha256sum backend/tests/fixtures/real_pdb_1a8o.pdb
b0b4bf0d52c313047d5d7a7b188bf07a66a57c480db8ec4680f29461f18cbe4e  backend/tests/fixtures/real_pdb_1a8o.pdb

$ grep -c "^ATOM" backend/tests/fixtures/real_pdb_1a8o.pdb  (ATOM lines)
524
$ grep "^HETATM" backend/tests/fixtures/real_pdb_1a8o.pdb | grep -c MSE  (MSE atom lines)
32
$ grep "^HETATM" backend/tests/fixtures/real_pdb_1a8o.pdb | grep MSE | cut -c23-26 | sort -u  (MSE residues)
 151
 185
 214
 215
-> 66 ATOM residues + 4 MSE HETATM residues = 70
```

## Test commands and results
```
$ cd backend && python -m pytest tests/test_hetatm_retention.py -q
7 passed in 0.11s

$ cd backend && python -m pytest tests -q
459 passed, 5 skipped, 1209 warnings in 8.24s
```

## Files changed
- backend/tos/ingestion/processor.py: MODIFIED_RESIDUE_PARENTS map + HETATM dispatch (MSE->MET etc.; HOH/ligands/ions excluded; D-amino acids deliberately unmapped)
- backend/tests/fixtures/real_pdb_1a8o.pdb: pinned real bytes (sha256 b0b4bf0d...)
- backend/tests/test_hetatm_retention.py: 7 tests

## Verdict changes
- Structures with MSE residues no longer lose them: backbone laws now measure across
  retained residues (e.g. 1A8O: LAW-110 sample grows by the MSE-adjacent bonds).
- No law threshold or meaning changed.
