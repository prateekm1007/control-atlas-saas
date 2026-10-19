# Item 11-14 evidence — AlphaFold acquisition repair

**Base SHA:** 07289bc · **Executed:** 2026-10-19T11:57:06Z

## Item 12 — field names verified against the LIVE API (2026-10-19)
The official schema page could not be fetched from this environment either; the live
API itself was queried instead (stronger evidence):
```
$ curl -s https://alphafold.ebi.ac.uk/api/prediction/Q5VSL9 | python -m json.tool | grep -E "modelEntityId|latestVersion|cifUrl|pdbUrl"
   "modelEntityId": "AF-Q5VSL9-F1",  "latestVersion": 6,  "allVersions": [1..6],
   "cifUrl":   "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v6.cif",
   "pdbUrl":   "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v6.pdb"
```
Confirmed coexistence of pdbUrl and cifUrl (the sunset-conflict note in the audit);
the module codes against cifUrl only, and logs fields_present per acquisition.

## Item 11 — metadata-first acquisition
- tos/generation/afdb.py: /api/prediction/{uniprot} -> select by modelEntityId
  (NEVER data[0]) -> max latestVersion -> cifUrl -> content-based format detection
  (JSON refused as a structure) -> AcquiredPrediction with URL/version/SHA-256/
  field names. Every failure raises AFDBUnavailable with a structured record.
- dispatcher._acquire_afdb rewired; hardcoded model_v4 URLs gone (the URLs that
  404'd for all 10 AF targets in v24).

## Item 13 — maintained mmCIF parser
- gemmi==0.6.7 pinned in requirements-lock.txt and requirements.txt.
- IngestionProcessor auto-detects format from content; mmCIF parsed via gemmi's
  CIF table API (no in-house parser), mirroring PDB-path semantics exactly:
  ATOM/HETATM-mapped dispatch, H/D exclusion, one-conformer-per-atom, MSE->MET, echo   label_seq_id ordering key (item 4 hook), source_residue_count, _exptl
  classification.

## Item 14 — fail-closed cohort cases (live-verified)
- AF-P0DTC2-F1: API serves modelEntityId AF-0000000365840314 (provider VR3D,
  ColabFold v1.5.2), not AF-P0DTC2-F1. Requesting the canonical ID fails closed
  with migrated_identifier recorded; UniProt entry itself still 'reviewed'.
  OWNER NOTE: whether the v25 cohort follows the migrated entity or replaces
  the target is a cohort decision (item 16).
- AF-P62988-F1: API HTTP 404; UniProt entryType 'Inactive' (merged entry).
  Recorded via uniprot_entry_type() on the acquisition record; removal from
  the cohort is justified by this status (final call: item 16, owner).

## Live smoke test (real network, real code path)
script: evidence/item-11/live_smoke_test.py — output: live_smoke_output.txt
```
AF-Q5VSL9-F1  ACQUIRED  v6  cif  776095 bytes  sha256 eb8fb805...
AF-P0DTC2-F1  FAIL-CLOSED  migrated_identifier=AF-0000000365840314
AF-P62988-F1  FAIL-CLOSED  HTTP 404  uniprot_entry_type=Inactive
AF-P69905-F1  ACQUIRED  v6  cif  137559 bytes
```

## End-to-end (acquisition -> gemmi CIF -> engine -> verdict)
evidence/item-11/end_to_end_measure.txt: AF-P69905-F1 verdict PASS,
coverage 99.3, LAW-120 1.59 deg, LAW-160 3.9 A, 1077 atoms.

## Unit tests (no network; injected transport)
```
$ cd backend && python -m pytest tests/test_afdb_acquisition.py -q
16 passed
$ cd backend && python -m pytest tests -q
506 passed, 5 skipped, 0 failed
```
