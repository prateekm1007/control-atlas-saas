# Evaluation Protocol v1.1 Addendum — Engine Audit Reconciliation

**Status:** REVIEWED ADDENDUM — reconciles v1 protocol with the 2026-10-23 engine audit.
**Effective:** 2026-10-23
**Parent protocol:** `evaluation_protocol_v1.md` (ACTIVE, `calibration_status=not_approved`, `supports_threshold_validation=false`)
**Scope:** Reconciliation of protocol claims with the engine correctness audit results.

---

## 0. What this addendum is and is not

- This addendum does NOT modify, supersede, or replace the v1 protocol. The v1 protocol remains the authoritative historical record.
- This addendum reconciles the protocol's claims with the engine audit findings and records what has changed.
- This addendum does NOT authorize any threshold change, cohort modification, or held-out evaluation.
- OD-16 remains in force: no v25 cohort execution, no held-out evaluation, no threshold derivation.

## 1. Engine audit summary

The engine correctness audit (branch `repair/od-scope-v25-engine-audit`, commit `05795f2`) implemented and verified the following:

### 1.1 Canonical ingestion
- Gemmi-based PDB and mmCIF parsing with identical coordinate extraction
- First-model selection with `excluded_models` accounting
- HETATM polymer residues included; waters, non-polymer ligands, and other HETATM excluded with counters
- Alternate-location selection by highest occupancy, tie-break blank then lexicographic
- Residue-level altloc selection (prevents artificial mixed-conformer residues)
- Insertion-code preservation through all downstream paths
- Malformed coordinates, empty inputs, and unsupported formats raise `IngestionError`
- `IngestionError` propagates through API (`HTTP 422`) and calibration runner (`MEASUREMENT_FAILED`)

### 1.2 OD-07 canonical identity
- Coordinate hash: canonical JSON encoding, sorted atom records, rounded coordinates, atom/residue identity
- Atom-order invariant (verified by test)
- Audit ID binds coordinate hash + `LAW_CANON_HASH` + modality-matrix hash via canonical JSON
- Component-hash sensitivity verified (each component independently changes the ID)
- Used consistently across API response, audit store, NKG, and PDF

### 1.3 Coverage
- Experimental coverage = complete-backbone residues (N, CA, C) / total residues (no fabricated 100%)
- Predicted coverage = complete-backbone residues with pLDDT >= 70 / total residues
- `LAW-105` exposes `coverage_basis`, `coverage_numerator`, `coverage_denominator`, `required_atoms`, `complete_backbone_residues`, `source_atom_records`, `selected_atoms`, `discarded_atoms`

### 1.4 Geometry
- `_is_sequential()` requires immediate insertion-code successor for same sequence; sequence+1 allows only blank→blank or insertion-code→blank
- `_has_peptide_bond()` additionally requires finite C–N distance in [0.5, 2.0]
- LAW-120: peptide-bond gating for inter-residue angles
- LAW-125: peptide-bond gating for phi/psi triples
- LAW-135: peptide-bond gating for omega
- LAW-160: uses continuity only (not C–N bond) to preserve break detection
- LAW-130: clash analysis with insertion-code-aware same-residue exclusion and canonical atom scope
- LAW-145/150/182: fail-closed on missing/invalid geometry (no favorable zeros)
- Degenerate geometry: all tier1 paths use `_valid` variants; invalid observations excluded from counts

### 1.5 Confidence provenance
- Experimental B-factors never labeled as pLDDT
- Unknown acquisition method returns `"unknown"` (treated as experimental, fail-closed)
- pLDDT observations are finite-filtered

### 1.6 NKG
- Cross-platform file locking (fcntl on Unix, msvcrt on Windows)
- Atomic write with fsync

## 2. What has NOT changed

- No thresholds modified
- No law definitions modified
- No cohort membership changed
- No golden artifacts regenerated
- No v25 cohort executed
- No held-out evaluation performed
- `calibration_status=not_approved` and `supports_threshold_validation=false` remain

## 3. Test evidence

- 570 tests passed, 5 skipped (gate: 548) under `PYTHONUTF8=1`
- 21 CP1252 failures are pre-existing (verified by stashing on base commit)
- Invariance tests: rotation, translation, atom reordering, chain relabeling, determinism, PDB/mmCIF parity
- Negative tests: altloc edge cases, ingestion errors, unknown acquisition, degenerate geometry, LAW-182 missing chain, coordinate hash sensitivity

## 4. Remaining gaps

- External engine validation (≥200 structures, MolProbity ≥95% agreement) — harness not yet built
- Append-only historical evidence repair (22 missing refs) — not yet implemented
- P-091-CONFIRM — requires actual CI workflow evidence from pushed SHA
- P-092-INTENT — terminal intent, to be introduced by final commit
- Protocol addendum does not authorize any new calibration evidence generation

## 5. Reconciliation with roadmap

This addendum reconciles the v1 protocol with the engine audit by:
1. Recording what engine contracts have been implemented and verified
2. Confirming no threshold, law, or cohort changes occurred
3. Documenting the test evidence
4. Preserving the v1 protocol's `not_approved` / `supports_threshold_validation=false` status
5. Not authorizing any new calibration evidence generation

OD-16 remains in force.
