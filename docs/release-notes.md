# Toscanini Release Notes

## v22.5.4-A.5-final-polish (2025-01-XX)

### Phase A.5 Complete: Surgical Residue-Level Diagnostics

**Core Features:**
- Chain-aware residue capture for 5 deterministic geometry laws:
  - LAW-125 (Ramachandran): residue-level outliers
  - LAW-150 (Rotamer): side-chain outliers
  - LAW-130 (Clashscore): residue-pair clashes
  - LAW-135 (Omega): peptide bond planarity outliers
  - LAW-145 (Chirality): D-amino acid violations

**PDF Enhancements:**
- Page 7: Unified "Residue-Level Diagnostics" header
- All 5 laws show specific residue lists with chain prefixes (e.g., A:12)
- Overflow truncation: max 10-20 items per law with "(+ N more)" notation
- Institutional formatting and explanatory text

**Remediation Package:**
- XML comments for all 5 laws in rosetta_relax.xml header
- Chain-aware residue lists guide Rosetta selector targeting
- Overflow logic enforced (max 10-15 residues per law)
- Deterministic ordering: sorted by chain, then residue number

**Technical:**
- Adjudicator propagates granularity + failing_residues + failing_residue_pairs
- Schema integrity: sample_size fix (no zero denominators)
- No governance drift (canon hash unchanged, thresholds unchanged)
- Determinism preserved throughout (ZIP byte-identical for same audit)

**Validation:**
- Tested with 4HHB (hemoglobin), AF-P01308-F1 (insulin), synthetic clash structures
- Overflow tested with 229 clash pairs → correctly truncates to 10 + "(+ 219 more)"
- All residue formats use chain:seq notation (A:3, B:12, etc.)

**Next Phase:**
Phase B (managed refinement orchestration) planning begins.

---

## v22.5.3 (Previous Release)
Phase A complete: Theater removal, deterministic remediation packages, before/after comparison.
