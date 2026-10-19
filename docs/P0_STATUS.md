# P0 Program Status (items 3–18 of the 2026-10-19 audit directions)

**Work branch:** `science/p0-correctness-2026-10-19` (from science head `296f10c`)
**Reading note:** `evidence/reading-note/READING_NOTE.md` (item 1 — done)
**PAT discipline:** item 2 — no token written to any file, commit, log, or
remote URL in this program. Push uses a runtime environment variable only.

Legend: `[x]` done, evidence linked · `[ ]` open · `[owner]` blocked on owner

## Engine correctness

- [x] **3. v24 marked superseded** — `backend/calibration/runs/v24/calibration_v24_SUPERSEDED.md` (JSON untouched). Evidence: `evidence/item-3/`
- [x] **4. File-order residue ordering** — evidence: `evidence/item-4/` (17 new tests; 452/5/0 full suite; 1PPE raw-bytes proof committed)
- [x] **5. HETATM retention incl. MSE→parent** — evidence: `evidence/item-5/` (1A8O chain A = 70 residues on pinned real bytes; 459/5/0 full suite)
- [x] **6. Altloc selection** — evidence: `evidence/item-6/` (real 3ERD fixture; 8 tests; 467/5/0 full suite)
- [x] **7. Dead `_law_155` removed** + uniqueness test — evidence: `evidence/item-7/` (472/5/0 × 3 runs; credits isolation defect also fixed)
- [x] **8. Honest statuses complete** — evidence: `evidence/item-8/` (inventory + adjudicator defaults; 481/5/0)
- [x] **9. Coverage denominator from source residues** — evidence: `evidence/item-9/` (incl. HETATM-mapped)
- [ ] **10. pLDDT scope verified for LAW-130/155/200** — evidence recorded; semantics change is an owner decision

## AlphaFold acquisition

- [ ] **11. Metadata-first lookup** (`/api/prediction/{uniprot}`, select by `modelEntityId` + `latestVersion`, `cifUrl`); no `data[0]`; content-based format detection; URL/version/SHA-256 recorded
- [ ] **12. Field names verified against official schema** (`alphafold.ebi.ac.uk/api-doc`)
- [ ] **13. Maintained mmCIF parser in lock file**
- [ ] **14. AF-P0DTC2-F1 fail closed + migrated identifier recorded; AF-P62988-F1 UniProt status checked**

## Runner (PR-3)

- [ ] **15. Runner hardening** — fail if run dir exists; content-based format detection (no hardcoded `origin.pdb`); record executing backend's commit; exit non-zero on any failure; `complete` only at 50/50

## Cohort (after 3–15)

- [ ] **16. Cohort rebuild per `protocol/evaluation_protocol_v1.md`** — 1I0T removed with recorded reason; PISCES-style quality filter; held-out split before threshold work. **[owner] all numeric criteria**

## Hygiene

- [ ] **17. Nested `control-atlas-saas/` removed in its own PR** (83 tracked files) + grep evidence
- [ ] **18. Junk files removed**

## Owner gates (not started by coder)

- [owner] Rotate the Gemini key (still in main's HEAD)
- [owner] Set EC2 secrets (`TOSCANINI_API_KEYS`, `INTERNAL_CALLBACK_SECRET`, `ALLOWED_ORIGINS`, `GEMINI_API_KEY`)
- [owner] Branch protection on main; require the Toscanini Stability Gate
- [owner] Governance file location (Gate 1) and remaining constitution decisions (Gate 2)
- [owner] Approve any verdict-changing PR before merge to main
