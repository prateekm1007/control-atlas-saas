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
- [owner] **10. pLDDT scope LAW-130/155/200** — verified + pinned; owner must choose all-atoms vs core (evidence: `evidence/item-10/DECISION_REQUIRED.md`)

## AlphaFold acquisition

- [x] **11. Metadata-first lookup** — evidence: `evidence/item-11/` (live-verified: Q5VSL9/P69905 v6 cif; no model_v4, no data[0], JSON refused)
- [x] **12. Field names verified** — against the LIVE API (schema page unreachable from this env); fields_present logged per acquisition
- [x] **13. gemmi==0.6.7 in lock file** — mmCIF via gemmi CIF API, PDB-path-identical semantics
- [x] **14. Fail-closed cohort cases** — P0DTC2 migrated id AF-0000000365840314 recorded; P62988 UniProt Inactive recorded (cohort disposition: item 16 owner)

## Runner (PR-3)

- [x] **15. Runner hardening** — evidence: `evidence/item-15/` (no overwrite; content format; engine-commit provenance + dirty-tree gate; complete only at 50/50)

## Cohort (after 3–15)

- [owner] **16. Cohort rebuild** — proposal written: `backend/calibration/protocol/cohort_v25_proposal.md` (1I0T removal verified + recorded; P0DTC2/P62988 dispositions verified live; ALL numeric criteria + split fraction await owner). Held-out split NOT written — blocked on owner per direction

## Hygiene

- [x] **17. Nested copy removed** — branch `chore/hygiene-2026-10-19` (86 deletions, grep evidence: zero references)
- [x] **18. Junk files removed** — render tarball, test_sieve.sh, run_full_atlas.ps1 (deferred-with-reasons list in evidence)

## Owner gates (not started by coder)

- [owner] Rotate the Gemini key (still in main's HEAD)
- [owner] Set EC2 secrets (`TOSCANINI_API_KEYS`, `INTERNAL_CALLBACK_SECRET`, `ALLOWED_ORIGINS`, `GEMINI_API_KEY`)
- [owner] Branch protection on main; require the Toscanini Stability Gate
- [owner] Governance file location (Gate 1) and remaining constitution decisions (Gate 2)
- [owner] Approve any verdict-changing PR before merge to main
