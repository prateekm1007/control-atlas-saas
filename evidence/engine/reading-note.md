# Constitution reading note - `engine/od-scope-v25`

**Status:** Preflight record for OD-03/04/05/20 implementation.
**Base:** `c78ef15f696a2b1e2658e7ac7b54461f70c6e040` on `engine/od-scope-v25`.
**Fresh fetch:** `origin/main` verified `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` at `2026-10-21T18:40:02Z`; re-fetched `2026-10-21T18:41:00Z`.
**Complete read / hash verification:** `2026-10-21T18:42:00Z`.
**First edit:** `2026-10-21T18:45:00Z` (this note).

## Governing sources read

| Source | Blob on `origin/main` | Notes |
|---|---|---|
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | Mandatory entry point; EC/OD hard stop. |
| `constitution.md` (at base) | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | OBS-153 and CRIT-154 are ACTIVE. |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | EC sec. 1-15 and Final Rule. |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | Execution rules; adds no authority. |
| `docs/decisions/OWNER_DECISIONS.md` | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` | Register rev 9; controlling text. |
| `docs/P0_STATUS.md` | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` | Status context only. |

## Rules enforced

- **EC sec. 2:** base SHA, branch, working-tree state, full constitution read, call-site map, and machine-readable preflight before edits.
- **EC sec. 3:** no credential in the tree, logs, evidence, or push URL. Secret-scan before push.
- **EC sec. 4:** distinguish measured, inferred, missing, heuristic, and advisory states. Do not convert absent/low-confidence measurements into PASS. Do not silently redefine a law ID; keep conflicts explicit.
- **EC sec. 5:** inspect duplicate implementations and alternate entry points before deletion. This directly applies to the tracked nested `control-atlas-saas/` copy.
- **EC sec. 6:** failing regression test before the fix, exact commands/SHA/counts, deterministic tests, and no skipped/xfail workaround.
- **EC sec. 7:** feature branch only; no force-push; hash-lock gate remains authoritative.
- **EC sec. 12-14:** evidence bound to the candidate, security/scientific impact note, rollback, changed-file manifest, and an honest blocked/unverified list.
- **OBS-153 / CRIT-154:** read and preserved. OBS-153 avidity does not govern pLDDT scope. CRIT-154 is addressed below as a pre-existing conflict, not silently changed.

## Canonical register vs. owner clarification

The user explicitly directed that the task brief governs these IDs, while keeping the existing IDs and creating no OD-22...OD-25. This note records that direction because the main register text is materially different:

- **OD-03 (register):** pLDDT advisory only for LAW-130/LAW-200; retire LAW-155; no new pLDDT threshold.
- **OD-03 (owner clarification):** core pLDDT `>=70` for LAW-100/105/110/120/125/130/135/145/150/155/200; experimental keeps all atoms; LAW-155/200 remain heuristic and cannot independently veto. LAW-130 remains deterministic per the modality matrix.
- **OD-04 (register):** apply OD-03 to disordered proteins; no invented disorder threshold.
- **OD-04 (owner clarification):** one cutoff at `70`; delete the old `50` boundary; a structure with `<70%` of **atoms** at pLDDT `>=70` is `INDETERMINATE`, not `PASS`.
- **OD-05 (register and owner):** retain LAW-190/210/220/230 as retired, not enforced; no implicit PASS. The owner explicitly retained the IDs; the brief's deletion is not applied.
- **OD-20 (register):** keep NMR separate from resolution-based strata; no pooling.
- **OD-20 (owner clarification):** LAW-182 is advisory-only and cannot veto; report it as heuristic.

**Disposition:** implement the owner clarification under the existing IDs on this feature branch. Do not edit or merge the canonical register. This is a recorded clarification, not a silent redefinition. The register conflict remains visible for review.

## Pre-existing CRIT-154 conflict

`constitution.md` CRIT-154 requires a veto for steric clash `<1.2 A`. The engine resolves LAW-130 as deterministic in the modality matrix, so the CRIT-154 veto remains consistent for that law. The task does not change CRIT-154 or the LAW-130 method class; it only changes the pLDDT measurement scope. No heuristic is used to satisfy CRIT-154.

## Scope exclusions

- No merge, deploy, v25 calibration, held-out evaluation, threshold/cohort change, or golden regeneration (OD-16).
- No OD-10/11/12/13/19/14 implementation.
- No CRIT-154/OBS-153 change.
- No edit to the canonical decision register.
