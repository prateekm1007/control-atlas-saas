# Owner clarification - OD-03 / OD-04 / OD-20

**Date:** 2026-10-21
**Branch:** `engine/od-scope-v25`
**Base:** `c78ef15f696a2b1e2658e7ac7b54461f70c6e040`
**Canonical register reviewed:** `origin/main:docs/decisions/OWNER_DECISIONS.md` at blob `12906bf5997a3e2c6b5c41e5d5c60fee14613274`.

## Direction recorded

The owner directed that the task brief governs these existing decision IDs and that no new decision IDs be created:

- **OD-03:** core pLDDT `>=70` for LAW-100/105/110/120/125/130/135/145/150/155/200; experimental structures retain all atoms; LAW-155/200 remain heuristic and cannot independently veto. LAW-130 remains deterministic per the modality matrix.
- **OD-04:** one cutoff at `70`; remove the old `50` boundary; a structure with `<70%` of atoms at pLDDT `>=70` is `INDETERMINATE`, not `PASS`.
- **OD-20:** LAW-182 is advisory-only, cannot veto, and must be reported as heuristic.

The owner also directed:

- Keep the existing OD-03, OD-04, OD-05, and OD-20 IDs.
- Do not create OD-22 through OD-25.
- Retain LAW-190/210/220/230 as retired/not enforced (OD-05); do not delete their IDs.
- Delete the tracked nested `control-atlas-saas/` duplicate tree.
- Fetch the pinned v24 coordinate URLs and verify their recorded SHA-256 values.
- Cite the `origin/main` governing blobs in the reading note; do not merge.

## Conflict disclosure

The canonical register text differs materially from the task brief for OD-03, OD-04, and OD-20. The owner reviewed that conflict and selected the brief text as controlling for this branch. This file is the clarification record. It does not amend the canonical register and must not be represented as a register edit.

The canonical OD-05 text already matches the owner's retain-retired-IDs direction.
