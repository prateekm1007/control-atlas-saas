# Auditor Decision — LAW-155 Classification Reconciliation (2026-10-24)

**Decision ID:** `AUDITOR-DECISION-2026-10-24-LAW155`
**Authority:** delegated auditor decision under the owner's direction (2026-10-24)
**Status:** ACTIVE (records the superseding classification direction)
**Supersedes:** the earlier retirement of LAW-155 in the decision register (OD-03
retirement language) and the `RETIRED_LAWS` entry for LAW-155 in
`backend/tos/governance/station_sop.py`, as of this date.

## Decision

For the present roadmap, **LAW-155 (Voxel Occupancy) is reclassified from
retired to heuristic/advisory: it is advisory and non-blocking pending
negative-control validation.** It must:

- NOT be marked retired again;
- NOT function as a deterministic veto;
- remain visible and labeled as advisory/heuristic in API and PDF output;
- retain its existing numerical threshold, operator, unit, and scientific
  meaning (no threshold change, no new scientific meaning).

## Supersession

This decision supersedes the OD-03 retirement language and the prior
`RETIRED_LAWS` classification for LAW-155. The historical record is not edited;
this document records the supersession. Hash locks affected by this
classification reconciliation are re-pinned in the implementing patch, with
old and new digests recorded there.

## Scope

The decision authorizes the *classification* reconciliation only. It does not
authorize changing LAW-155's threshold, operator, unit, scope, or meaning, and
it does not authorize any calibration, cohort, held-out, or threshold change.
OD-16 remains in force.

## Digests (old -> new)

- Canonical law hash: `f24ec581b3e04732` -> `6a9cd4b4349b81de`
- Law-semantic hash: `d2b3b8501b1c73f67e7839662f70b68596d3c1e25b51d824d6b8a78325816dd5`
  -> `f825f38b9493a9409aa9678c1f98593e7715ba6b9e01aefba808096571a27750`
- Modality-matrix hash: unchanged (locked value `ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4`).
