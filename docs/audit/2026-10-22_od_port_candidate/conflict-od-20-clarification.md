# Conflict Disclosure: OD-20 Clarification File vs Canonical Register Amendment

Round-26 directive 5. Written on candidate
`review/od-scope-v25-2026-10-22` @ `f962c6a`.

## The conflict

1. On the frozen delivery branch `engine/od-scope-v25` @ `97d62a8` there is a
   clarification note `evidence/engine/owner-clarification-od-03-04-20.md`
   (dated 2026-10-21, base `c78ef15`). It restates OD-03/OD-04/OD-20 guidance
   (core `pLDDT >= 70`; LAW-155/200 and LAW-182 remain heuristic and cannot
   independently veto; no new decision IDs).
2. Separately, the branch's commit history claims a "LAW-182 advisory-only"
   registry edit, while the branch's own registry diff
   (`governance/law_registry.yaml` @ `360373d`) contains NO LAW-182 row changes —
   only LAW-105/LAW-130/LAW-155/LAW-200 hunks exist.
3. The canonical register (`main:docs/decisions/OWNER_DECISIONS.md`, blob
   `12906bf5…`) records OD-20's calibration separation (NMR structures
   separated from resolution-based strata). No amendment to the canonical
   register exists for LAW-182.

## Resolution on the candidate (per directive 5)

- The clarification file is treated as a CLARIFICATION (interpretation/scope
  note), not a register amendment. It cannot, by itself, establish a
  LAW-182 "advisory-only" registry row.
- The canonical OWNER_DECISIONS OD-20 content remains the authority for NMR
  separation and modality tiering; the ported/governance code aspects
  (core scoping, admissibility of heuristic failures) come from the reviewed
  registry hunks and code, which the candidate ports.
- Therefore: **no LAW-182 registry row edit on the candidate**; the candidate
  achieves the same OUTCOME through code-side status text only:
  `tier1_measurements.py` reports `FAIL (Advisory)` for heuristic LAW-155/LAW-182/
  LAW-200 and the adjudicator treats heuristic failures as non-VETO contributions.
- The original clarification file is preserved verbatim on the frozen branch
  `97d62a8` (not modified, not deleted, not ported mechanically into
  governance docs).

## Round-26 regeneration note (delivered tests re-certify the outcome)

The reviewed commit `360373d` rewrites `backend/tests/test_law_182.py`
(+35 lines) to pin the advisory OUTCOME (tier-1 status `FAIL (Advisory)`,
adjudicator `method == heuristic` + non-vetoing row contribution
`assert result.verdict == "PASS"`) as test expectations, while the tier-1
registry hunk itself contains no LAW-182 row. This record stays as written
above: the registry gap is a delivery-branch bookkeeping issue, and the
candidate implements the certified outcome code-side only, without
inventing a register row. The clarification file's OD-20 content remains
a clarification, not an amendment; only the owner can amend the canon.

## Branch commit-message overclaim (disclosed)

`engine/od-scope-v25` commit history claims registry registration of
"LAW-182 advisory-only", which the branch diff does not contain. This is
recorded here so the future audit pass can distinguish:

- intended-but-missing registry content (LAW-182), and
- the clarification-note artifact which should NOT be promoted to a registry row
  without an owner register amendment.
