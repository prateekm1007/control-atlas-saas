# Reconciliation 2026-10-21 - provenance of the archived-branch confirmations

Canonical record for the Round-26 archive-tag provenance reclassifications. It is the
`review_ref` target of the `push-not-proven` records that follow, so it names every reclassified
entry explicitly: `P-011-CONFIRM`, `P-012-CONFIRM`, `P-013-CONFIRM`, `P-014-CONFIRM`.

## Archived-branch provenance reconciliation

The reclassified entries are `P-011-CONFIRM` (branch `proposal/ec16-ec26-enforcement`), `P-012-CONFIRM` (branch `docs/owner-decisions`), `P-013-CONFIRM` (branch `proposal/ec16-ec26-enforcement`) and `P-014-CONFIRM` (branch `governance/canonical-docs`); each is reclassified `push-not-proven` with a `review_ref` to this section.


Round-26 D-V4 pruned the repository's branches after archiving each one as an ANNOTATED tag
`archive/<branch>` at its final tip. The push-log checker resolves those tags as remote evidence
for a pruned branch. Two distinct situations arise, and they are dispositioned differently:

1. **The recorded tip IS the archived tip, or an ancestor of it.** The confirmation is preserved:
   the tag's peeled commit is the remote proof of the recorded commit. This covers the Round-24
   restore set and the R24 candidate (for example the confirmation of
   `repair/r23-governance-recovery` at `781175d18f94b65ee09635f2e4f38bf67d545b45`).
2. **The recorded tip is in NEITHER the archived tip nor its history.** The branch was later
   REBUILT: its early tip is not contained in the final history, so neither the live ref nor the
   archive tag can prove that the recorded push ever happened. Commit existence may still be
   provable; push provenance is not. Those confirmations are reclassified `push-not-proven` - they
   stay visible in the log, they are disclosed as NOT PROVEN on every verification run, and they
   are never counted as evidence for any declaration. The reclassified entries are `P-011-CONFIRM`, `P-012-CONFIRM`, `P-013-CONFIRM`, `P-014-CONFIRM`.

What was searched for each reclassified entry: the live remote ref `refs/heads/<branch>` (absent -
the branch was pruned), the annotated archive tag `refs/tags/archive/<branch>` (present, but its
peeled commit is not the recorded SHA and does not contain it in its history), the local object
store of the repository history (the recorded commit is not reachable from the archived tip), and
the GitHub compare API between the recorded SHA and the archived tip (the recorded SHA is not an
ancestor of the archived tip). No evidence-back record exists for these entries, because the
commits they name are on a lineage that the branch's rebuild discarded.

This record is written by the coder in the same recovery series; it is the canonical location
the `review_ref` mechanism resolves to, and it is disclosed as such in the Round-26 reading note.
