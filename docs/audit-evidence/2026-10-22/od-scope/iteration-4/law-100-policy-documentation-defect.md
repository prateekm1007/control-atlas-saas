# LAW-100 policy-documentation defect record

Recorded: 2026-10-22
Branch: `review/od-scope-v25`
Engine candidate: `3029e21196dedd7cecd3927598639c9bd1a22e6c`
Evidence/archive candidate: `e71cec36d8189125b74568343a80fff1eb9637cd`

## Finding

The owner ruling keeps LAW-100 on the canonical percentage metric for
`(|z| > 4 sigma)` bonds and the `<= 5.0%` decision rule in every source
regime. The current modality matrix still reclassifies LAW-100 as
`advisory_experimental`, and its rationale still reads:

> Bond RMSZ depends on resolution and refinement protocol. False vetoes
> expected across all experimental methods.

That rationale refers to the previous RMSZ metric, which is no longer the
runtime basis for the LAW-100 comparison. The method classification and
hash are preserved for this patch; only the stale rationale is recorded
here as a policy-documentation defect for a separately reviewed update.

## Preserved values

- Modality matrix SHA-256: `ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4`
- Method classification preserved: `advisory_experimental`
- No edit was made to `backend/tos/governance/modality_matrix.json`
- No silent revision of the matrix or its hash is implied

## Required follow-up

Any future correction of the stale rationale must be submitted as a
separately reviewed policy-documentation change, re-pinned in the matrix
and registry, and verified on an exact candidate before activation.
