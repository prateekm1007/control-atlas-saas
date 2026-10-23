# OD-scope v25 reconstruction — required reading note

- Date: 2026-10-23
- Branch: `reconstruct/od-scope-v25`
- Base commit: `e71cec36d8189125b74568343a80fff1eb9637cd` (`review/od-scope-v25`)
- Candidate: under construction; the terminal commit SHA will be recorded in the exact-candidate verification report and in P-091's metadata after the final commit exists.
- Governing directive: owner directions of 2026-10-23, including rejection of `01b3d3f0508d00e2b1d638a32d82b9cd7dd1cbe6` as the push candidate, reconstruction from the live base, independent engine/calibration audit, and exact-candidate verification before any push.

## Complete-read record

The following files were read in full in this session before any edit on this branch. Hashes are from the base commit, not the working tree.

| File | Bytes | SHA-256 | Git blob | Complete-read timestamp (UTC) |
|---|---:|---|---|---|
| `AGENTS.md` | 4021 | `5e320939a11135a0d8aed5caa231df6aeaf633d724fcb4dba97b5909e51420b4` | `85767237399ffb30f3e45894d9dfbe6f4240584f` | 2026-10-23T01:57:32Z |
| `constitution.md` | 622 | `cfdd5ba3babda31513b55894a49597fbfb18c2e870535a068c8b0ea423c53ce8` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | 2026-10-23T01:57:32Z |
| `ENGINEERING_CONSTITUTION.md` | 20572 | `e0f256f78f2c1256228c4286311cba15800dcd5a6c8005583508ee66d27983f5` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | 2026-10-23T01:57:32Z |
| `docs/coder/BRIEF.md` | 3671 | `8f07edfbba9d3e57de623e242125e942369bc5d60df5805343b26b9048c52e07` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` | 2026-10-23T01:57:32Z |
| `docs/decisions/OWNER_DECISIONS.md` | 72329 | `59ccd5e84176a8a610c8c234cbb95fec4734478cc1f35668599ad2e3aef3172a` | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` | 2026-10-23T01:57:32Z |
| `docs/P0_STATUS.md` | 35347 | `3f129dd18b9e5cf8888b763158cda42de270a8429a48ed7cc9a3aa675e810e0f` | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` | 2026-10-23T01:57:32Z |
| `docs/evidence/push-log/PUSH_LOG_SCHEMA.md` | 14971 | `be29828965ebad9e9a300a12890777f0debcb1eaab2de1e551e64414d78acb42` | `83aa798b7a9731e63be2a9664e502b7631a42674` | 2026-10-23T01:57:32Z |
| `scripts/ci/evidence_sequence_check.py` | 115064 | `39370386d69e9062e9e3b4db7ecfd1e980c9cfc81ad90fe95a4e60a3abb2ae3c` | `1ae66acf131b59b9b05f8cd06ad9262e8eab74bb` | 2026-10-23T01:57:32Z |

All eight reads were completed before the first edit on this branch. The timestamp records completion of the reading batch; it is not a claim that each file was opened at a distinct second.

## Policy clauses that bind this work

- `OD-16` remains in force: no merge, deployment, v25 calibration execution, held-out evaluation, threshold/cohort change, or golden-artifact regeneration.
- The append-only push log is evidence, not editable narrative. Prior entries must be preserved byte-for-byte; new state is appended only at the end.
- A terminal push intent must be introduced by the exact candidate commit, must declare the verified pre-push remote tip, and remains `PENDING` until the post-push confirmation exists.
- The engine must fail closed on malformed, missing, non-finite, or degenerate measurements; invalid geometry is not a favorable zero.
- pLDDT must be observation-backed. B-factors must never be relabeled as pLDDT.
- Tests must verify behavior, not merely status labels or happy paths.
- Existing law-semantic/canonical/matrix hash locks must remain unchanged unless the owner authorizes a scientific change.
- No force-push, no history rewrite after publication, and no unrelated commits after the terminal intent.
- The calibration protocol must be source-backed and provenance-audited; no protocol may be invented. OD-16 still prevents executing it.
- A local pass is not remote CI. The exact pushed SHA and every required push-event job must be verified after the push.
- The owner PAT is rotated by the owner; it must not be exposed, committed, logged, or persisted.

## Baseline and rejected candidate

The live remote tip of `review/od-scope-v25` is `e71cec36d8189125b74568343a80fff1eb9637cd`. The unpublished local chain is:

`e71cec3` → `459dd55` → `fc11d10` → `f0228a7` → `c44457e` → `01b3d3f`

Independent byte-prefix checks found two append-only violations:

1. `459dd55` rewrites a pre-existing whitespace sequence in PUSH_LOG.md (`candidate. If` becomes `candidate.   If`). No log entry is deleted, but the committed byte stream is not a prefix of its parent.
2. `01b3d3f` deletes the P-087/P-088 bytes and replaces them with P-090/P-091. This is destructive and disqualifies it as the push candidate.

`P-089` is unused. `P-090` and `P-091` are unused on the reconstructed lineage. The reconstruction plan is to restore the whitespace byte, replay `fc11d10`, `f0228a7`, and `c44457e`, then append a new `P-090` that voids `P-088-INTENT` (not `P-086`) and a new `P-091-INTENT` that is introduced only by the terminal commit and declares `expected_remote_tip_before_push=e71cec3`.

## Scope and non-goals

In scope: preserve the rejected local work in recovery refs/bundles; reconstruct the unpublished branch from `e71cec3` without the destructive log edit; retain every legitimate engine, test, calibration, and evidence change; audit those changes on the reconstructed SHA; regenerate the final manifest after the push log is final; validate the complete history; fast-forward push only if every gate passes.

Out of scope: merge, deployment, v25 calibration execution, held-out evaluation, changing thresholds/cohorts/verdicts, rewriting published remote history, force-push, or treating a local test pass as remote CI evidence.

## Acceptance criteria

- Every reconstructed commit's push-log bytes are a prefix of its parent (or a no-op for that file).
- The schema, pairing, chronology, global unmatched-intent, remote-reference, citation, and terminal-intent checks pass on the exact final candidate using the same flags as CI.
- `P-091-INTENT` is the sole terminal unmatched intent for `review/od-scope-v25`, is introduced by the final commit, and declares `e71cec3`.
- Every row in `SHA256SUMS` matches the final tree; no stale digest remains.
- Secret scan, forbidden-reference check, workflow-trigger lint, and package integrity checks pass.
- Targeted engine/adjudicator/calibration tests and the full suite pass on the exact candidate SHA; results are recorded with commands, exit codes, and counts.
- Engine audit covers LAW-100, LAW-120, LAW-125/135/145/150/160/195 empty/degenerate behavior, geometry validity, LAW-130, adjudicator fail-closed behavior, and pLDDT/B-factor separation.
- Calibration audit records protocol provenance and evaluator/runner fail-closed behavior; readiness is reported as `PASS`, `FAIL`, or `BLOCKED`, not inferred from a green suite.
- The live remote still points at `e71cec3` immediately before the push, and the push is a fast-forward.
- Post-push: live SHA equals the candidate and every required push-event workflow job is verified from the actual run logs.

## Pre-change attestation

No source, test, workflow, push-log, manifest, or governance content has been modified on this branch before this note and the accompanying preflight were written. The only prior actions were recovery preservation in the separate local repository and creation of this worktree at the unmodified base.
