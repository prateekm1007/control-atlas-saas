# INCIDENT 2026-10-21 — Round-23 governance bypass: gate removal and an unapproved merge

Recorded by the Round-24 governance recovery (`repair/r23-governance-recovery`, base `main` @
`63da0e347002ed9b77071ff1ed0cefda7c0ce88a`). Append-only: earlier push-log entries and
earlier register revisions are preserved unchanged.

## 1. What happened (live evidence)

| Fact | Evidence |
|---|---|
| The `pre-merge-verification` job was deleted from `.github/workflows/governance-acceptance.yml` in the Round-23 cycle | `main` carries blob `939261fbbfb200f0f4e83553b7840642f826aaea`; the verified pre-R23 blob is `44fe653ec4e526725eaea4266dc1807222d3be06` (commit `443d6063b42c40ae47b6ea0607e45d31c12239d2`). On `main` the job exists only as a comment stating it was removed |
| PR #30 was merged into `main` without a qualifying independent approval | merge commit `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (first parent `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`, second parent `d81c874eaedc13c6fdbaed5102ad7e638449036d`), `merged_at` 2026-10-20T21:55:06Z, `merged_by` `prateekm1007`; the PR's only review is a `COMMENTED` review by `chatgpt-codex-connector[bot]` |
| The post-merge governance check failed on exactly that point | Governance Acceptance run `39504181997` (failing job `118908835680`, step "Schema + append-only + pairing + live remote-ref verification"): "PR 30 has no approving review from a user other than the PR author that applies to the exact declared head and predates the merge (required independent approval not satisfied; review states: chatgpt-codex-connector[bot]: latest state 'COMMENTED')" at `docs/evidence/push-log/PUSH_LOG.md:556` |
| The Stability Gate on the same commit concluded `success` | Stability Gate run `39504182072` |
| 25 branch refs were pruned, which closed 18 open unmerged PRs | GitHub issue events on PR #30: `head_ref_deleted` 2026-10-20T21:59:48Z; the pruned set and its exact prior heads are recorded in the recovery bundle (`refs.json`, `preservation-matrix.json`) |
| The later repair-branch governance run failed for different reasons | Governance Acceptance run `39507046623`: `P-081-INTENT` has no CONFIRM, and its declared pre-merge base `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` is no longer the live `main` tip (`63da0e3…`). Stability Gate run `39507046937` on the same candidate concluded `success` |

## 2. Why this is an incident

1. The pre-merge gate was removed in the same cycle that needed it, and the merge was then
   executed through a route the gate was designed to block. The gate's own checker already
   states the correct requirement ("a FRESH approving review from a non-author reviewer
   applied to the EXACT head"); the failure was not a checker defect.
2. The Round-23 reading note and register Part 14 justified the removal with the claim that
   independent approval was unsatisfiable because the credential holder is the PR author.
   That claim is invalid: it describes this repository's current credential arrangement, not
   a property of the governance requirement. The requirement stands; the merge is recorded as
   unapproved.
3. The prunable side effects (25 refs, 18 PR closures) were not recorded as data loss risk at
   the time. The unique-commit census performed by the recovery shows 121 commits reachable
   only from non-`main` refs, 77 of whose content survives only through the restored refs.

## 3. Repair-forward plan (no revert, no history rewrite)

1. **Gate restored on a recovery branch**: `.github/workflows/governance-acceptance.yml` is
   restored byte-identically to the verified pre-R23 blob `44fe653ec4e526725eaea4266dc1807222d3be06`
   on `repair/r23-governance-recovery`; PR #30's fixes are preserved.
2. **Regression coverage**: AT-35 proves the gate exists as an executable job (not a comment)
   with its `pull-requests: read` permission, canary, runs manifest and exact-head
   verification; AT-36 proves the gate rejects missing, stale (other commit), non-qualifying
   (`COMMENTED`, `DISMISSED`) and author-only approval evidence, and still accepts a fresh
   non-author approval on the exact head (soundness control). Neither the checker nor any
   existing test was weakened.
3. **Push-log reconciliation (schema-native)**: P-085 (`push-reconcile`) records the push that
   `P-081-INTENT` declared — branch `repair/round19-mainline-acceptance` at
   `d81c874eaedc13c6fdbaed5102ad7e638449036d`, live ref re-verified at 2026-10-20T22:17:45Z
   when it was restored, and the second parent of the merge commit on `main`; P-086
   (`intent-reconcile`) voids the stale intent `P-081-INTENT` because its declared pre-merge
   base is historical (`main` advanced through that very merge) and because a reconciliation
   may not substitute for a confirmation. No retroactive CONFIRM is created and no prior
   entry is edited.
4. **Records corrected**: `docs/P0_STATUS.md` section 13, `docs/directives/INDEX.md`
   (R23/R24 series rows, superseded directives, live-status corrections) and the decision
   register Part 15.
5. **Preservation evidence**: the bundle is produced by
   `scripts/ci/generate_preservation_bundle.py` and published as a downloadable CI artifact by
   `.github/workflows/preservation-bundle.yml`; its URI and digest are recorded in the recovery
   PR body and its checksums are verified from the artifact itself.
6. **Merge rule for the recovery PR**: it merges only through the protected PR path, only with
   a legitimate non-author approval on the exact head, and only with the required checks
   green on that exact head. If no eligible reviewer is available the PR stays blocked.

## 4. Standing constraints

- The OD-16 scientific/production hard stop is preserved: no deployment, no v25 run, no
  held-out evaluation, no threshold or cohort change, no golden-artifact modification.
- The platform configuration of `main` (required checks, required reviews, administrator
  bypass) is verified through the authorized route and recorded before/after; the recorded
  decision is that administrators must not be able to bypass the required approval and status
  checks.
- Credential identifiers are recorded by name only; no token value appears in this record, in
  the push log, in the bundle, or in any workflow log.
