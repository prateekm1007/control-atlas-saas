# INCIDENT-2026-10-21-R27 — Round-23 governance bypass: gate removal and an unapproved merge

> **Ported record (Round-27, auditor decision D-W2).** This incident note was
> written on the `repair/r23-governance-recovery` line (`781175d`) and is
> ported here so the record survives the closure of that branch (D-W3). The
> identifier is NEW: `INCIDENT-2026-10-21-R27-GOVERNANCE-BYPASS` is not used
> anywhere else in the tree. The incident FACTS are preserved verbatim;
> section 3 is adjusted to this line and section 5 records exactly what was
> and was not imported. The original identifier
> `INCIDENT-2026-10-21-R23-GOVERNANCE-BYPASS` remains on the archived
> branch and on its tag, and is not used in this tree.

Recorded by the Round-24 governance recovery (`repair/r23-governance-recovery`, base `main` @
`63da0e347002ed9b77071ff1ed0cefda7c0ce88a`); ported to the single recovery line
(`repair/round24-claim-head-binding`) in Round-27. Append-only: earlier push-log
entries and earlier register revisions are preserved unchanged.

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

The plan below is the one carried out on the `781175d` line; the equivalent
state on THIS line (the single recovery line) is:

1. **Gate restored**: `.github/workflows/governance-acceptance.yml` on this line
   is the verified pre-R23 blob `44fe653ec4e526725eaea4266dc1807222d3be06`
   with ONLY the intended Round-24/Round-25 additions re-applied (the
   fail-closed live PR-body retrieval, the explicit `--expect-head-sha`, the run
   identity in the runs manifests, the acceptance-suite label, the
   preservation-bundle jobs and the restoration comment). The exact diffs are
   recorded in `evidence/r25/gate-diff.txt` and
   `evidence/r25/gate-diff-after-restore.txt`; the file is blob
   `358dea6e4828f8ff77cce4dfdfcf0705d7dba13d`. PR #30's fixes are preserved:
   the merge is NOT reverted and no useful code is removed.
2. **Regression coverage**: AT-39 proves the gate exists as an executable job
   (not a comment) with its `pull-requests: read` permission, its canary, its
   runs manifest and its exact-head verification, and that the Round-23 removal
   claim is absent; AT-40 proves the gate rejects missing, stale (other
   commit), non-qualifying (`COMMENTED`, `DISMISSED`) and author-only approval
   evidence, and still accepts a fresh non-author approval on the exact head.
   (The `781175d` side carried these as AT-35/AT-36; they are renumbered here
   because this line already registers AT-35..AT-38 for different controls. No
   identifier is reused and no existing test was weakened.) AT-41 proves that a
   pruned branch keeps its content in an annotated archive tag and that the
   push-log checker resolves that tag as remote evidence.
3. **Push-log reconciliation (schema-native)**: the `d81c874` push is recorded
   in THIS line's log as `P-081-CONFIRM` (the second parent of the merge
   commit on `main`). No retroactive CONFIRM is created and no prior entry is
   edited.
4. **Records corrected**: `docs/P0_STATUS.md` (sections 13-16),
   `docs/directives/INDEX.md`, the decision register (Parts 16, 16.9 and 17),
   and `docs/evidence/push-log/INCIDENT-2026-10-21-R25-PR-IDENTITY.md` for the
   pull-request identity facts and the push-log identifier analysis.
5. **Preservation evidence**: the bundle is produced by
   `scripts/ci/generate_preservation_bundle.py` (byte-identical to the
   `781175d` version) and published as a downloadable CI artifact; its digest
   is recorded and every member's SHA-256 is re-verified from a fresh
   checkout. The prescribed port verification is in `evidence/r27/port-verify.txt`.
6. **Merge rule for the recovery PR**: it merges only through the protected PR
   path, only with a legitimate non-author approval on the exact head, and only
   with the required checks green on that exact head. If no eligible reviewer is
   available the PR stays blocked.

## 4. Standing constraints

- The OD-16 scientific/production hard stop is preserved: no deployment, no v25 run, no
  held-out evaluation, no threshold or cohort change, no golden-artifact modification.
- The platform configuration of `main` (required checks, required reviews, administrator
  bypass) is verified through the authorized route and recorded before/after; the recorded
  decision is that administrators must not be able to bypass the required approval and status
  checks.
- Credential identifiers are recorded by name only; no token value appears in this record, in
  the push log, in the bundle, or in any workflow log.

## 5. What was imported from the 781175d line, and what was not (auditor D-W2)

**DEVIATION, DISCLOSED.** D-W2 says the `781175d` push-log entries are NOT
imported and that their facts are recorded in the incident note instead. The
published Round-26 fold on this branch (`3d0b5ab`) DID import them,
RENUMBERED: the `781175d` line's `P-082-INTENT`/`P-082-CONFIRM` (the push of
`repair/r23-governance-recovery` at `781175d`) are present here as
`P-090-INTENT` and `P-090-CONFIRM`, and the record bodies say so explicitly
("this record is folded into PR #37 in this round and is renumbered from
P-082-INTENT, which this head had already used"). No earlier entry was edited
and no duplicate ENTRY-ID exists, so the append-only rule and the schema check
hold; but the deviation from D-W2's letter is recorded here rather than hidden.
The import is a record of a real push of a branch that existed, with the run
evidence of that push (Stability Gate run 39532002306, Governance Acceptance
run 39532002077, both `success`).

What was ported from `781175d` into this line:

- `scripts/ci/generate_preservation_bundle.py` - byte-identical (the prescribed
  `git diff 781175d <head> -- scripts/ci/generate_preservation_bundle.py` is
  empty; `evidence/r27/port-verify.txt`);
- the two acceptance tests that came with it, renumbered AT-35/AT-36 ->
  AT-39/AT-40 (this line already registers AT-35..AT-38);
- this incident note, under a new identifier.

What was NOT imported: the `781175d` push-log entries under their ORIGINAL
identifiers (they would duplicate this line's P-082..P-086 space), and the
`781175d` branch's own `.github/workflows/preservation-bundle.yml` workflow
file - the bundle publication on this line is the live `preservation-bundle`
job of `governance-acceptance.yml`, which runs the same (ported) generator and
adds a fresh-checkout digest verification. The `781175d` head is preserved as
the tag `archive/repair-r23-governance-recovery`, so nothing it contains is
lost.
