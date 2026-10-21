# INCIDENT-2026-10-21-R25-PR-IDENTITY — the two recovery pull requests and the
# containment check that did not hold

**Recorded:** 2026-10-21T04:00:00Z (Round-25 verification round, auditor direction 3).
**Status:** OPEN. This note records FACTS and the containment difference; it
does not decide anything and it does not close anything.

## 1. The pull-request facts (live API, `2026-10-21T04:00:00Z`)

| PR | author login | head branch | head SHA | base | state | merged | merge_commit_sha | created |
|---|---|---|---|---|---|---|---|---|
| #37 | `prateekm1007` | `repair/round24-claim-head-binding` | `bc6fdb5d92a237d421570d1a462199e4cdcabf28` | `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` | open | no | `6466029d9b68388875d8536855c1e0cdc4940a92` (GitHub's computed merge commit, not yet realised) | 2026-10-20T23:49:20Z |
| #38 | `prateekm1007` | `repair/r23-governance-recovery` | `781175d18f94b65ee09635f2e4f38bf67d545b45` | `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` | open | no | `d282259a8ae41b4e8854ec5fb1cc13365a8f3851` (GitHub's computed merge commit, not yet realised) | 2026-10-21T02:10:07Z |

Both pull requests were opened by `prateekm1007` - the credential holder and
therefore the PR author for both. No other account appears as an author, a
reviewer or a committer on either pull request. This is the fact behind the
auditor's D-U3: a non-author approval is required for PR #37, and this coder
cannot produce one (self-approval is forbidden and the credential holder is the
author).

## 2. The containment check (auditor direction 4) DID NOT HOLD

```
git merge-base --is-ancestor 781175d18f94b65ee09635f2e4f38bf67d545b45 bc6fdb5d92a237d421570d1a462199e4cdcabf28   -> exit status 1 (NOT an ancestor)
git merge-base 781175d18f94b65ee09635f2e4f38bf67d545b45 bc6fdb5d92a237d421570d1a462199e4cdcabf28                  -> 63da0e347002ed9b77071ff1ed0cefda7c0ce88a
```

The two recovery branches are **siblings**: they diverge from the same merge
base (`main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`) and neither is an ancestor of the other. The
prescribed mechanical test therefore fails, and the direction's own condition
("if it is not, stop and report the difference") applies. Nothing was closed,
tagged, deleted or merged on the strength of a containment claim.

Commits only in PR #38's head (not in the R25 head):

- `6dcf851` fix(R24): restore the pre-merge gate, add AT-35/AT-36 coverage,
  publish the preservation bundle, correct the live records
- `781175d` fix(R24): correct the tracked manifest row, make the bundle
  checksum format sha256sum-compatible, reconcile the historical P-081
  push-log state

Commits only in the R25 head (not in PR #38's head):

- `6ff3402` fix(R24): bind the PR-body claim gate to the exact current head
- `2ce77e5` fix(R24-2): pin the package manifest over the committed bytes
- `7ffe471` fix(R25): fail-closed live PR body, required-workflow CI evidence,
  committed-tree manifest, restored pre-merge gate, preservation bundle, A-04
- `bc6fdb5` docs(R25): complete the push-log protocol for the Round-25
  candidate

## 3. The content difference, in both directions

Only in PR #38's head (19 files differ in total; the deletions and additions
below are the material ones):

- `scripts/ci/generate_preservation_bundle.py` - a DIFFERENT bundle generator
  with capabilities the R25 generator does not have: the 25-ref accounting in
  `refs.json`, unique commits with **patch-ids**, a per-branch merge base and
  ahead/behind with renames, a preservation disposition **per unique commit**,
  a runtime dependency map that names the references/workflows depending on the
  restored branches, PR claim extraction with verification outcomes, and a
  **fail-closed** exit when a required API read fails (the R25 generator
  records an unavailable read and continues).
- `docs/evidence/push-log/INCIDENT-2026-10-21-R23-GOVERNANCE-BYPASS.md` - the
  R23 governance-bypass incident note (the R25 head records the same facts in
  register Part 16.5 and P0_STATUS section 14, but not as a push-log incident
  note under this identifier).
- a different reconciliation of the historical P-081 push-log state and a
  different tracked manifest row.

Only in the R25 head:

- `scripts/ci/fetch_pr_body.py` - the fail-closed live PR-body retrieval (the
  event-payload fallback removed) that PR #38 does not have.
- the required-workflow identity in the claim gate, the committed-tree
  (HEAD) manifest semantics, the two preservation-bundle CI jobs with their
  fresh-checkout digest verification, and the adversarial tests AT-36, AT-37 and
  AT-38.

**Consequence:** the two push logs are MUTUALLY EXCLUSIVE. The entry headings
have the grammar `## <UTC timestamp> | <ENTRY-ID> | <verb> | <subject>` - the
identifier is the SECOND field, so a search for `## P-08x` finds nothing; the
identifiers are found with `^## \S+ \| (P-\S+) \|`. With that pattern the
P-08x identifiers present on the two sides are:

- on the `781175d` side: `P-080`, `P-081-INTENT`, `P-082-INTENT`,
  `P-082-CONFIRM`, `P-083-INTENT`, `P-084`, `P-085`, `P-086`;
- on this side: `P-080`, `P-081-INTENT`, `P-081-CONFIRM`, `P-082-INTENT`,
  `P-082-CONFIRM`, `P-083-INTENT`, `P-083-CONFIRM`, `P-084-INTENT`,
  `P-084-CONFIRM`, `P-085`, `P-086-INTENT`, `P-086-CONFIRM`, `P-087`,
  `P-088`, `P-089-INTENT`, `P-089-CONFIRM`, `P-090-INTENT`, `P-090-CONFIRM`,
  `P-091`, `P-092`, `P-093-INTENT`, `P-093-CONFIRM`, `P-094-INTENT`, `P-095`;
- the IDENTIFIERS COMMON TO BOTH SIDES: `P-080`, `P-081-INTENT`,
  `P-082-INTENT`, `P-082-CONFIRM`, `P-083-INTENT`, `P-085` - and they denote
  DIFFERENT branches on the two sides (`P-082-INTENT` is the first push of
  `repair/r23-governance-recovery` there and of
  `repair/round24-claim-head-binding` here).

So the two logs cannot be merged: the published log would carry duplicate
ENTRY-IDs and the schema check would fail. Exactly one of the two heads can
become the recovery line, which is consistent with the auditor's D-U1 (one
recovery PR). (The same collision existed in the acceptance-test registry:
AT-35/AT-36 meant the pre-merge-gate tests on the `781175d` side and the
exact-head/live-body tests here.)

**RESOLUTION (Round-27).** The published Round-26 fold (`3d0b5ab`) moved that
content into this branch under RENUMBERED identifiers: the `781175d` entries
`P-082-INTENT`/`P-082-CONFIRM` are present here as `P-090-INTENT`/
`P-090-CONFIRM`, the pre-merge tests are `AT-39`/`AT-40`, and the incident
note is under the new identifier `INCIDENT-2026-10-21-R27-GOVERNANCE-BYPASS`.
The deviation from D-W2's "do not import" is recorded in that note's section
5 rather than hidden. Exactly one of the two heads can become the recovery line, which
is consistent with the auditor's D-U1 (one recovery PR) - but D-U1's
precondition ("its content is confirmed contained in the R25 head") is NOT yet
satisfied, because PR #38's bundle generator and its R23-bypass incident note
exist only on its side.

## 4. What is needed to close PR #38 (D-U1)

Either

1. the auditor rules that the R25 head is the recovery line and accepts that
   PR #38's `generate_preservation_bundle.py` and its R23-bypass incident note
   are superseded rather than ported (the R25 bundle already publishes the nine
   required members, is verified from a fresh checkout, and its checksum
   format is sha256sum-compatible - the CI verification job runs
   `sha256sum -c SHA256SUMS` successfully on the R25 head); or
2. the R25 head absorbs the outstanding #38-only content first (the
   fail-closed bundle behaviour, the patch-id unique-commit records, the
   per-commit preservation disposition and the incident note), after which the
   containment is true by content and the direction's mechanical test can be
   replaced by a content comparison.

Until one of those is decided, PR #38 stays OPEN and its branch is NOT deleted
(a branch deletion without an archive tag is forbidden, and the archive tag
belongs to the closure step that is conditioned on the containment check).
