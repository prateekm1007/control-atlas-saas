# META amendment — Directive lifecycle (PROPOSED — PENDING OWNER APPROVAL)

- **Status:** PROPOSED. Not in force. The two lists below are the coder's DRAFT for owner review; the auditor's verbatim text governs in case of conflict.
- **Source:** auditor message, 2026-10-20, section "Two patterns from my own directions need correcting too" (verbatim below).

## Auditor text (verbatim)

> Two patterns from my own directions need correcting too. I issued new directive series while earlier ones were still open, which made the coder unsure which applied. And I wrote "wait for my response" into several items, which stalled work that could have proceeded. The constitution should say which items may proceed without a response, and which need one.

## Proposed constitutional text

**R-DL1 One active series.** Exactly one directive series is ACTIVE at a time. A new series may be opened only after the previous one is COMPLETE, or by explicitly marking it SUPERSEDED in `docs/directives/INDEX.md` with a pointer to its replacement. Work performed against a non-active series is void and must be re-based on the active series.

**R-DL2 Series identity in every report.** Every report, PR, and evidence folder names its directive series ID (EC-19) and lists the items it closes.

**R-DL3 Proceed-without-response.** The following item classes MAY proceed without a further auditor/owner response once the directive issuing them is active. A directive may narrow this list; it may not widen it.

1. Pushes to non-`main` branches (per the standing push authorization, OD-15).
2. Append-only evidence updates: evidence folders, the push log (EC-24), receipts, sequence logs (EC-18).
3. Drafting proposals and their inert scaffolding (this package is the first instance).
4. Documentation of observed state: verification runs, GET-only remote checks, lint outputs, status recordings.
5. Bundle/manifest creation (`SHA256SUMS`, format-patches) for completed items.

**R-DL4 Response-required.** The following item classes require an explicit recorded response before execution, regardless of any standing authorization:

1. Anything that merges to `main`, or deploys, or changes a verdict outside a recorded diff.
2. Threshold, protocol, cohort-definition, or registry changes (C3–C11 remain closed).
3. Opening a new directive series while one is ACTIVE (R-DL1).
4. Credential operations: adding, rotating, or revoking credentials; anything that spends money.
5. Any deviation from a ruling's text without a written waiver (EC-20).

**R-DL5 Default when unspecified.** If a directive is silent about whether an item waits, the item falls on the R-DL3 list by default. A "wait for my response" clause is valid only if it cites which R-DL4 class the item belongs to; otherwise the wait clause is void and the item proceeds.

## Coder notes

- R-DL5 is the direct answer to the auditor's observed failure ("I wrote 'wait for my response' into several items, which stalled work that could have proceeded") and is the coder's draft for owner review.
- The ACTIVE series at the time of writing: the EC-16..EC-26 proposal packaging (this branch). The prior six-direction series (Round-12) is COMPLETE per worklog Task 35; the post-Task-35 remediation commits (c02f774 / 5149435 / 0773f0b) closed the Round-13/14 findings — see `docs/directives/INDEX.md` (rows RECONSTRUCTED).
