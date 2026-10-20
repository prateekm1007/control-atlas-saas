# Push log schema v4 — two-phase records + post-hoc reconciliation + remote verification + binding/merge

(Defined 2026-10-20, EC-18/EC-24 auditor revisions; v3 added the `push-reconcile`
verb and the remote-verification layer per the Round-20 auditor direction 4; v4
adds the Round-21 identifier-pairing, terminal-intent exact-candidate binding,
historical-proof (`evidence-back`, `intent-reconcile`) and protected-branch
merge-transition (`merge-intent`) rules per auditor directions 4–6.)

This document defines the record format for `PUSH_LOG.md`. It is a separate file so the
log itself stays append-only under the parent-prefix check (`scripts/ci/evidence_sequence_check.py
--parent-prefix`): defining the schema here never requires editing prior log entries.

## Defect being corrected (auditor, item 3.6)

The P-011 entry recorded `sha=self-referential-tip` — not a SHA — because an entry shipped
inside the commit it describes cannot contain that commit's hash. It also stated a push time
that could not be known until the push returned. Both are fixed by splitting each push into
two truthful, separately-timestamped records:

## Record phases

1. **INTENT (pre-push).** Written and committed BEFORE the push is issued. May contain only
   facts knowable at that moment:

       ## <UTC> | P-<n>-INTENT | push-intent | branch=<name> head_at_intent=<40hex|DEFERRED> expected_remote_tip_before_push=<40hex|ABSENT> credential=<name> authorized_by=<reference>

   `expected_remote_tip_before_push` is the verifiable pre-push state: the exact 40-hex tip
   the remote is expected to hold (v4: it must equal the newest CONFIRM/RECONCILE recorded
   for the branch before this intent, and be an ancestor of the pushed candidate), or
   `ABSENT` for the first push of the ref name (no prior branch record may contradict it).

2. **CONFIRM (post-push).** Appended in a LATER commit, after the push returned, with its own
   truthful `date -u` timestamp and the remote evidence source:

       ## <UTC> | P-<n>-CONFIRM | push-confirm | branch=<name> remote_sha=<40hex> credential=<name> evidence="git ls-remote origin refs/heads/<name>"

   `remote_sha` MUST be the full 40-hex SHA observed on the remote (placeholder markers such as
   `self-referential-tip` are rejected by `--pushlog-schema`). CONFIRM entries for a branch's
   own tip are published by a SUBSEQUENT push of that branch; the very last confirm of a
   session is recorded in the round report/worklog and appended to this log by the next
   legitimate push (disclosed limitation, never backdated).

   **v4 identifier pairing (direction 4):** `P-058-CONFIRM` consumes exactly the unpaired
   `P-058-INTENT` (or legacy placeholder push `P-011`) for the SAME branch. A confirmation
   whose identifier has no matching intent FAILS; it can never consume another intent.

3. **RECONCILE (post-hoc, v3).** For pushes that already happened without a preceding INTENT
   — the Round-19 follow-up-push class found by the Round-20 audit — a truthful, dated record
   that reconciles the push without pretending an INTENT existed:

       ## <UTC> | P-<n> | push-reconcile | branch=<name> remote_sha=<40hex> credential=<name> evidence="git ls-remote <url> refs/heads/<name> -> <sha> (re-verified <UTC>); <run citations>"

   A RECONCILE never pairs with or consumes INTENT records (it is post-hoc by design), never
   backdates, and its evidence must contain both the recorded SHA and the branch name, plus
   the live-ref verification performed at reconciliation time and (when available) check-run
   citations for the exact head.

4. **EVIDENCE-BACK (v4, direction 5).** A reviewed record proving that a HISTORICAL
   confirmation whose SHA left the branch lineage describes a real pushed commit:

       ## <UTC> | P-<n> | evidence-back | backs=<confirmed-entry-id> branch=<name> remote_sha=<40hex> credential=<name> evidence="<machine-verifiable facts: commits-API existence, run citations with head_branch, preservation refs, the reviewed rebuild direction>"

   The checker requires the backed SHA to be a real remote commit (local object store or the
   commits API) and verifies every run the record cites (exact head_sha; head_branch connects
   the SHA to the claimed branch). A textual transcript alone never satisfies the proof.

5. **INTENT-RECONCILE (v4, direction 4).** Voids a named stale/superseded intent so it does
   not count as an excess unmatched intent for its branch:

       ## <UTC> | P-<n> | intent-reconcile | voids=<intent-entry-id> branch=<name> credential=<name> evidence="<why the intent is void>"

   It never converts an intent into a confirmation and never conceals missing evidence: the
   voided intent stays visible in the log; only its pairing obligation is discharged.

6. **MERGE-INTENT (v4, direction 6).** Pre-declares a future merge of this branch into a
   protected branch through the reviewed PR mechanism:

       ## <UTC> | P-<n> | merge-intent | target=<protected-branch> pr=<number> source_branch=<name> expected_pre_merge_base=<40hex> expected_head=<40hex|DEFERRED> credential=<name> authorized_by=<reference>

   `expected_head=DEFERRED` means the head is the commit that introduces this record (a
   record cannot embed its own commit's SHA); the checker resolves it via the introducing
   commit. On the source branch, the active merge-intent must bind the exact candidate and
   the declared base must still be the live tip of the target. After the merge, a push event
   on the protected branch is DECLARED only by a VERIFIED merge-intent: the candidate must
   be a merge commit whose first parent is the declared pre-merge base and whose second
   parent is the declared head, the merge-intent must ride the merged PR branch history,
   and the PR must be merged with an approving review from a user other than the PR author.
   Raced, unrelated, or unapproved updates FAIL. The merge must be executed as a MERGE
   COMMIT through the protected PR path — never a direct push, never a bypass.

## Remote verification (v3/v4, `--remote-verify`)

`scripts/ci/evidence_sequence_check.py --remote-verify [--branch <under-test>
[--expect-tip <sha>]] [--runs-manifest FILE] [--repo <slug>]` verifies, fail-closed:

- every CONFIRM/RECONCILE's recorded `remote_sha` against the **live remote ref** for its
  claimed branch (`git ls-remote`): equal => the confirmation is CURRENT; different => it is
  HISTORICAL — valid for the push it describes, and explicitly NOT a confirmation of the
  newer tip;
- that the `evidence=` field actually **supports** the recorded SHA: the SHA and the claimed
  branch must appear in the evidence (or the entry body), and any cited Actions run ID whose
  manifest record shows a different head_sha FAILS (fabricated support);
- **(v4) that every run cited anywhere in an entry (evidence field or body) verifies**:
  found in the runs manifest, or fetched through the Actions API (`--repo` + token) when
  outside the manifest window; an unverifiable citation FAILS the gate — a manifest omission
  is never an approval. Failed or cancelled runs are recorded honestly and cannot prove a
  successful push;
- **(v4) that a HISTORICAL confirmation is proven real**: either the recorded SHA is an
  ancestor of the branch's live tip (it is on the current lineage), or — when the SHA left
  the lineage — a reviewed `evidence-back` record exists whose SHA is a real remote commit
  (object store or commits API) and whose run citations verify. A different SHA alone, or a
  transcript naming a branch and SHA, is not proof;
- that the branch under test's CURRENT tip is DECLARED: by a bound unpaired terminal INTENT
  (reported PENDING — the honest steady state for the newest push of a session; v4 requires
  the terminal intent to have been introduced by the exact candidate under test and to
  declare the verified pre-push tip, with all excess unmatched intents voided), by an
  exact-tip CONFIRM/RECONCILE, or — on a protected branch — by a VERIFIED merge-intent. A
  tip whose push was never declared FAILS (an ls-remote failure or an absent remote ref
  also fails closed — missing remote verification is never a pass).

## Truthfulness rules

- INTENT and CONFIRM are separate entries with separate IDs; an INTENT entry never claims
  post-push knowledge, a CONFIRM entry never rewrites its INTENT.
- RECONCILE records are visibly post-hoc: they disclose the missing INTENT rather than
  manufacturing one, and they never backdate a push time. EVIDENCE-BACK and
  INTENT-RECONCILE records are likewise append-only additions; no prior entry is edited.
- Timestamps come from `date -u` at the moment the record is written. Historical entries
  whose timestamps could not be observed at action time are labeled `RECONSTRUCTED`.
- Prior entries are never edited or deleted (enforced by the parent-prefix check against
  the previous commit's version of the log).
- A CONFIRM entry reconciles (in its body) any difference between the INTENT/entry timestamp,
  the git commit metadata, and the remote observation; differences are explained, not averaged
  away, and unprovable claims are marked as such.
- Credential identifiers are recorded by NAME only; token values are never written anywhere.
