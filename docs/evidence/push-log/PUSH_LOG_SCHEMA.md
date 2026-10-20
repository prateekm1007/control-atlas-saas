# Push log schema v5 — two-phase records + post-hoc reconciliation + remote verification + binding/merge + claim separation

(Defined 2026-10-20, EC-18/EC-24 auditor revisions; v3 added the `push-reconcile`
verb and the remote-verification layer per the Round-20 auditor direction 4; v4
added the Round-21 identifier-pairing, terminal-intent exact-candidate binding,
historical-proof (`evidence-back`, `intent-reconcile`) and protected-branch
merge-transition (`merge-intent`) rules per auditor directions 4–6; v5 adds the
Round-22 claim-separation rules per auditor directions 3–5: `push-not-proven`
reclassification, `review_ref=` adjudication references, recorded `push_event=`
citations, PushEvent before/after pre-push-state binding, the global
unmatched-intent invariant, and the hardened merge/pre-merge verification.)

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
   v5 rules: the target must be a PRECEDING eligible intent; the branch scope must
   match (voiding a push-intent requires the reconcile's `branch=` to equal the
   intent's branch; voiding a merge-intent requires it to equal the intent's
   `source_branch` — unambiguous source/target semantics, never a free-form
   branch value); an intent is voided exactly once (duplicate or conflicting
   voids FAIL); voiding an intent that a preceding CONFIRM already completed
   FAILS (a reconciliation must not consume or substitute for a confirmation),
   and a CONFIRM naming a voided intent FAILS. The unmatched-intent invariant
   is GLOBAL (v5): no branch — including branches other than the one under
   test — may hold more than one unvoided, unconsumed intent; only one valid
   terminal intent per branch may remain pending.

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
   and the PR must carry a FRESH approving review from a non-author reviewer applied
   to the EXACT declared head and submitted before the merge (v5: latest review state
   per reviewer; dismissed or superseded approvals do not count); the PR's
   `merge_commit_sha` must equal the protected branch's live candidate, the source
   ref, merged state, and `merged_at` must match the declaration, and the merge
   PushEvent's `before` must equal the declared base. Raced, unrelated, unapproved,
   stale-approved, or API-unverifiable updates FAIL. The merge must be executed as a
   MERGE COMMIT through the protected PR path — never a direct push, never a bypass.
   On the source branch, prior merge-intents on the same lineage are SUPERSEDED via
   their introducing commits; anything else must be voided.

7. **PUSH-NOT-PROVEN (v5, direction 3).** The honest reclassification of a
   historical confirmation whose branch push cannot be proven:

       ## <UTC> | P-<n> | push-not-proven | reclassifies=<confirm-entry-id> branch=<name> remote_sha=<40hex> credential=<name> review_ref=<canonical-record-ref> evidence="<what was searched, what was found, why proof is insufficient>"

   Commit existence may be proven; push provenance is not. The reclassified entry
   stays visible in the log, is disclosed as NOT PROVEN on every verification run,
   and is never counted as evidence for any declaration. `review_ref` MUST resolve
   into the canonical auditor-adjudicated record (repo-relative path + `#` heading
   anchor) whose section mentions the reclassified entry — "reviewed" is never
   merely a string written by the same record's author. Exactly one
   reclassification per confirmation; reclassifying a non-CONFIRM/RECONCILE entry
   or a mismatched branch/SHA fails. This record never edits the original entry
   (append-only).

8. **EVIDENCE-BACK claim separation (v5 hardening of the v4 record).** The
   evidence-back grammar is unchanged, but the checker now verifies three
   SEPARATE claims: (i) commit existence (object store or commits API);
   (ii) push provenance — a push-event Actions run (event=push, head_branch,
   head_sha), a live Events-API PushEvent (ref + before/after), a RECORDED
   `push_event=<id> ref=<ref> before=<40hex|0{40}> after=<40hex>` citation,
   or a RECORDED prior-tip citation `push_event=<id> ref=<ref>
   prior_tip=<40hex>` (a LATER push on the same ref whose `before` field
   records the backed SHA — GitHub's own record that the ref pointed there;
   used when the push's own event rolled off the retained window and the
   feed elides its `after`) — each adjudicated by a valid `review_ref=` on
   the back; and (iii) CI success —
   only successful runs count, with failed/cancelled runs disclosed separately
   (a failed run may evidence the push event itself but is never successful
   gate evidence). A commits-API HTTP 200 or a free-text rebuild narrative is
   NOT push provenance; without provenance or a push-not-proven
   reclassification the gate fails closed (local no-source tool mode discloses
   the gap and never silently passes it).

## Remote verification (v3/v4/v5, `--remote-verify`)

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
  declare the verified pre-push tip, with all excess unmatched intents voided — v5 makes
  the unmatched-intent invariant global across ALL branches), by an
  exact-tip CONFIRM/RECONCILE, or — on a protected branch — by a VERIFIED merge-intent. A
  tip whose push was never declared FAILS (an ls-remote failure or an absent remote ref
  also fails closed — missing remote verification is never a pass);
- **(v5, direction 4) that the terminal intent's declared pre-push state is bound to
  TRUSTED push-event evidence** when the events source is configured (the CI gate always
  configures it). ADAPTED to this repository's live feed (verified 2026-10-20: the feed
  ELIDES PushEvent.after — only `before` is reliable — and ref updates made through the
  git-data API generate NO PushEvent at all), the accepted evidence, in preference
  order: (1) a PushEvent on the branch's ref with `after` == the candidate AND `before`
  == the declared `expected_remote_tip_before_push` (all-zero for an ABSENT
  declaration); (2) a PushEvent on the ref whose `before` == the declared pre-push tip
  (after elided) PLUS a verified push-event Actions run on the exact candidate
  (event=push, head_branch, head_sha); (3) a verified push-event Actions run on the
  exact candidate alone (the API-driven-ref-update class, disclosed). An ABSENT
  declaration is backed by the ref's absent prior state (the all-zero before), not
  merely by the absence of older log records. A prior unlogged ref update, a false
  ABSENT, a raced update, no push-event evidence at all, or a PR-event run substituting
  for push-event provenance FAILS;
- **(v5, direction 3) that every cited run's EVENT TYPE and HEAD_BRANCH are validated
  against the particular claim**: a push-event run on a different branch contradicts the
  entry and FAILS; push-event provenance and CI success are verified and reported as
  separate claims (a failed or cancelled run is disclosed and can never be represented as
  successful gate evidence).

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
