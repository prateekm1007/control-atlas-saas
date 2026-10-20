# Push log schema v3 — two-phase records + post-hoc reconciliation + remote verification

(Defined 2026-10-20, EC-18/EC-24 auditor revisions; v3 adds the `push-reconcile`
verb and the remote-verification layer per the Round-20 auditor direction 4.)

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

       ## <UTC> | P-<n>-INTENT | push-intent | branch=<name> head_at_intent=<40hex> credential=<name> authorized_by=<reference>

2. **CONFIRM (post-push).** Appended in a LATER commit, after the push returned, with its own
   truthful `date -u` timestamp and the remote evidence source:

       ## <UTC> | P-<n>-CONFIRM | push-confirm | branch=<name> remote_sha=<40hex> credential=<name> evidence="git ls-remote origin refs/heads/<name>"

   `remote_sha` MUST be the full 40-hex SHA observed on the remote (placeholder markers such as
   `self-referential-tip` are rejected by `--pushlog-schema`). CONFIRM entries for a branch's
   own tip are published by a SUBSEQUENT push of that branch; the very last confirm of a
   session is recorded in the round report/worklog and appended to this log by the next
   legitimate push (disclosed limitation, never backdated).

3. **RECONCILE (post-hoc, v3).** For pushes that already happened without a preceding INTENT
   — the Round-19 follow-up-push class found by the Round-20 audit — a truthful, dated record
   that reconciles the push without pretending an INTENT existed:

       ## <UTC> | P-<n> | push-reconcile | branch=<name> remote_sha=<40hex> credential=<name> evidence="git ls-remote <url> refs/heads/<name> -> <sha> (re-verified <UTC>); <run citations>"

   A RECONCILE never pairs with or consumes INTENT records (it is post-hoc by design), never
   backdates, and its evidence must contain both the recorded SHA and the branch name, plus
   the live-ref verification performed at reconciliation time and (when available) check-run
   citations for the exact head.

## Remote verification (v3, `--remote-verify`)

`scripts/ci/evidence_sequence_check.py --remote-verify [--branch <under-test>
[--expect-tip <sha>]] [--runs-manifest FILE]` verifies, fail-closed:

- every CONFIRM/RECONCILE's recorded `remote_sha` against the **live remote ref** for its
  claimed branch (`git ls-remote`): equal => the confirmation is CURRENT; different => it is
  HISTORICAL — valid for the push it describes, and explicitly NOT a confirmation of the
  newer tip;
- that the `evidence=` field actually **supports** the recorded SHA: the SHA and the claimed
  branch must appear in the evidence (or the entry body), and any cited Actions run ID whose
  manifest record shows a different head_sha FAILS (fabricated support);
- that the branch under test's CURRENT tip is DECLARED: by an unpaired INTENT (reported
  PENDING — the honest steady state for the newest push of a session) or by a exact-tip
  CONFIRM/RECONCILE. A tip whose push was never declared FAILS (an ls-remote failure or an
  absent remote ref also fails closed — missing remote verification is never a pass).

## Truthfulness rules

- INTENT and CONFIRM are separate entries with separate IDs; an INTENT entry never claims
  post-push knowledge, a CONFIRM entry never rewrites its INTENT.
- RECONCILE records are visibly post-hoc: they disclose the missing INTENT rather than
  manufacturing one, and they never backdate a push time.
- Timestamps come from `date -u` at the moment the record is written. Historical entries
  whose timestamps could not be observed at action time are labeled `RECONSTRUCTED`.
- Prior entries are never edited or deleted (enforced by the parent-prefix check against
  the previous commit's version of the log).
- A CONFIRM entry reconciles (in its body) any difference between the INTENT/entry timestamp,
  the git commit metadata, and the remote observation; differences are explained, not averaged
  away, and unprovable claims are marked as such.
- Credential identifiers are recorded by NAME only; token values are never written anywhere.
