# Push log schema v2 — two-phase records (defined 2026-10-20, EC-18/EC-24 auditor revisions)

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

## Truthfulness rules

- INTENT and CONFIRM are separate entries with separate IDs; an INTENT entry never claims
  post-push knowledge, a CONFIRM entry never rewrites its INTENT.
- Timestamps come from `date -u` at the moment the record is written. Historical entries
  whose timestamps could not be observed at action time are labeled `RECONSTRUCTED`.
- Prior entries are never edited or deleted (enforced by the parent-prefix check against
  the previous commit's version of the log).
- A CONFIRM entry reconciles (in its body) any difference between the INTENT/entry timestamp,
  the git commit metadata, and the remote observation; differences are explained, not averaged
  away, and unprovable claims are marked as such.
- Credential identifiers are recorded by NAME only; token values are never written anywhere.
