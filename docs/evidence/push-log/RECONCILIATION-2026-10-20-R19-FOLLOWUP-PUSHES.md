# RECONCILIATION — Round-19 follow-up pushes (Round-20 auditor direction 4, 2026-10-20)

## What happened

During the Round-19 PAT stage, each of the eight Round-19 branches received a
second push that appended its `P-0xx-CONFIRM` entry (completing the two-phase
protocol for the FIRST push) and, for the repair branch, the two workflow
fixes. Those second pushes were **not pre-declared by their own INTENT
records**: the confirmations they carried describe the earlier first-push
tips, not the tips they created.

The Round-20 audit correctly found that the then-current PR heads
(`f44bf4d`, `bd8d86a`, `32d272d`, `151b8c1`, `97067a0`, `fa854f8`,
`28ee983`, `df98f5e`) were therefore **unreconciled in the evidence**: the
confirmations were valid for the pushes they described, but nothing declared
or confirmed the pushes that carried them. The Round-19 report's claim that
the push protocol was fully closed was not supported.

## What this record does (and does not do)

- The eight `push-reconcile` entries **P-050 … P-057** (appended to
  `PUSH_LOG.md` in this same commit) reconcile those pushes truthfully: dated
  at the reconciliation window `2026-10-20T15:58:48Z`, with live `ls-remote`
  re-verification of every tip performed at that window, and check-run
  citations for each exact head.
- **Nothing is backdated.** No INTENT is manufactured for a push that already
  happened; no prior entry is edited or deleted (the append-only
  parent-prefix check stays in force).
- The protocol gap itself is recorded here rather than smoothed over. The
  `push-reconcile` verb is defined in `PUSH_LOG_SCHEMA.md` (schema v3) and is
  validated by `scripts/ci/evidence_sequence_check.py` (AT-24): it carries
  the same required fields as a CONFIRM, never pairs with or consumes INTENT
  records, and is classified CURRENT/HISTORICAL against the live remote ref
  exactly like a CONFIRM.
- Going forward, the hardened checker **fails closed** on this class: with
  `--remote-verify --branch <name> --expect-tip <sha>`, a branch whose live
  tip is declared by neither an unpaired INTENT (PENDING, disclosed) nor an
  exact-tip CONFIRM/RECONCILE FAILS the check. The Governance Acceptance
  workflow's `push-log-integrity` job now runs exactly that verification on
  every push and pull-request event.

## Reconciled pushes (all eight verified live at 2026-10-20T15:58:48Z)

| Entry | Branch | Reconciled tip | Checks cited on the exact head |
|---|---|---|---|
| P-050 | `repair/round19-mainline-acceptance` | `f44bf4d93f8d69a4bd651396f383e8f3ad2b707a` | gate run 39443515998; governance-acceptance run 39443516042 (PR-event, auditor-cited) |
| P-051 | `activation/ec16-ec19-ec20-r19` | `bd8d86a3a93c1717921587f9585f2b97baf22836` | push-event run 39443519540 (Round-19 report) |
| P-052 | `activation/ec21-r19` | `32d272d3243f62224248c9d29eb5dcb0c8a20e5a` | push-event run 39443538818 (Round-19 report) |
| P-053 | `activation/ec26-r19` | `151b8c12d6e5286cc3153461677b193b6d57495d` | push-event run 39443554590 (Round-19 report) |
| P-054 | `activation/ec18-r19` | `97067ad0320968bfb73f7e045408fcb60ccd5c7d` | push-event run 39443566175 (Round-19 report) |
| P-055 | `activation/ec23-r19` | `fa854f8f3b7fcff9c83db5570784970f42df9742` | push-event run 39443578786 (Round-19 report) |
| P-056 | `activation/ec24-r19` | `28ee983e0f98723c096ce4b13892f5c16bfbb22e` | push-event run 39443597533 (Round-19 report) |
| P-057 | `docs/owner-decisions-r19` | `df98f5e3671055d67961d5b9d14495c5b037ce19` | push-event run 39443611965 (Round-19 report) |

The run citations above are corroborating evidence; the primary evidence for
each entry is the live-ref verification quoted in its `evidence=` field and
re-performed by the checker on every run. Where a runs manifest is available
(CI), a cited run whose `head_sha` contradicts the recorded SHA fails the
check (AT-24 case 9).

## Identifier census

New identifiers P-050…P-058 were chosen after a fresh census of every
published remote branch's `PUSH_LOG.md` at the reconciliation window: the
highest previously used identifier is **P-049** (on `fix/mainline-r18`, the
concurrent writer's series). No identifier is reused across any published
log (the P-015-CONFIRM collision class does not recur).
