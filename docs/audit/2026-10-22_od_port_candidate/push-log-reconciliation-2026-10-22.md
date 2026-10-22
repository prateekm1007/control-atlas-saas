# Push-log reconciliation and governance-run record - Round-26 review branch (2026-10-22)

## 1. Why this record exists

The round-26 candidate `e2cd17d` (audit package: LAw-182 per-chain fix +
three-engine C5 evidence regeneration + port rationale + OD-20 conflict
record + OD-14 dossier) was pushed without appending any push-log rows. The
Governance Acceptance checker (EC-18/EC-24) on that tip reported the branch
UNDECLARED: "the live tip of review/od-scope-v25-2026-10-22 (e2cd17d0f9c5...)
is UNDECLARED - no unpaired INTENT declares its push and no CONFIRM/RECONCILE
matches it". In fact, NO push on this branch was ever declared in
`docs/evidence/push-log/PUSH_LOG.md`: the branch was created with a wrong head,
force-corrected, and advanced, all without log rows (the port rationale, section
2, discloses the correction). This push appends the missing rows in the same
push they describe, per the log's two-phase protocol.

## 2. Every push on this branch (now declared)

| # | Push (full SHA) | Time (UTC) | Ref transition | Declared by |
|---|---|---|---|---|
| 1 | `97d62a89f0efb2917e80a82808a930d34af7f492` (delivery tip, wrong refspec) | 2026-10-21T21:16:46Z (RECONSTRUCTED, bound by run creation) | ref created -> 97d62a8 | P-082 |
| 2 | `63da0e347002ed9b77071ff1ed0cefda7c0ce88a` (audited base) | 2026-10-21T21:17:01Z (RECONSTRUCTED, bound by run creation) | 97d62a8 -> 63da0e3 (force; only this ref touched) | P-083 |
| 3 | `e2cd17d0f9c5231c489c42b8dd90d74055d1f648` (round-26 audit package) | 2026-10-22T12:34:27Z (bound by run creation) | 63da0e3 -> e2cd17d | P-084 (push-reconcile, live-verified at write time) |
| 4 | `9c4bb1637474c199d1238cec10b35a6ca88d1efb` (reconciliation rows + run record) | 2026-10-22T18:16:37Z (bound by run creation) | e2cd17d -> 9c4bb16 | P-085-INTENT (two-phase, deferred) |
| 5 | `37698cd448cc351a25a7340d34599dac6b4a4f45` (SHA256SUMS digest refresh) | 2026-10-22T18:19:37Z (bound by run creation) | 9c4bb16 -> 37698cd | P-085-CONFIRM (records the two-phase pair opened by P-085-INTENT; the refresh push is the legitimate push that carries the confirm) |
| 6 | (the tip this push creates) | 2026-10-22 | 37698cd -> (this push's tip) | P-086-INTENT (two-phase; P-086-CONFIRM lands with the next legitimate push of the branch) |

## 3. Verification runs cited by the appended rows

| Run ID | Workflow | Head SHA | Event | Conclusion |
|---|---|---|---|---|
| 39677471359 | Toscanini Stability Gate | 97d62a8 | push | success |
| 39677508627 | Governance Acceptance | 63da0e3 | push | failure (pre-existing dead-ref set below + undeclared tip at that time) |
| 39677508732 | Toscanini Stability Gate | 63da0e3 | push | success |
| 39779457004 | Governance Acceptance | e2cd17d | push | failure (22 pre-existing dead-ref entries + 1 new undeclared-tip finding) |
| 39779457093 | Toscanini Stability Gate | e2cd17d | push | success |
| 39829905845 | Governance Acceptance | 9c4bb16 | push | failure (22 pre-existing dead-ref entries + 1 package-manifest digest mismatch) |
| 39829905992 | Toscanini Stability Gate | 9c4bb16 | push | success |
| 39830332863 | Governance Acceptance | 37698cd | push | failure (22 pre-existing dead-ref entries + terminal-intent binding failure carried over from the undeclared 9c4bb16 push) |
| 39830332865 | Toscanini Stability Gate | 37698cd | push | success |

The runs for THIS push's tip are recorded after it completes; their IDs are
available in the workflow history of the tip that P-086-INTENT declares
(this record was updated before those runs concluded).

## 4. Pre-existing EC-18/EC-24 dead-ref findings (NOT introduced by round-26)

The checker verifies every confirmed push against the LIVE remote ref. Branches
that were pruned post-merge (per the sanctioned Round-23 branch prune, and the
earlier round-19/20 branch cleanups) no longer exist, so their historical
CONFIRM/RECONCILE rows fail closed. Twenty-two such rows fail on EVERY governance
run of this branch, including the pre-port base commit `63da0e3` (run
39677508627, 2026-10-21):

| Line | Entry | Branch (pruned/deleted) |
|---|---|---|
| 51 | P-011-CONFIRM | proposal/ec16-ec26-enforcement |
| 88 | P-012-CONFIRM | docs/owner-decisions |
| 106 | P-013-CONFIRM | proposal/ec16-ec26-enforcement |
| 115 | P-014-CONFIRM | governance/canonical-docs |
| 137 | P-015-CONFIRM | proposal/ec16-ec26-enforcement-rebased |
| 161 | P-017-CONFIRM | repair/round19-mainline-acceptance |
| 171 | P-050 | repair/round19-mainline-acceptance |
| 176 | P-051 | activation/ec16-ec19-ec20-r19 |
| 181 | P-052 | activation/ec21-r19 |
| 186 | P-053 | activation/ec26-r19 |
| 191 | P-054 | activation/ec18-r19 |
| 196 | P-055 | activation/ec23-r19 |
| 201 | P-056 | activation/ec24-r19 |
| 206 | P-057 | docs/owner-decisions-r19 |
| 218 | P-058-CONFIRM | repair/round19-mainline-acceptance |
| 251 | P-059-CONFIRM | repair/round19-mainline-acceptance |
| 340 | P-060-CONFIRM | repair/round19-mainline-acceptance |
| 379 | P-066-CONFIRM | repair/round19-mainline-acceptance |
| 421 | P-069-CONFIRM | repair/round19-mainline-acceptance |
| 444 | P-070-CONFIRM | repair/round19-mainline-acceptance |
| 481 | P-073-CONFIRM | repair/round19-mainline-acceptance |
| 542 | P-078-CONFIRM | repair/round19-mainline-acceptance |

These are out of scope for round-26: they pre-date this candidate (the base run
39677508627 failed the identical set, before any round-26 change existed), and
repairing them would require either resurrecting pruned refs or editing the
checker/workflow, both outside the OD-16 hard-stop mandate. They are disclosed,
not repaired.

## 5. What the appended rows do under the checker

- **P-082 / P-083** - legacy `push` rows (RECONSTRUCTED): bookkeeping entries
  for the two pushes that preceded the round-26 package. No pairing
  obligation.
- **P-084** - `push-reconcile` for `e2cd17d`. At the time of this push's CI
  run, the live tip will be the new push (T5), so P-084's recorded SHA
  (`e2cd17d`) is on the branch's current lineage: it classifies
  SUPPORTED-BY-ANCESTRY (a confirmation on an ancestor of the live tip is
  valid for the push it describes).
- **P-085-INTENT** - terminal intent for this push (DEFERRED two-phase).
  The checker binds it when: `expected_remote_tip_before_push` equals the
  newest recorded branch state (`e2cd17d`, P-084) and the intent's
  introducing commit is the exact candidate under test. Both hold by
  construction. Its CONFIRM is delivered by the next legitimate push of the
  branch (PENDING-disclosed, never backdated).

## 6. Scope statement

This push modifies ONLY: `docs/evidence/push-log/PUSH_LOG.md` (append), this
record file, and one port-rationale subsection. No source, checker,
governance, registry, or test content is modified. No merge, deploy,
calibration, held-out, threshold, cohort, or golden change is made (OD-16).
The C5 evidence files and engine behavior are unchanged.
