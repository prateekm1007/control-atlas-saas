# Incident record — P-015-CONFIRM identifier collision + force-update deviation (2026-10-20)

**Recorded:** 2026-10-20 (Round-19 repair, auditor direction 4).
**Kind:** process deviation + evidence-identifier collision. No data loss; no
secret exposure; no main history rewrite; all prior entries preserved unchanged.

## 1. What happened

During Round-18, two writers using the same sandbox git identity (`Z User
<z@container>`, shared container; attribution discipline per A-09 — git author
metadata alone does NOT establish which credential or person executed a push)
each executed the same auditor directive ("rebuild the proposal package onto
the release main") in parallel, on different branches, each carrying its own
copy of `docs/evidence/push-log/PUSH_LOG.md`:

- **Writer A** (proposal branch, unmerged): rebased `proposal/ec16-ec26-enforcement`
  from `a0e321c` to `c4a0fcc` and pushed with **`--force-with-lease`**
  (`a0e321c` → `c4a0fcc`, remote update observed by ls-remote; commit
  timestamps `2026-10-20T11:38:48Z`).
- **Writer B** (new branch, merged): created
  `proposal/ec16-ec26-enforcement-rebased` @ `0f1e9a8` (11:27:20Z), pushed it
  (11:27:51Z), appended `P-015-CONFIRM` in `3745cd2` (11:28:06Z), and merged
  into `main` as `f535a71` (11:32:26Z) via PR #17 (gate run 39415713844 =
  success, but the gate runs no round-specific acceptance checks — the defects
  below passed unnoticed; see direction 7).

## 2. The identifier collision (P-015-CONFIRM, one ID — two meanings)

| History | Entry | Records | Remote SHA cited |
|---|---|---|---|
| `main` (via `proposal/ec16-ec26-enforcement-rebased` @ `3745cd2`, PR #17) | `P-015-INTENT` + `P-015-CONFIRM` | push of the NEW branch `proposal/ec16-ec26-enforcement-rebased` | `0f1e9a8856bafec592c657270e5799734ed61d6a` |
| `proposal/ec16-ec26-enforcement` @ `c4a0fcc` (unmerged, preserved as evidence) | `P-015-CONFIRM` | confirmation of the EARLIER tip of the ORIGINAL branch | `a0e321c0ef08cd3430e7910936075994791fd8df` |

Both entries are individually truthful about the push each describes, but they
reuse the ENTRY-ID `P-015-CONFIRM` for different branches/SHAs. Merging the two
logs would create a duplicate identifier and violate the "one ID, one meaning"
discipline; the hardened `evidence_sequence_check.py --pushlog-schema` (this
repair) now rejects exactly this class.

## 3. The force-update deviation (recorded honestly, not smoothed over)

Register rev 7 Part 8 authorized "continued focused, auditable corrections on
`proposal/ec16-ec26-enforcement` — **appended commits only, no force-push**".
The reported update of that branch from `a0e321c` to `c4a0fcc` was executed as
a **force-with-lease update** (rebase of the branch's own commits onto
`eb6d6a9`), which is a force-push of a published branch and therefore a
deviation from that authorization. Facts recorded:

- Before: `a0e321c0ef08cd3430e7910936075994791fd8df` (tip since 10:32:13Z).
- After: `c4a0fcc3357930b7eba61e056489c804b16d2277` (commit time 11:38:48Z;
  the push itself was issued immediately after; the exact remote acceptance
  moment is not separately evidenced — the branch's own `P-016-INTENT` entry
  discloses the force-with-lease plan, and no `P-016-CONFIRM` was ever
  appended: the two-phase protocol on that branch is INCOMPLETE/PENDING).
- Effect on the audit trail: the pre-rebase commit chain (`a0e321c` and its
  ancestors) remains reachable on the remote ONLY through
  `proposal/ec16-ec26-enforcement`'s reflog server-side and any clones; the
  branch ref no longer points at it, so history-anchored evidence that cited
  "tip of proposal/ec16-ec26-enforcement" as of 11:00–11:38Z now resolves
  differently than at citation time. The earlier tip was preserved locally and
  in this repository's object store; nothing was fabricated.
- No push to `main` occurred from writer A at any point in Round-18
  (single-writer check stopped the merge attempt).

## 4. Resolution (auditor direction 2 + direction 4 — repair forward)

1. **Current `main` (`f535a71`) is the baseline.** It is repaired forward
   through a reviewed, narrowly scoped PR (this branch); it is NOT replaced by
   or merged with `c4a0fcc`. The adjudication is the auditor's Round-19
   direction 2.
2. **Mainline log entries are preserved unchanged.** No existing entry in
   `PUSH_LOG.md` is edited, renumbered, or deleted. The collision is recorded
   HERE (dated, evidence-backed) instead of being erased.
3. **The `c4a0fcc` push-log copy is never merged wholesale.** The reviewed
   content it carried (provider-token removals, README corrections) was carried
   into this repair by content; its push-log entries (`P-015-CONFIRM` for
   `a0e321c`, `P-016-INTENT` unconfirmed) remain only on that branch, preserved
   as evidence of what happened, per direction 3A's "do not copy its push log
   wholesale".
4. **New entries use globally unused IDs.** A census of every published
   branch's `PUSH_LOG.md` (main, both proposal branches, all eight activation
   branches) was taken before selecting IDs: used IDs are `P-001`…`P-011`,
   `P-011-CONFIRM`, `P-012/013/014-INTENT`+`-CONFIRM`, `P-015-INTENT`+
   `-CONFIRM` (mainline), `P-015-CONFIRM`+`P-016-INTENT` (c4a0fcc), and the
   activation-cycle `P-INTENT-ec16/18/19/20/21/23/24/26` forms. New mainline
   numbering therefore starts at **`P-017`**.
5. **Confirmations record exact live remote refs** (`git ls-remote` evidence
   strings), never local commit hashes or workflow success, per direction 4.4.
6. **Pending is called pending.** The activation branches' `P-INTENT-ec<NN>`
   entries have no CONFIRM records; the hardened validator reports them as
   PENDING (not complete, not erased). Their reconciliation happens on their
   own follow-up commits before any of those branches merge.

## 5. Verification commands (repeatable)

    git ls-remote origin | grep -E "ec16-ec26|activation"
    git show c4a0fcc:docs/evidence/push-log/PUSH_LOG.md | grep -n "P-015\|P-016"
    git show 3745cd2:docs/evidence/push-log/PUSH_LOG.md | grep -n "P-015"
    python3 scripts/ci/evidence_sequence_check.py docs/evidence/push-log/PUSH_LOG.md --pushlog-schema
