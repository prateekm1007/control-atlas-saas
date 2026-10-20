# Push-log reconciliation — collision and supersession record (Round-18, 2026-10-20)

This document explains the two competing push-log histories and the superseded SHAs that
the push-log protocol checker
(`scripts/ci/evidence_sequence_check.py --pushlog-schema --remote-slug …`) verifies against
the machine-readable record in `PUSH_LOG_RECONCILIATION.json` (same directory).

## 1. The two histories (both preserved; neither copied onto the other)

| History | Where it lives | Push-log content | P-015-CONFIRM meaning |
|---|---|---|---|
| Released mainline | `main` @ `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` (merged from `proposal/ec16-ec26-enforcement-rebased` @ `3745cd2af2e0b5b909cb4c3b49c515afdec4cce8`) | P-001…P-015 (INTENT/CONFIRM pairs for the Round-11→17 pushes) | confirms `proposal/ec16-ec26-enforcement-rebased` @ `0f1e9a8856bafec592c657270e5799734ed61d6a` |
| Divergent proposal history | branch `proposal/ec16-ec26-enforcement` @ `c4a0fcc3357930b7eba61e056489c804b16d2277` (advanced by the "round-18 adjustment" push disclosed in its own `P-016-INTENT`) | its own copy of P-001…P-014 plus its own `P-015-CONFIRM` and `P-016-INTENT` | confirms `proposal/ec16-ec26-enforcement` @ `a0e321c0ef08cd3430e7910936075994791fd8df` |

The same entry ID `P-015-CONFIRM` denotes two different confirmations in two different
lineages. Per the Round-18 audit direction: the two records are NOT combined, the
`c4a0fcc` push log is NOT copied onto `main`, and both histories remain untouched.
The checker validates the log on the CURRENT history and never merges records across
histories.

## 2. Superseded SHAs (verified on the remote, but no longer on the current branch lineage)

Each row below was verified via the GitHub commits API on 2026-10-20 (HTTP 200 — the object
exists on the remote) and via `git merge-base --is-ancestor` against the current remote tip
(NOT an ancestor; the objects are not even in a fresh full clone). They are classified
`SUPERSEDED-BY-DOCUMENTED-REBUILD` with a pointer to the documenting evidence; nothing was
rewritten.

| Entry | Branch | Recorded SHA | Current remote tip | Why not on lineage | Pointer |
|---|---|---|---|---|---|
| P-006 | `docs/owner-decisions` | `0dcc5a10756871f2b2dd239b1858cf8aae6de41a` | `afd445cc3c8eafe7950a04df95e4f16b303b194a` | divergent parallel-writer register superseded by the reconciliation merge (`b421237`) and the Round-17 rebuild | `INCIDENTS.md` INC-2; `PUSH_LOG.md` P-012-CONFIRM body; `docs/directives/INDEX.md` SERIES-2026-10-20-R17 row |
| P-010 | `docs/owner-decisions` | `0773f0b1ca698829ae20b59660798dc8b6504675` | `afd445cc3c8eafe7950a04df95e4f16b303b194a` | register rev 6 tip; superseded when the register was rebuilt and merged to `main` as `eb6d6a9` (Round-17 direction 5) | `docs/directives/INDEX.md` SERIES-2026-10-20-R17 row |
| P-012-CONFIRM | `docs/owner-decisions` | `78b147d55660dc79c2f5aada368d8499339bc339` | `afd445cc3c8eafe7950a04df95e4f16b303b194a` | same rebuild event as P-010 | `docs/directives/INDEX.md` SERIES-2026-10-20-R17 row |
| P-011-CONFIRM | `proposal/ec16-ec26-enforcement` | `17eb88184c9f461c6bc9ee592f916f7ff4c921c4` | `c4a0fcc3357930b7eba61e056489c804b16d2277` | old proposal history replaced by the D-4 rebase onto `eb6d6a9` (rebased branch merged as `f535a71`); the original branch was later advanced to `c4a0fcc` by the push disclosed in its own `P-016-INTENT` | `INCIDENTS.md` INC-2; `docs/directives/INDEX.md` SERIES-2026-10-20-R17 row; the `c4a0fcc` history's own push log |
| P-013-CONFIRM | `proposal/ec16-ec26-enforcement` | `dab33d946fe31891efbd6e4917e360abcb1872c1` | `c4a0fcc3357930b7eba61e056489c804b16d2277` | same replacement event as P-011-CONFIRM | as P-011-CONFIRM |
| P-014-CONFIRM | `governance/canonical-docs` | `b12f30e95255b0c5d31a70cb9b6bb42776aac9a6` | `d51ad7c7fcbe521c116b8f32be4904ae6a57ab41` | old governance history superseded by the Round-17 docs-only rebase (`d51ad7c` merged to `main` as `b851e01`) | `docs/directives/INDEX.md` SERIES-2026-10-20-R17 row |

## 3. Rules going forward (binding)

1. New push-log events use globally unused IDs: P-017 upward (highest ID used anywhere in
   any preserved history is P-016, on the `c4a0fcc` history). Incident IDs use `INC-`.
2. A branch whose history must be rebuilt is NEVER force-pushed; the rebuild goes to a NEW
   branch name (the D-4 pattern), and the old branch is preserved as historical evidence.
3. A push INTENT without a matching CONFIRM is PENDING — never silently treated as
   complete — until the confirm is appended by a later legitimate push, or the intent is
   confirmed live by the checker against the actual remote ref, or the intent is reconciled
   here with a pointer to documenting evidence.
