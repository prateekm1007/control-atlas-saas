# Incident Record — Unexplained Remote Pushes Under This Session's Git Identity (push-2)

**Incident ID:** incident-push-2
**Recorded:** 2026-10-20 (Round 14, direction 1)
**Severity:** credential/authorship anomaly + coder report error
**Status:** recorded; owner token-history review requested (A-09)

## 1. What was observed

1. The Round-14 audit found `feature/od07-audit-id` on GitHub at
   `08a9403a7dd376da08831d14f1fa4fa26b15430e`, pushed ~2026-10-20
   05:07 UTC, author/committer `Z User <z@container>` — the same git
   identity this coder session uses. My Round-13 report had stated that
   all pushes were blocked and had reported a different OD-07 tip
   (`73ef84c1…`) as local-only.
2. During this session's fresh fetch, TWO MORE remote branches appeared
   that did not exist at the audit's check:
   - `feature/n8-altlocs-per-residue` @
     `b5afd6ca01a345432852d5865a79999e464636b4` (commits `eb7e565`
     authored 04:11:14 UTC and `b5afd6c` authored 05:14:45 UTC, both
     `Z User <z@container>`)
   - `evidence/round11-fix-2026-10-20` @ `8301c53` (authored 04:21:10
     UTC, `Z User <z@container>`), carrying a second, divergent
     register text and test evidence.
   The remote therefore gained content from a writer using this
   session's git identity while my session was running and after it
   ended.

## 2. Verified facts (this session, fresh fetch)

- `08a9403` is based on `science/p0-correctness-2026-10-19` @
  `ebc1fd79` and contains TWO commits (`0c97094` reading
  note/preflight; `08a9403` OD-07 implementation). Its implementation
  keeps `governance.audit_id = sha256(coord_hash)` with canon/matrix in
  a separate `governance.audit_identity` object — the deviation the
  auditor identified (same coordinates → same ID under different
  canon/matrix).
- None of these remote commits exist in this repository's reflog; they
  entered only via `git fetch`. Their content (file names
  `test_od07_audit_identity.py`, `docs/audit-evidence/round12/`;
  functions `compute_legacy_audit_id`, `resolve_audit_ref`) matches
  neither my reported work (`backend/tos/governance/audit_identity.py`,
  `test_ids.py`) nor my reported SHAs.
- This session's reflog (all times UTC, from the Round-13 session):
  - `hygiene/h1-remove-gemini`: branch created 04:18:47 from
    `origin/main`; commits 04:23:54 / 04:24:13 (amend) → `a761d72`
  - `feature/od07-audit-id`: branch created 04:32:38 from science;
    commit 04:40:07 → `73ef84c`
  - `feature/n8-altlocs-per-residue`: branch created 04:40:38 from
    `73ef84c`; commit 04:49:45 → `af2c164`
  - (provenance/census/governance/register branches committed
    04:5x–05:0x; worklog last write 05:00; durability bundles 04:59)
- The parallel commits' author dates (04:11:14, 04:21:10, 05:07:10/21,
  05:14:45) interleave with — and in one case precede — my session's
  branch creations. `eb7e565` (04:11:14) predates my first branch
  creation (04:18:47).
- This session issued exactly ONE push command in Round 13
  (`hygiene/h1-remove-gemini`, ~04:5x UTC), which FAILED (see §3). No
  push of `feature/od07-audit-id` or any other branch was issued by
  this session before the Round-14 token re-test.

## 3. Coder report error (root cause, corrected)

My Round-13 report claimed the PAT was rejected ("401 Bad credentials")
and that pushes were therefore blocked. **That claim was wrong.**
Re-testing with the export and the use of the variable inside a SINGLE
shell invocation returns HTTP 200 (`login: prateekm1007`). Root cause:
the Round-13 session exported `CA_PAT` in one shell invocation and
consumed `${CA_PAT}` in separate invocations; the variable did not
propagate, so the push URL and the API header carried an EMPTY
credential, which GitHub rejects as `Bad credentials`. The token itself
was valid throughout. This error blocked my push series; it does NOT
explain the remote pushes, because this session never issued them.

## 4. Credential name used for the remote pushes

Not determinable from this environment. The only credential held by
this session is the owner-provided classic repository PAT (prefix
`ghp_`, value never written to any file, artifact, or log; provided
again by the owner in the Round-14 message with the instruction to push
and not to request revocation). GitHub does not expose token-level
attribution through the public API; the push-event feed showed no
PushEvent entries for today's pushes on the first page queried.
**Attribution requires the owner's token-history review (A-09).** The
account-level actor for the pushes could not be recovered from this
environment either; the commits' author/committer identity is the
sandbox git identity (`Z User <z@container>`), which is not
account-proof.

## 5. Possible explanations (recorded, not concluded)

- A parallel automated session using the same sandbox git identity and
  a working credential, executing Round-12-era directions (the pushed
  commit messages reference "Round-12 Direction 1" and "Direction 3
  (Round-12)").
- The owner or a party authorized by the owner pushing work authored
  under this identity.
- This coder session cannot distinguish between these from available
  evidence and does not guess.

## 6. Handling (per direction 1 and the constitution)

- No force-push; no branch deletion. `feature/od07-audit-id` @
  `08a9403` is preserved exactly as pushed.
- OD-07 compliance is restored by a NEW commit on top of `08a9403`
  (direction 2), making `audit_id` the SHA-256 of the canonical
  serialization binding coordinate hash, canon hash, and matrix hash,
  with the ruling test added and the superseded assertion rewritten.
- The parallel `feature/n8-altlocs-per-residue` was reviewed against
  R-4 (per-residue default, per-atom opt-in, 11 passing policy tests,
  before/after pinned-cohort verdict diff showing zero verdict changes)
  and is left untouched; my independent N8 implementation remains
  preserved as patch evidence.
- The divergent register text on `evidence/round11-fix-2026-10-20` is
  recorded in the canonical register (rev 6) as a divergence; the
  canonical location remains `docs/decisions/OWNER_DECISIONS.md` on
  `docs/owner-decisions`.
- Owner actions requested: A-09 token-history review of the 05:07 and
  subsequent pushes; A-02 issuer-side revocation confirmation for the
  removed provider key.

## 7. Machine-readable annex

See `annex-git-logs.txt` in this directory: the pushed branches' full
commit metadata, the local reflog excerpts, and the fetch output
showing the three new remote branches.
