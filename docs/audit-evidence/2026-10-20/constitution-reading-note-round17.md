# Constitution Reading Note — Round 17 (2026-10-20)

**Reader:** the coder (this session), before any edit this round.
**Scope:** Round-17 directions 1–6: reconcile concurrent-writer moves; build
`release/integration-to-main` from `integration/p0-hardening-2026-10-19`
@ `0eae6d5…`; dispatch-only deploy; Gemini-wiring removal; gate on the exact
SHA; merge to `main`; rebase + gate + merge the docs branches.

## Documents read, in full, before any edit

| Document | Blob (re-verified this round) | Where verified |
|---|---|---|
| `constitution.md` (domain, v6.1-Hardened) | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | `origin/main` ls-tree |
| `ENGINEERING_CONSTITUTION.md` (engineering, 1.0.0) | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | `origin/policy/engineering-constitution` and `origin/governance/canonical-docs` ls-tree (identical) |

Both were read end-to-end from the pinned blobs, not from memory or from a
working tree of uncertain provenance. Amendment proposals **EC-16 through
EC-26** (branch `proposal/ec16-ec26-enforcement`, now at `a0e321c…`) are
**PROPOSALS AND NOT IN FORCE**; no rule from that package is enforced or
relied upon in this round.

## Domain rules in force (read, and untouched this round)

- **OBS-153** (ACTIVE): bipodal success measured by Contact Density
  ρ(Bipod) ≥ 2.5 × ρ(Monopod), unique heavy-atom contacts within 4.5 Å.
- **CRIT-154** (ACTIVE): Two-Stage Warhead-RMSD < 1.0 Å vs clinical
  benchmark PDB 5NIU; veto on non-physical deformation (steric clash < 1.2 Å).

Round 17 performs **no verdict, threshold, cohort, or law-semantic work**
(direction 6). The deterministic engine that computes OBS-153/CRIT-154
outcomes is not modified. Hash locks asserted by the gate (canon
`6a9cd4b4349b81de`, matrix `ed1731ace09a…`) must remain bit-identical —
they are asserted on the exact release SHA in CI.

## Engineering rules enforced this round (EC 1–15, in force)

- **EC 1 (purpose/authority):** conflicts are recorded, not silently
  resolved. The Round-16 audit finding (stale reported tips; concurrent
  merges carrying unevaluated content) is recorded in the reconcile report
  with per-commit file lists, not argued away.
- **EC 2 (pre-coding protocol):** exact base SHA confirmed
  (`0eae6d5…`, fetched and matched against `origin`); both constitutions
  read; the deploy workflow, gate workflow, requirements, lock, call sites
  (`backend/main.py`, nested legacy `main.py`), consumers
  (`pdf_generator.py` — zero references to `witness_reports`/
  `ai_model_used`; `dashboard/dashboard.py` — `.get()` with built-in
  degraded state), and deployment files (compose, helm) were all read
  before editing.
- **EC 3 (security):** the provider credential removed from the tree is
  treated as compromised until its issuer confirms revocation (owner action
  A-02, open). No secret is committed, echoed, or logged by this round's
  changes; the replacement deploy workflow references no secrets at all.
- **EC 5 (architecture/dependencies):** import graph checked before the
  compiler deletion (`enrichment/__init__.py` empty; only `backend/main.py`
  and the nested legacy `main.py` import the module; the dossier generator
  never reads its output). `google-auth`/`google-genai` pins removed only
  after confirming no other lock consumer.
- **EC 6 (testing/evidence):** local acceptance recorded with exact
  commands, SHA, and counts (below); CI evidence reported by run ID read
  from the Actions API.
- **EC 7 (CI/CD):** the gate must run on the exact proposed merge candidate;
  the replacement deploy workflow has **no push trigger and no
  `workflow_run` trigger**, so no push or gate completion can start a
  deployment — an independent push-triggered deploy can no longer bypass
  the quality gate. No gate step is weakened (EC 13).
- **EC 11 (change/release discipline):** deployment controls are a
  high-risk surface; the deploy replacement is fail-closed (the manual job
  performs no deployment and references no hosts or secrets).
- **EC 12 (definition of done):** base and proposed SHAs recorded here and
  in the commit message; the final report lists evidence, blocked checks,
  and unverified items.
- **EC 13 (prohibited shortcuts):** no gate weakened or bypassed to obtain
  green; no merging a commit that did not pass its checks (directions 3–5
  gate every merge candidate on its exact SHA before merge).
- **EC 14 (task artifacts):** this note, the commit message manifest, and
  the gate run records on the exact SHAs are the task artifact set.

## Single-writer discipline (auditor decision, this round)

Before every push: fetch; read the current tip; push only from the tip just
read. If another writer has moved the branch, stop and report. No merges of
unevaluated content. This is applied to `release/integration-to-main`, to
`main` (merge), and to both docs branches (rebase + merge).

## Binding precedent applied

The Gemini removal ports the **accepted H-1 ruling** (auditor direction,
Round 14; accepted on `hygiene/h1-remove-gemini` @ `5149435…`, gate run
39376518935 GREEN): the narrative compiler is **deleted, not replaced**
(the dossier generator provably never reads `witness_reports`; the
dashboard renders its built-in pending state; response fields stay
optional/backward-compatible). The same file set is swept here on the
integration lineage.

## Local acceptance evidence (this branch, pre-push)

Recorded in the commit message with exact commands and counts: acceptance
grep (`git grep -nE "GEMINI|generativeai|AIzaSy"` outside evidence/security
folders → empty), stricter case-insensitive sweep → empty, backend pytest
counts vs the pristine integration base, and canon/matrix hash verification.
