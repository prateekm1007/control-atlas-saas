# Constitution Reading Note — Round-12 remediation session (2026-10-20)

**Task:** Round-12 audit directions 1–6 (OD-07 audit identity; N8 altloc default;
H-1 Gemini removal; authoritative decision register; evidence-README corrections).

## 1. Documents read in full, before any code edit

| Document | Source | Git blob (verified) | Extent |
|---|---|---|---|
| `ENGINEERING_CONSTITUTION.md` (SaaS Engineering Constitution v1.0.0) | branch `policy/engineering-constitution` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | all 178 lines |
| `constitution.md` (Toscanini MDI Constitution v6.1-Hardened) | same branch; identical blob present at `ebc1fd7` root | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | all 11 lines |

Both blobs were resolved with `git ls-tree` / `git cat-file` from the fetched
repository and match the SHAs cited in the audit direction. Both documents were
read in full before the first code edit of this session.

Per EC §1 this note records the rules relevant to the task; per EC §2 the
session preflight is `PREFLIGHT_OD07_2026-10-20.md` in this directory.

## 2. EC IDs enforced this session, and how

- **EC §1 (Purpose and authority)** — Both constitutions read; `constitution.md`
  governs domain/scientific claims, the Engineering Constitution governs this
  engineering work. Where the audit's earlier claims conflicted with observed
  repository state, the conflict is recorded (evidence-branch README
  correction, Direction 6), not silently reconciled. Ambiguity in OD-07's
  "combined with" phrasing is resolved conservatively (identity triple, no
  change to the SHA-256 formula) and the interpretation is stated in the
  commit, not guessed silently.
- **EC §2 (Mandatory pre-coding protocol)** — Base commit recorded
  (`ebc1fd79…`, science branch `science/p0-correctness-2026-10-19`; tree clean
  in a fresh worktree); current implementation, contracts, tests and callers
  read before editing (`backend/main.py`, `backend/tos/ingestion/processor.py`,
  `backend/tos/storage/audit_store.py`, `backend/tos/schemas/response_v1.py`,
  `backend/tos/engine/batch_processor.py`, `backend/tests/test_ids.py`);
  machine-readable preflight recorded (see preflight file).
- **EC §3 (Security is a release requirement)** — *Authorize objects, not just
  requests*: the OD-07 lookup is key-scoped and fail-closed (records without
  key attribution are invisible to any key; cross-tenant access is tested and
  denied). *Fail closed*: unknown/invalid reference formats resolve to
  not-found, never to a best-effort match. No new secrets introduced; no
  secret values printed in evidence.
- **EC §4 (Product and scientific integrity)** — Provenance binding made
  explicit: `governance.audit_identity` binds the canonical input digest to
  the canon hash and matrix hash in force. No law meaning, threshold, or
  verdict semantics touched (OD-14 discipline: the only behavior changes trace
  to fixed defects — random ingestion-side identifiers, unscoped audit-object
  access, and the altloc default defect on the N8 branch).
- **EC §5 (Architecture and maintainability)** — Smallest complete change;
  API contract changes are strictly additive and optional
  (`audit_identity`, `legacy_audit_id` on `Governance`), so existing consumers
  stay compatible. Import graph and call sites inspected before removing the
  `uuid4()` identifier from the ingestion layer (all `Structure` consumers
  audited).
- **EC §6 (Testing and evidence)** — Regression tests written for every new
  behavior: identity triple, determinism, legacy resolution (owner + deny),
  fail-closed anonymous records, invalid refs, PDB legacy resolution, batch
  attribution. Exact command, branch commit and pass/fail/skip counts recorded
  in the preflight/evidence file; a test-order interference defect found in
  the first full-suite run was fixed in the test (per-test key fixture), not
  by weakening any assertion. Pre-existing base failures are disclosed, never
  hidden: the identical 9-failure set occurs on clean `ebc1fd7`
  (4× `ModuleNotFoundError: gemmi` — sandbox cannot install gemmi; 4×
  `test_b6_regression::TestUsageWidget`; 1× `test_runner_hardening`).
- **EC §7 (CI/CD and supply chain)** — Work developed on a feature branch;
  no push to `main`; no force-push (OD-15). CI (Toscanini Stability Gate)
  runs on the pushed branch SHA; the push is not claimed to be review or
  approval.
- **EC §11 (Change, release, and incident discipline)** — Each direction has a
  problem statement (audit finding), acceptance criteria (the direction text),
  tests, and a rollback plan (git revert of the single commit). No history
  rewrite anywhere.
- **EC §12 (Definition of done)** — Exact base and proposed SHAs recorded;
  changed-file manifests in commit messages; "not run" is never reported as
  "passed" (gemmi-limited tests disclosed above); nothing is claimed
  production-ready, nothing merged.
- **EC §13 (Prohibited shortcuts)** — No gate bypassed or weakened; no test
  skipped/xfail/deleted to obtain green CI; no secret committed; no law
  meaning, threshold, or ID altered; nothing merged or deployed that did not
  pass checks; no task claimed complete while material checks are unrun
  without stating the limitation.
- **EC §12 cross-reference (definition of done)** and **EC §13** apply to the
  follow-up branches as well (`feature/n8-altlocs-per-residue`,
  `hygiene/h1-remove-gemini`, `docs/owner-decisions`); each carries its own
  preflight note and evidence in its commits.
- **EC §14 (Required task artifacts)** — Preflight + baseline SHA, scope map,
  test results bound to the commit, changed-file manifest, impact note, and
  final checklist are delivered in `PREFLIGHT_OD07_2026-10-20.md` and the
  commit messages of this branch.

## 3. Domain constitution IDs (constitution.md, v6.1-Hardened)

- **OBS-153 (Avidity Quantification — ACTIVE)** — Read in full. NOT engaged by
  Round-12 directions 2–6: no bipodal/monopodal contact-density measurement,
  threshold, or report is touched. Recorded as out-of-scope rather than
  silently ignored; no change this session may alter its metric or standard.
- **CRIT-154 (Structural Admissibility — ACTIVE)** — Read in full. NOT engaged
  by Round-12 directions 2–6: no warhead-RMSD validation, clinical-benchmark
  comparison, or veto semantics touched. Same out-of-scope recording applies.

## 4. Governance decisions applied (register `docs/owner-decisions`)

- **OD-07** (implemented on this branch), **OD-08** (key scope), **OD-14**
  (defect-traceability discipline), **OD-15** (non-`main` branches only; no
  force-push; credential names only), **OD-16 → H-6/H-11**, **H-1** (separate
  hygiene branch from `main`), **H-8** (not issued — no action invented for
  it).

## 5. Session sequencing disclosure (EC §2 honesty)

The constitutions were fetched, hash-verified and read in full before the
first code edit. The preflight file was finalized immediately before the code
commit rather than before exploratory file reading; the exploration was
read-only (`git show`/`Read`), so no unrecorded mutation preceded the
preflight. This note is the reading-note artifact required by Round-12
Direction 1 and is committed as the first commit on this branch.
