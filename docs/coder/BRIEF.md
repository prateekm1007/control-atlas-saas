# Coder Brief — Execution Rules for control-atlas-saas

**Status:** Canonical execution brief (OD-01 target layout)
**Audience:** any coder or automated agent executing work in this repository
**Authorities:** `constitution.md` (domain), `ENGINEERING_CONSTITUTION.md`
(engineering), `docs/decisions/OWNER_DECISIONS.md` (decisions). This brief
adds no new authority; it operationalizes theirs.

## 1. Before any work

1. Read `AGENTS.md`, both constitutions, and the decision register. Record
   the read in a reading note committed to your feature branch: cite every
   EC, OBS, and CRIT ID you enforce, the constitution blob hashes, and the
   fetch / complete-read / hash-verification / first-edit timestamps.
2. Pin your base: record the exact branch and commit SHA you start from.
   `main` moves; verify with a fresh fetch.
3. Write a preflight (EC §2 item 5): expected behavior, relevant risks,
   files expected to change, regression tests, rollout/rollback strategy,
   and acceptance criteria. Attach it to the task/PR/CI run.
4. Read every caller and consumer of what you will touch. Search for
   duplicate implementations and alternate entry points before modifying
   behavior (EC §5).

## 2. While working

- One focused patch per branch. Prefer the smallest complete change
  (EC §5); scope expansion needs separate justification.
- Never: weaken or delete a test to obtain green CI; convert an error,
  absent metric, timeout, or failed simulation into a normal verdict;
  silently relabel a metric; rename or redefine a law ID (EC §13).
- Synthetic/mock data is test-only. It must never receive an authoritative
  verdict or enter the production negative-knowledge store (EC §4).
- Credentials: use runtime environment variables; never write them to
  files, configs, logs, commit messages, or push URLs that persist.
  Secret-scan the working tree before every push.
- Identifiers: `OD-xx` one meaning forever; no `D-xx` in user-visible
  files; audit IDs follow OD-07 (deterministic full SHA-256, combined with
  the canon hash and matrix hash).

## 3. Verification duties (EC §6, EC §12)

- Reproduce bugs with failing regression tests before fixing when feasible.
- Run targeted tests, security tests, and the full suite on the exact
  candidate commit. Record: command, SHA, runtime versions, pass/fail
  counts. Reproduce the unmodified baseline in the same environment when
  pre-existing failures are suspected, and report both runs.
- Verify the hash locks where applicable (canon, matrix, law-semantic)
  and state that your change leaves them bit-identical — or stop and
  report, because law semantics drifted.
- Every patch states its **verdict impact** (or "none"), the commands run,
  and the output. "Verdict impact: none" requires an argument, not just
  the phrase.

## 4. Delivery

- Push feature branches; report branch name, SHA, and CI output. Do not
  push to `main`; do not force-push; do not merge or deploy (hard stop,
  OD-16, until the owner lifts it).
- Durable evidence: prefer repository-tracked, repo-relative artifacts
  under `docs/audit-evidence/<date>/…` with recorded input digests. No
  absolute local filesystem paths in evidence bundles.
- Final report: changes, evidence, unresolved risks, blocked checks, and
  what remains unverified. Nothing more, nothing hidden.

## 5. Standing blockers (do not cross without a decision)

- v25 calibration run; held-out evaluation; threshold or cohort changes;
  golden-artifact regeneration; merges; deploys; retire/redefine of any
  law; certificate issuance for non-PASS states (OD-06 restricts);
  publication of any unverified feature as production-ready.
