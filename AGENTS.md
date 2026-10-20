# AGENTS.md — Mandatory Entry Point

**Status:** Canonical governance entry point for `prateekm1007/control-atlas-saas`
**Established by:** OD-01 (governance files), auditor-decided 2026-10-20
**Read this file first.** Every coder, reviewer, and automated agent working in
this repository MUST read this entry point and the documents it links before
making any change.

## 1. Governing documents (read order)

1. **`constitution.md`** (repository root) — the domain/scientific policy.
   Currently states OBS-153 (avidity quantification, contact-density
   standard) and CRIT-154 (structural admissibility, warhead-RMSD and
   steric-clash veto). Both are ACTIVE.
2. **`ENGINEERING_CONSTITUTION.md`** (repository root) — the engineering,
   security, testing, and evidence policy (EC §1–§15 plus the Final Rule).
   This is the byte-exact canonical copy of the engineering constitution
   (blob `b6bb35e84956d84a8c9d9e65734c212baf3ed075`).
3. **`docs/coder/BRIEF.md`** — the execution brief: how work is actually
   done here (preflight, branches, evidence, gates).
4. **`docs/decisions/OWNER_DECISIONS.md`** — the controlling decision
   register (OD-01…OD-21 and A-01…A-06). Until the governance series is
   merged, the canonical copy lives on the `docs/owner-decisions` branch.
5. **`docs/P0_STATUS.md`** — implementation status, blockers, test
   evidence, and exact candidate SHAs for the P0 program.

## 2. Non-negotiable rules

- **Constitution-first.** No write, patch, commit, or test-driven code
  change before both constitutions have been read in full (EC §2). Record
  the read (paths, blob hashes, timestamps) in a reading note.
- **The hard stop is in force (OD-16).** No merges to `main`, no deploys,
  no v25 calibration run, no held-out evaluation, no threshold or cohort
  changes, and no golden-artifact regeneration without an explicit,
  recorded decision.
- **Feature branches only (EC §7).** All changes are developed on feature
  branches and reviewed through pull requests. Direct pushes to `main`
  are prohibited. Never force-push a published branch.
- **No verdict or law-semantics changes without a decision.** Thresholds,
  operators, scopes, law IDs and their meanings are frozen (OD-02, OD-03,
  OD-05, OD-14). The law canon hash, matrix hash, and law-semantic hash
  must remain bit-identical unless a recorded decision authorizes a
  reviewed change.
- **No secrets, ever (EC §3, EC §13).** No credential, key, token, or
  password may be committed, logged, echoed, or embedded in artifacts.
  A credential found in history is treated as compromised until its
  issuer confirms revocation.
- **Evidence or it did not happen (EC §6, EC §12).** Every claimed test
  result includes the exact command, commit SHA, environment, and counts.
  "Not run" is never reported as "passed."
- **One ID, one meaning.** `OD-xx` decision IDs have exactly the meaning
  fixed in the register and are never reused. Retired internal `D-xx`
  identifiers do not appear in user-visible files.

## 3. Change procedure (summary)

1. Read the governing documents above.
2. Pin the exact base SHA; never assume `main` has not moved.
3. Record a preflight: expected behavior, risks, files to change,
   regression tests, acceptance criteria, rollback plan.
4. Implement the smallest complete change on a feature branch.
5. Run targeted tests, security tests, and the full suite on the exact
   candidate; record commands and output.
6. Secret-scan the tree; verify the hash locks where applicable.
7. Push the feature branch and report its SHA and CI output.
8. Await review; do not merge or deploy without authorization.

## 4. Pointers

- CI: `.github/workflows/toscanini-gate.yml` (Stability Gate: secret scan,
  unit tests, canon/matrix/law-semantic hash locks).
- Deployment: `.github/workflows/deploy.yml` (main-only; blocked by the
  hard stop).
- Evidence conventions: `docs/audit-evidence/<date>/…`.
- Historical v23 cohort quarantine: see OD-17 and OD-10 in the register.
