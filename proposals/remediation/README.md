# proposals/remediation — EC-25 remediation package (INERT until the approved merge path)

Auditor item 3.7: activating the no-deploy-on-push lint alone would introduce a check that
FAILS on the existing `main` baseline (its `deploy.yml` triggers a deploy on every push).
The full remediation is therefore prepared here as an inert, reviewed patch set. Nothing in
this directory is a live workflow; nothing here merges or deploys anything.

## STATUS (2026-10-20, Round-18 correction — read before citing anything below)

**The gate/defect half of this package is HISTORICAL AND SUPERSEDED.** It was written
against the pre-release mainline `main @ 2195a048e87dcec4aa63e8dc7981b3bb6d237ca0`
(gate blob `b86274ef3bd29f4eb4c534e71327a5d4a72fb812`, deploy blob
`838627b460a374b2c789f73cdbcad3c8c38abc4d`). The released mainline
(`d292a38…`, merged as `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`) already repaired
both defects, verified on the live tree:

- `.github/workflows/toscanini-gate.yml` on `main @ f535a71` carries
  `push.branches: ["**"]` AND `pull_request.branches: ["**"]` — the `["*"]` / `ain]`
  defects below no longer exist on main. The two files
  `toscanini-gate.branch-triggers.yml` and `fixtures/toscanini_gate_AIN_DEFECT.yml`
  are retained ONLY as historical evidence of the pre-release defect (exact source
  blobs cited above); they are not an apply target anymore.
- `.github/workflows/deploy.yml` on `main @ f535a71` is `workflow_dispatch`-only and
  intentionally performs NO deployment (inert by design, auditor Round-17 direction 2a).
  `fixtures/deploy_on_push_VIOLATION.yml` remains a byte-exact copy of the OLD
  push-triggered deploy (blob `838627b4…`) kept so AT-14 can prove the lint flags it.

**What is actually in force vs actually running on `main` today (do not conflate):**

| Thing | Status | Evidence |
|---|---|---|
| EC-25 (no ungated deploy on push) | In force as POLICY (auditor D-5, 2026-10-20) | register rev 7+, `docs/decisions/OWNER_DECISIONS.md` |
| `deploy.yml` on main | Dispatch-only, inert — no deploy logic, no push/workflow_run trigger | live blob on `main @ f535a71` |
| Toscanini Stability Gate | The ONLY live automated enforcement running on every push/PR (`["**"]`): secret scan, unit tests, canon/matrix hash locks | live blob on `main @ f535a71` |
| `no-deploy-on-push-lint.yml` and the other `proposals/ci/*.yml` drafts | INERT — under `proposals/`, GitHub Actions never reads them | path + `proposals/ci/README.md` move-map |
| The deploy-trigger lint (`scripts/ci/workflow_deploy_trigger_lint.py`) | Script exists; wired into live CI only by `.github/workflows/acceptance-gate.yml` (Round-18). Before Round-18 it was NOT invoked by any live workflow | live workflow on the Round-18 branch |

Nothing in this directory merges or deploys anything, and no file here is a live
enforcement point by virtue of existing. An `ACTIVE`/in-force policy status never
implies that the automation named in these drafts is running (EC-26 discipline:
policy status and enforcement status are recorded separately in
`docs/directives/INDEX.md`).

## Contents

| File | Purpose |
|---|---|
| `deploy.dispatch-only.yml` | Proposed replacement for `.github/workflows/deploy.yml`: `workflow_dispatch`-only, restricted to the approved ref (`main`), exact-SHA validation against the remote head, required-checks-green gate on that exact SHA, `environment: production`, `concurrency` group, deploy secrets referenced only inside the protected job, post-deploy health-signal step (a successful deploy command is not health, EC §8). |
| `toscanini-gate.branch-triggers.yml` | Proposed `on:`-block fix for `.github/workflows/toscanini-gate.yml`: `push.branches: ["**"]` (matches nested branch names; the current `["*"]` does not) and `pull_request.branches: [main]` (repairs the malformed `ain]` fragment). Job bodies byte-identical to main's. |
| `fixtures/deploy_on_push_VIOLATION.yml` | Byte-exact copy of main's current `deploy.yml` (blob `838627b460a374b2c789f73cdbcad3c8c38abc4d`), used by AT-14 to prove the lint FLAGS it. |
| `fixtures/toscanini_gate_AIN_DEFECT.yml` | Byte-exact copy of main's current `toscanini-gate.yml` (blob `b86274ef3bd29f4eb4c534e71327a5d4a72fb812`, line 7 `branches: ain]`), used by AT-14. |

## Activation plan (owner-approved path only, per the move-map)

1. Apply `toscanini-gate.branch-triggers.yml` (replace main's gate `on:` block) so nested
   branch pushes run the gate.
2. Apply `deploy.dispatch-only.yml` (replace main's `deploy.yml`).
3. THEN copy `proposals/ci/no-deploy-on-push-lint.yml` into `.github/workflows/`. The lint
   goes active only after steps 1–2, or it would fail on the un-remediated baseline.
4. Re-run the lint against the remediated tree: expected PASS.

## GitHub-settings prerequisites (verify, never presume from YAML)

Per GitHub's documented deployment-protection model, the following must be configured in
repository settings and verified through the settings/API — a workflow's `environment:` key
alone protects nothing:

- **Environment `production`** must EXIST with: required reviewers (deployment approval) and
  a deployment-branch policy restricted to `main`.
- **Required checks** ("Unit Tests + Hash Lock", "Secret Scan (gitleaks)") configured on
  `main` (A-04) — the dispatch gate step also verifies them at run time, but the
  branch-protection configuration is the durable control.
- **Direct pushes to `main`** disabled for normal development (A-04).
- **Secrets scoping**: `EC2_HOST` / `EC2_USER` / `EC2_SSH_KEY` must remain referenced only in
  the protected deploy job (they are environment/job-scoped references in the proposed file).

## Evidence corrections recorded with this package (append-only; 2026-10-20)

1. `docs/evidence/ec16-ec26/EC25_LINT_ON_MAIN.txt` — correction appended: the earlier line
   implied the deploy trigger was `push: branches: ain]`. Actual: deploy.yml's push trigger
   is `branches: [main]` (no typo); the `ain]` fragment is in `toscanini-gate.yml` line 7,
   blob `b86274ef…` (pull_request trigger).
2. The amendment document's EC-25 section is corrected in the same direction, and the
   `[ain]` claim is now backed by exact tracked-file evidence (path, blob SHA, line) instead
   of a worklog citation.
