# proposals/remediation — EC-25 remediation package (INERT until the approved merge path)

Auditor item 3.7: activating the no-deploy-on-push lint alone would introduce a check that
FAILS on the existing `main` baseline (its `deploy.yml` triggers a deploy on every push).
The full remediation is therefore prepared here as an inert, reviewed patch set. Nothing in
this directory is a live workflow; nothing here merges or deploys anything.

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
