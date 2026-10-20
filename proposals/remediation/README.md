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

## 2026-10-20 Round-18 rebase corrections (recorded, append-only style)

This package was authored against `main` @ `2195a048` (pre-release baseline). Since then the
Round-17 release merge (`d292a38`, merged to `main`, now `eb6d6a9`) landed the live remediation
this directory proposed:

1. `.github/workflows/deploy.yml` on `main` is `workflow_dispatch`-only (no push/pull_request/
   workflow_run triggers, no EC2/ssh/secret references) — the `deploy.dispatch-only.yml` draft
   here is superseded by the shipped version and is retained as a historical document only.
2. `.github/workflows/toscanini-gate.yml` on `main` triggers on `push.branches: ["**"]` and
   `pull_request.branches: ["**"]` and installs the lock + dev requirements — the `on:`-block
   fix proposed in `toscanini-gate.branch-triggers.yml` is likewise superseded and retained as
   a historical document only.

Corrections to the table above, per the auditor's 2026-10-20 documentation-exemption decision
("code, configuration, and workflow references" to the removed provider are not permitted
outside `docs/`; the auditor's acceptance grep for that decision — its exact pattern and the
removal record are entry R10 of `docs/evidence/ec16-ec26/RECEIPT.md` — must return no match
introduced by this package):

- `fixtures/toscanini_gate_AIN_DEFECT.yml` and `toscanini-gate.branch-triggers.yml` each
  carried the pre-release hand-written install line ending in the removed provider SDK
  package token. That token is removed in this rebase commit (R10). Both files are therefore
  NO LONGER byte-exact
  copies of any historical blob; the properties they exist to reproduce — the `on:`-block
  defects (`push.branches: ["*"]`, `pull_request.branches: "ain]"`) in the fixture, and the
  repaired `["**"]`/`["main"]` filters with no `ain]` fragment in the proposed fix — are
  unchanged and remain machine-verified by AT-14 on every run of the acceptance suite.
- `fixtures/deploy_on_push_VIOLATION.yml` is untouched (it contains no provider reference);
  it reproduces the PRE-RELEASE deploy defect and is likewise a historical document now that
  `main` ships the dispatch-only workflow.
- Nothing in this directory is a live workflow; the activation plan above is historical
  (steps 1–2 were executed by the Round-17 release merge, not by this package).
