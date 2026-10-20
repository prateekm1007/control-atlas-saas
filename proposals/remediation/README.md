# proposals/remediation — EC-25 remediation package (SUPERSEDED / HISTORICAL)

> **Status (Round-19 correction, 2026-10-20): every patch in this directory is a
> SUPERSEDED HISTORICAL DOCUMENT.** None of it describes the live repository state,
> and none of it is pending application. The remediation it proposed was shipped to
> `main` by the Round-17 release merge (`d292a38`, gate run 39411457483 = success).
> This file is retained as evidence of the reviewed remediation design and of the
> pre-release defects it corrected.

## Current live state (verified on `main` @ `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`)

| Workflow | Live trigger state (current, not historical) | Live blob |
|---|---|---|
| `.github/workflows/toscanini-gate.yml` | `push.branches: ["**"]` and `pull_request.branches: ["**"]` — nested branch names DO run the gate; installs the locked/dev requirements | `ff3d66d919666fe234220ed922d636b8d775f249` |
| `.github/workflows/deploy.yml` | `workflow_dispatch` ONLY; the single job checks out nothing, connects to no host, reads no secrets and performs NO deployment (inert by design) | `6b32942dd4c8fc73aace304ecb5595eaf6b83c9f` |

Neither defect this directory was authored to fix exists on `main` anymore. Earlier
versions of this README described the trigger defects and the push-triggered
deployment as if they were the live state; that description is corrected here
(Round-19 direction 3B) and must not be reintroduced.

## Historical contents (what each file was, and the exact defect it reproduced)

| File | Historical purpose | Historical source blob it reproduced at authoring time |
|---|---|---|
| `deploy.dispatch-only.yml` | Proposed replacement for the PRE-RELEASE `deploy.yml`, which triggered a deploy on every gate completion (`workflow_run`) and referenced EC2 hosts/credentials | reproduces pre-release `deploy.yml` blob `838627b460a374b2c789f73cdbcad3c8c38abc4d` (at `main` @ `2195a048`, superseded by blob `6b32942d…` from `d292a38`) |
| `toscanini-gate.branch-triggers.yml` | Proposed `on:`-block fix for the PRE-RELEASE gate, whose `push.branches: ["*"]` missed nested branch names and whose line 7 carried the malformed `pull_request.branches: ain]` fragment | reproduces pre-release `toscanini-gate.yml` blob `b86274ef3bd29f4eb4c534e71327a5d4a72fb812` (at `main` @ `2195a048`, superseded by blob `ff3d66d9…` from `d292a38`) |
| `fixtures/deploy_on_push_VIOLATION.yml` | Byte-exact copy of the PRE-RELEASE `deploy.yml` used by AT-14 to prove the lint flags it | byte-exact copy of blob `838627b460a374b2c789f73cdbcad3c8c38abc4d` (unchanged) |
| `fixtures/toscanini_gate_AIN_DEFECT.yml` | Byte-exact copy of the PRE-RELEASE gate used by AT-14 | copy of blob `b86274ef3bd29f4eb4c534e71327a5d4a72fb812`; NO LONGER byte-exact after the Round-18/19 provider-token removals on its pip-install line (the reproduced `on:`-block defects are unchanged and still machine-verified by AT-14) |

Both `pip install` lines in this directory previously carried the removed provider
SDK package token; it was removed from `toscanini-gate.branch-triggers.yml` and
`fixtures/toscanini_gate_AIN_DEFECT.yml` (Round-18 reviewed edit, carried to `main`
by the Round-19 repair). The token never appears in any live workflow.

## Policy status vs active enforcement (kept distinct, per Round-19 direction 3B)

- **Policy status:** EC-25 (no ungated deploy on push) is ACTIVE policy — recorded
  in `docs/directives/INDEX.md` and the decision register upon the Round-17 D-2
  merge (`b851e01`). A rule being in force is not proof that its automated checks
  run; the two are recorded separately.
- **Active enforcement (before the Round-19 repair):** none of the checks named in
  this package ran in CI — `proposals/ci/*.yml` are inert drafts (GitHub Actions
  reads only `.github/workflows/`), and the live gate ran no deploy-trigger lint,
  so the dispatch-only state of `deploy.yml` was policy-compliant but unenforced.
- **Active enforcement (from the Round-19 repair):** the live workflow
  `.github/workflows/governance-acceptance.yml` runs
  `scripts/ci/workflow_deploy_trigger_lint.py` against the actual current
  `.github/workflows/` tree on every push and pull request, so a reintroduced
  push/workflow_run-triggered deploy job now FAILS CI instead of passing unnoticed.

## Activation plan (SUPERSEDED — recorded as history, not a pending procedure)

The plan below is the ORIGINAL text, kept for the audit trail. Steps 1–2 were
executed by the Round-17 release merge (`d292a38`), not by this package; step 3's
lint activation is superseded by the Round-19 live governance workflow described
above.

1. Apply `toscanini-gate.branch-triggers.yml` (replace main's gate `on:` block) so
   nested branch pushes run the gate. — DONE via `d292a38` (superseded).
2. Apply `deploy.dispatch-only.yml` (replace main's `deploy.yml`). — DONE via
   `d292a38` (superseded; the shipped version is inert-by-design, stricter than
   this draft).
3. THEN copy `proposals/ci/no-deploy-on-push-lint.yml` into `.github/workflows/`.
   — SUPERSEDED by `.github/workflows/governance-acceptance.yml` (Round-19).
4. Re-run the lint against the remediated tree: expected PASS. — Executed by the
   live job on every run from the Round-19 repair forward.

## GitHub-settings prerequisites (still relevant only if deployment is ever reintroduced)

Per GitHub's documented deployment-protection model, the following must be
configured in repository settings and verified through the settings/API — a
workflow's `environment:` key alone protects nothing. These remain OPEN owner
actions (A-03/A-04 in the decision register) and are NOT claimed as configured:

- **Environment `production`** with required reviewers and a deployment-branch
  policy restricted to `main` (does not exist while no deployment runs).
- **Required checks** on `main` (A-04) — the exact required-check configuration
  remains unverified.
- **Direct pushes to `main`** disabled for normal development (A-04).
- **Secrets scoping**: deployment secrets must remain referenced only inside a
  protected deploy job if one is ever reintroduced.

## Evidence corrections recorded with this package (append-only; 2026-10-20)

1. `docs/evidence/ec16-ec26/EC25_LINT_ON_MAIN.txt` — correction appended: the
   earlier line implied the deploy trigger was `push: branches: ain]`. Actual:
   the pre-release deploy.yml triggered on `workflow_run`; the `ain]` fragment
   was in the pre-release `toscanini-gate.yml` line 7, blob `b86274ef…`
   (pull_request trigger).
2. The amendment document's EC-25 section is corrected in the same direction, and
   the `[ain]` claim is backed by exact tracked-file evidence (path, blob SHA,
   line) instead of a worklog citation.
