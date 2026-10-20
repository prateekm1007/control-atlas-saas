# proposals/ci — INERT enforcement scaffolding for EC-16…EC-26

**Nothing in this directory runs anywhere.** GitHub Actions only reads `.github/workflows/`; the PR template only takes effect at `.github/PULL_REQUEST_TEMPLATE.md`. Everything here is a draft that activates only when the owner approves the corresponding amendment and the move-map below is executed on the approved path.

## Move-map (on owner approval, per batch)

| Adoption batch | Amendment | Activate by | Becomes |
|---|---|---|---|
| 1 | EC-17 | (no CI file) | Practice + audit checklist; push log (`docs/evidence/push-log/PUSH_LOG.md`) already live on this branch |
| 1 | EC-25 | `cp no-deploy-on-push-lint.yml ../../.github/workflows/` | Active lint: fails any PR adding a deploy job on `push`/`workflow_run`. **Note:** `main`'s existing `deploy.yml` violates EC-25 and is remediated on the approved merge path (recorded in the main amendment doc) |
| 2 | EC-22 | `cp status-vocabulary-lint.yml ../../.github/workflows/` | Active report lint |
| 2 | EC-23 | `cp gitleaks-all-branches.yml ../../.github/workflows/` | Active gitleaks on every branch + weekly `main` full-history scan |
| 3 | EC-16 + EC-19 + EC-20 | `cp verify-claimed-shas.yml ../../.github/workflows/` and `cp PULL_REQUEST_TEMPLATE.md ../../.github/PULL_REQUEST_TEMPLATE.md` | PR-body field checks (Directive-ID, Claimed-SHAs, Base-SHA, Acceptance-Tests) + remote-SHA verification + acceptance-tests-in-diff |
| 3 | EC-18 | `cp evidence-sequence-check.yml ../../.github/workflows/` | Active sequence checker over `docs/evidence/**` logs |
| 3 | EC-21 | `cp check-declared-base.yml ../../.github/workflows/` | Active declared-base ancestry check |
| 3 | EC-26 | `cp reading-note-map-check.yml ../../.github/workflows/` | Active reading-note mapping check |

## Checker scripts (shared)

All jobs call `scripts/ci/*.py`. The scripts are stdlib + PyYAML only, self-tested by `scripts/ci/tests/run_all.py` (AT-01…AT-07, all passing at packaging time — see the receipt).

## Design constraints honored by every workflow draft

1. **Trigger-safe:** every draft uses `branches: ["**"]` for push triggers — the documented `main` gate defect (`["*"]` + `[ain]` typo, worklog Task 35) must not be reproduced.
2. **EC-25-clean:** no draft contains a deploy job; they are checks. The lint itself (`no-deploy-on-push-lint.yml`) runs on push by design — a check is not a deploy.
3. **Least tokens:** jobs use `${{ github.token }}`; no PAT in CI.
4. **Fail-closed:** unknown/missing fields fail the check (the EC-22 vocabulary and EC-16 SHA checks reject empty input rather than passing it).
