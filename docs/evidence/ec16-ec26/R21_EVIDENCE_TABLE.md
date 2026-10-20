# Round-21 evidence table (GENERATED FROM LIVE EVIDENCE)

- **Generated:** 2026-10-20T18:14:03Z by `scripts/ci/generate_evidence_table.py` (committed; re-runnable).
- **Sources:** anonymous `git ls-remote` against `https://github.com/prateekm1007/control-atlas-saas.git` at the generation timestamp; the GitHub Actions API (`/repos/prateekm1007/control-atlas-saas`) for every run citation found in the push log; the push log itself for the entry inventory. Nothing in this table is hand-copied.
- **Scope note:** the exact final candidate SHA of the Round-21 repair is the commit carrying the terminal `P-060-INTENT` (this table is part of that commit, so it cannot embed the commit's own SHA); the candidate is bound by the push-log terminal intent and PR #30's body after the push. Run rows below are the verified pre-candidate history.

## 1. Live branch tips (anonymous ls-remote at generation time)

| Branch | Live tip SHA |
|---|---|
| `main` | `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` |
| `repair/round19-mainline-acceptance` | `62fab8a5ad0c07c71a508cb34b93e1a4af24e73e` |
| `activation/ec16-ec19-ec20-r19` | `bd8d86a3a93c1717921587f9585f2b97baf22836` |
| `activation/ec21-r19` | `32d272d3243f62224248c9d29eb5dcb0c8a20e5a` |
| `activation/ec26-r19` | `151b8c12d6e5286cc3153461677b193b6d57495d` |
| `activation/ec18-r19` | `97067ad0320968bfb73f7e045408fcb60ccd5c7d` |
| `activation/ec23-r19` | `fa854f8f3b7fcff9c83db5570784970f42df9742` |
| `activation/ec24-r19` | `28ee983e0f98723c096ce4b13892f5c16bfbb22e` |
| `docs/owner-decisions-r19` | `df98f5e3671055d67961d5b9d14495c5b037ce19` |
| `docs/owner-decisions` | `afd445cc3c8eafe7950a04df95e4f16b303b194a` |
| `governance/canonical-docs` | `d51ad7c7fcbe521c116b8f32be4904ae6a57ab41` |
| `proposal/ec16-ec26-enforcement` | `c4a0fcc3357930b7eba61e056489c804b16d2277` |
| `proposal/ec16-ec26-enforcement-rebased` | `3745cd2af2e0b5b909cb4c3b49c515afdec4cce8` |

## 2. Push-log entry inventory and remote-ref classification

| Entry | Verb | Branch | Recorded remote_sha | Live tip | Classification |
|---|---|---|---|---|---|
| `P-001` | push | `evidence/round9-r2-2026-10-20` | `-`… | `-`… | legacy push record |
| `P-002` | push | `evidence/round11-fix-2026-10-20` | `-`… | `-`… | legacy push record |
| `P-003` | push | `feature/od07-audit-id` | `-`… | `-`… | legacy push record |
| `P-004` | push | `feature/n8-altlocs-per-residue` | `-`… | `-`… | legacy push record |
| `P-005` | push | `hygiene/h1-remove-gemini` | `-`… | `-`… | legacy push record |
| `P-006` | push | `docs/owner-decisions` | `-`… | `afd445cc3c8e`… | legacy push record |
| `P-007` | push | `evidence/round11-fix-2026-10-20` | `-`… | `-`… | legacy push record |
| `P-008` | push | `feature/od07-audit-id` | `-`… | `-`… | legacy push record |
| `P-009` | push | `hygiene/h1-remove-gemini` | `-`… | `-`… | legacy push record |
| `P-010` | push | `docs/owner-decisions` | `-`… | `afd445cc3c8e`… | legacy push record |
| `P-011` | push | `proposal/ec16-ec26-enforcement` | `-`… | `c4a0fcc33579`… | legacy push record |
| `P-011-CONFIRM` | push-confirm | `proposal/ec16-ec26-enforcement` | `17eb88184c9f`… | `c4a0fcc33579`… | HISTORICAL (valid for the push it describes) |
| `P-012-INTENT` | push-intent | `docs/owner-decisions` | `-`… | `afd445cc3c8e`… | INTENT (terminal pending disclosed by the checker) |
| `P-013-INTENT` | push-intent | `proposal/ec16-ec26-enforcement` | `-`… | `c4a0fcc33579`… | INTENT (terminal pending disclosed by the checker) |
| `P-014-INTENT` | push-intent | `governance/canonical-docs` | `-`… | `d51ad7c7fcbe`… | INTENT (terminal pending disclosed by the checker) |
| `P-012-CONFIRM` | push-confirm | `docs/owner-decisions` | `78b147d55660`… | `afd445cc3c8e`… | HISTORICAL (valid for the push it describes) |
| `P-013-CONFIRM` | push-confirm | `proposal/ec16-ec26-enforcement` | `dab33d946fe3`… | `c4a0fcc33579`… | HISTORICAL (valid for the push it describes) |
| `P-014-CONFIRM` | push-confirm | `governance/canonical-docs` | `b12f30e95255`… | `d51ad7c7fcbe`… | HISTORICAL (valid for the push it describes) |
| `P-015-INTENT` | push-intent | `proposal/ec16-ec26-enforcement-rebased` | `-`… | `3745cd2af2e0`… | INTENT (terminal pending disclosed by the checker) |
| `P-015-CONFIRM` | push-confirm | `proposal/ec16-ec26-enforcement-rebased` | `0f1e9a8856ba`… | `3745cd2af2e0`… | HISTORICAL (valid for the push it describes) |
| `P-017-INTENT` | push-intent | `repair/round19-mainline-acceptance` | `-`… | `62fab8a5ad0c`… | INTENT (terminal pending disclosed by the checker) |
| `P-017-CONFIRM` | push-confirm | `repair/round19-mainline-acceptance` | `2dcbbd781a2c`… | `62fab8a5ad0c`… | HISTORICAL (valid for the push it describes) |
| `P-050` | push-reconcile | `repair/round19-mainline-acceptance` | `f44bf4d93f8d`… | `62fab8a5ad0c`… | HISTORICAL (valid for the push it describes) |
| `P-051` | push-reconcile | `activation/ec16-ec19-ec20-r19` | `bd8d86a3a93c`… | `bd8d86a3a93c`… | CURRENT |
| `P-052` | push-reconcile | `activation/ec21-r19` | `32d272d3243f`… | `32d272d3243f`… | CURRENT |
| `P-053` | push-reconcile | `activation/ec26-r19` | `151b8c12d6e5`… | `151b8c12d6e5`… | CURRENT |
| `P-054` | push-reconcile | `activation/ec18-r19` | `97067ad03209`… | `97067ad03209`… | CURRENT |
| `P-055` | push-reconcile | `activation/ec23-r19` | `fa854f8f3b7f`… | `fa854f8f3b7f`… | CURRENT |
| `P-056` | push-reconcile | `activation/ec24-r19` | `28ee983e0f98`… | `28ee983e0f98`… | CURRENT |
| `P-057` | push-reconcile | `docs/owner-decisions-r19` | `df98f5e36710`… | `df98f5e36710`… | CURRENT |
| `P-058-INTENT` | push-intent | `repair/round19-mainline-acceptance` | `-`… | `62fab8a5ad0c`… | INTENT (terminal pending disclosed by the checker) |
| `P-058-CONFIRM` | push-confirm | `repair/round19-mainline-acceptance` | `c1680eaf3d33`… | `62fab8a5ad0c`… | HISTORICAL (valid for the push it describes) |
| `P-059-INTENT` | push-intent | `repair/round19-mainline-acceptance` | `-`… | `62fab8a5ad0c`… | INTENT (terminal pending disclosed by the checker) |
| `P-059-CONFIRM` | push-confirm | `repair/round19-mainline-acceptance` | `62fab8a5ad0c`… | `62fab8a5ad0c`… | CURRENT |
| `P-060-INTENT` | push-intent | `repair/round19-mainline-acceptance` | `-`… | `62fab8a5ad0c`… | INTENT (terminal pending disclosed by the checker) |
| `P-061` | merge-intent | `repair/round19-mainline-acceptance` | `-`… | `62fab8a5ad0c`… | MERGE-INTENT |
| `P-062` | evidence-back | `proposal/ec16-ec26-enforcement` | `17eb88184c9f`… | `c4a0fcc33579`… | EVIDENCE-BACK (reviewed historical support) |
| `P-063` | evidence-back | `docs/owner-decisions` | `78b147d55660`… | `afd445cc3c8e`… | EVIDENCE-BACK (reviewed historical support) |
| `P-064` | evidence-back | `proposal/ec16-ec26-enforcement` | `dab33d946fe3`… | `c4a0fcc33579`… | EVIDENCE-BACK (reviewed historical support) |
| `P-065` | evidence-back | `governance/canonical-docs` | `b12f30e95255`… | `d51ad7c7fcbe`… | EVIDENCE-BACK (reviewed historical support) |

## 3. Cited Actions runs (GitHub API at generation time)

| Run ID | Workflow | Event | Head branch | Head SHA | Conclusion |
|---|---|---|---|---|---|
| 39408223554 | Toscanini Stability Gate | push | docs/owner-decisions | `78b147d55660…` | failure |
| 39408297643 | Toscanini Stability Gate | push | governance/canonical-docs | `b12f30e95255…` | failure |
| 39438362505 | Toscanini Stability Gate | push | repair/round19-mainline-acceptance | `2dcbbd781a2c…` | success |
| 39438362567 | Governance Acceptance | push | repair/round19-mainline-acceptance | `2dcbbd781a2c…` | failure |
| 39438788810 | Toscanini Stability Gate | pull_request | repair/round19-mainline-acceptance | `2dcbbd781a2c…` | success |
| 39438788822 | Governance Acceptance | pull_request | repair/round19-mainline-acceptance | `2dcbbd781a2c…` | failure |
| 39443515998 | Toscanini Stability Gate | pull_request | repair/round19-mainline-acceptance | `f44bf4d93f8d…` | success |
| 39443516042 | Governance Acceptance | pull_request | repair/round19-mainline-acceptance | `f44bf4d93f8d…` | success |
| 39443519540 | Governance Acceptance | push | activation/ec16-ec19-ec20-r19 | `bd8d86a3a93c…` | success |
| 39443538818 | Governance Acceptance | push | activation/ec21-r19 | `32d272d3243f…` | success |
| 39443554590 | Governance Acceptance | push | activation/ec26-r19 | `151b8c12d6e5…` | success |
| 39443566175 | Governance Acceptance | push | activation/ec18-r19 | `97067ad03209…` | success |
| 39443578786 | Governance Acceptance | push | activation/ec23-r19 | `fa854f8f3b7f…` | success |
| 39443597533 | Toscanini Stability Gate | push | activation/ec24-r19 | `28ee983e0f98…` | success |
| 39443611965 | Governance Acceptance | push | docs/owner-decisions-r19 | `df98f5e36710…` | success |
| 39457834846 | Toscanini Stability Gate | push | repair/round19-mainline-acceptance | `c1680eaf3d33…` | success |
| 39457834926 | Governance Acceptance | push | repair/round19-mainline-acceptance | `c1680eaf3d33…` | success |
| 39457843044 | Toscanini Stability Gate | pull_request | repair/round19-mainline-acceptance | `c1680eaf3d33…` | success |
| 39457843064 | Governance Acceptance | pull_request | repair/round19-mainline-acceptance | `c1680eaf3d33…` | success |
| 39458583993 | Governance Acceptance | push | repair/round19-mainline-acceptance | `62fab8a5ad0c…` | success |
| 39458583995 | Toscanini Stability Gate | push | repair/round19-mainline-acceptance | `62fab8a5ad0c…` | success |
| 39458592968 | Toscanini Stability Gate | pull_request | repair/round19-mainline-acceptance | `62fab8a5ad0c…` | success |
| 39458593293 | Governance Acceptance | pull_request | repair/round19-mainline-acceptance | `62fab8a5ad0c…` | success |
| 39459283515 | Toscanini Stability Gate | pull_request | repair/round19-mainline-acceptance | `62fab8a5ad0c…` | success |
| 39459283617 | Governance Acceptance | pull_request | repair/round19-mainline-acceptance | `62fab8a5ad0c…` | success |

## 4. Remaining pending state (disclosed, never smoothed over)

- `P-060-INTENT` is the terminal unpaired intent for `repair/round19-mainline-acceptance`: the Round-21 candidate push it declares is disclosed PENDING until a next legitimate push delivers `P-060-CONFIRM` (labeled PENDING, never COMPLETE; no push is manufactured merely to eliminate it, per the Round-21 adjudication).
- `P-061` (merge-intent) is pre-declared and PENDING until the reviewed merge of PR #30 executes through the protected path; the post-merge `main` push event is declared by it once the checker verifies the merge-commit parents, the merged PR, and the approving non-author review.
- PR #30 remains UNAPPROVED: merge awaits the auditor's review, the full gate on the exact final candidate, the demonstrated post-merge transition, and the required independent approval.
- Standing OPEN items (not closed by this round): independent r2 evidence-archive verification; portable 40-row cohort metadata table; A-02 issuer-side revocation; EC2 deployment secrets (deployment stays blocked).
- The PAT-hygiene sweep is a LOCAL coder claim; GitHub CI does not independently verify it. No token value appears in this table, the repository, or any worklog.
