# P0 Hardening Evidence — 2026-10-19

Repository-published evidence set for the P0 hardening chain of
`prateekm1007/control-atlas-saas` (directive: "durable evidence — put the
results in GitHub"). Sanitized: no secret values, no patch bodies (hashes +
reproduction commands only), no customer data.

**Chain:** baseline `2195a04` → WP-1…WP-10 (`014b155`) → WP-12 `fd59206` → WP-13 `2c156b5` → WP-14 `c8faf2f` (this branch's integration tip includes the same commits).

## Contents

| File | What it is |
|---|---|
| `preflight.json` | Machine-readable preflight for WP-13/WP-14 + P0 containment (base SHAs, both policy-document SHA-256 hashes, caller maps, risks, acceptance tests, rollback) |
| `p0_credential_containment.json` | Full-history + tracked-tree + worktree gitleaks scans (tool 8.30.1, redacted); issuer-side rotation recorded **BLOCKED** (SEC-018 pending provider) |
| `compliance_matrix.json` | 71-requirement retrospective audit against SaaS Engineering Constitution v1.0.0 (satisfied 36 / partial 23 / not_satisfied 7 / N.A. 2 / not_verified 3; gaps prioritized) |
| `compliance_report.md` | Human-readable rendering of the matrix |
| `security_findings.json` | Register rev 1.1 — 22 findings: 15 fixed_in_pr, 4 open, 1 corrected_false_positive (SEC-016), 1 pending_provider (SEC-018), 1 unverified (SEC-020) |
| `wp_scope_manifest.json` | WP-1…WP-14 scope + per-WP `git diff --name-status` manifest |
| `patch_hashes.json` | Per-WP diff SHA-256 fingerprints + clean-checkout reproduction commands (bodies intentionally not committed) |
| `test_summary.json` | Exact pytest counts bound to each WP head (244 baseline → 383 at WP-14; 0 failed/0 errors from WP-7 on), CI run references, protected-hash verification log |
| `MANIFEST.sha256` | SHA-256 of every file in this directory (self-exclusive) |

## Clean-checkout verification

```bash
git clone https://github.com/prateekm1007/control-atlas-saas
cd control-atlas-saas
git checkout <commit containing this directory>   # integration candidate head

# 1. every evidence hash must validate:
cd docs/audit-evidence/2026-10-19 && sha256sum -c MANIFEST.sha256

# 2. per-WP patch hashes must reproduce from the commit chain:
#    (ranges in wp_scope_manifest.json; e.g. WP-13)
git diff fd59206..2c156b5 | sha256sum      # == patch_hashes.json WP-13 entry

# 3. the full chain reproduces the audited tree (first 10 WPs):
#    git am --keep-cr combined-chain-2195a04..014b155.mbox  (mbox in the local
#    evidence package; resulting tree b2019910... == 014b155 exactly)

# 4. re-run the suite at this head:
cd backend && python -m pytest tests/ -q    # expect 383 passed / 5 skipped
```

## Policy documents in force for this work

- Domain: `constitution.md` @ baseline — SHA-256 `cfdd5ba3babda31513b55894a49597fbfb18c2e870535a068c8b0ea423c53ce8`
- Engineering: `ENGINEERING_CONSTITUTION.md` @ `4993cb3` (PR #11, proposed) — SHA-256 `e0f256f78f2c1256228c4286311cba15800dcd5a6c8005583508ee66d27983f5`

## Disposition checklist at publication time

- **Completed (evidenced):** P0 scans archived; SEC-015 fixed (WP-13, fail-closed + 12 tests); deployment gating policy implemented + asserted (WP-14, 16 assertions); suite 383/5/0/0 at chain tip; protected hashes verified at every head; register rev 1.1.
- **Blocked (owner/provider):** SEC-018 issuer-side key revocation; PR merges (#1–#14, then this integration PR); branch-protection/required-checks enforcement (configured/verified separately via API); actual EC2 deployment with secrets + health record.
- **Unverified by design:** live deployment behavior (no secrets in sandbox — policy proven via assertions); Docker builds (SEC-020); LLM narrative paths; accessibility; SLOs.
- **Open findings:** SEC-015 ✅ fixed · SEC-016 ✅ corrected (false positive) · SEC-017/019/021 open (gpu_worker test debt, CI scope, hermeticity) · SEC-018 blocked on provider · SEC-020 unverified.
