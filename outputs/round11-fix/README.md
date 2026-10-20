# ROUND-11 FIX DELIVERY — README

> ## SUPERSESSION (2026-10-20, Round-12 Direction 5)
>
> **The register on this branch is superseded.** The partial register
> reconstruction carried on this branch — `docs/decisions/OWNER_DECISIONS.md`
> (rev 5) and the preserved `outputs/OWNER_DECISIONS.md` (rev 4) — is left
> UNEDITED here, as directed, and is **superseded by the authoritative
> reconciled decision register (45 rows: OD-01..OD-21, A-01..A-08, G-1..G-3,
> H-1..H-13)** on branch **`docs/owner-decisions` @ `0dcc5a1`** (pushed
> 2026-10-20; auditor-delegated authority; every row with status, authority,
> rationale, dependent work items, override path, and a separate
> implementation-status column). The "Standing governance state" section
> below ("All 21 ODs remain OPEN") describes the superseded rev-5 view and
> is likewise superseded — decision statuses now live in the authoritative
> register. The reconstruction gaps listed in "What is NOT in this delivery"
> were closed by the auditor's Round-12 reconciliation message and are
> resolved in that register, not here.

**Branch:** `evidence/round11-fix-2026-10-20`
**Base:** `origin/main` @ `2195a048e87dcec4aa63e8dc7981b3bb6d237ca0`,
carrying the Round-11 delivery tip
`7e3a3395be2d0719facf5c61ac532a0d347a49d3`
(`evidence/round9-r2-2026-10-20`) unchanged, plus the corrections below.
**Push scope:** this new branch only; no existing branch, tag, or ref
moved; nothing merged or deployed; PR #16 untouched.

**Purpose.** Execute the remediation directed by the Round-11 audit
(verdict NOT ACCEPTED, four failures). Two audit checks had already
PASSED (the Round-11 evidence branch exists at `7e3a339` based on
`main`; the lost tips `94b0f7c`/`e9985c6`/`df44aa6` are absent from
GitHub) and are not re-litigated here.

## What this delivery contains

| Path | What it is | Audit failure addressed |
|---|---|---|
| `docs/decisions/OWNER_DECISIONS.md` | **The canonical register, rev 5** — first time the register exists at its repository home on GitHub. Carries rev 4's Parts 1–5 forward verbatim in substance, adds Part 6 (audit reconciliation: passed checks, four-failure remediation table, reconstruction-gap register, SUPERSEDED ROWS table) and the A-07/A-08 rows | Failure 1 (structural part) |
| `outputs/round11-fix/RETRACTION_SECRETS_SCAN_CLAIM.md` | **Written retraction** of the Round-11 delivery README's clean-scan claim; documents the actual state (revoked Gemini literal at `system_diagnosis.py:8` inherited from `main`'s tree on every commit of the delivery branch) with the key value redacted | Failure 2 |
| `outputs/round11-fix/patches-n8-per-residue-default-v2/` | **N8 v2 patch series** (cover letter + 0001): flips `ALTLOC_POLICY` to `per_residue` as the ENGINE DEFAULT exactly as the audit directs; `per_atom` kept as explicit opt-in (original P0 item 6 semantics, byte-for-byte); selection logic unchanged from v1 | Failure 3 |
| `outputs/round11-fix/test-evidence/` | Raw test-run evidence for N8 v2: targeted altloc runs (19 passed / 3 failed, all 3 = `ModuleNotFoundError: gemmi`), full-suite runs v2 vs base `ebc1fd79` (508/12/5/1 vs 497/9/5/1), sorted failure-list diff (v2 = base + exactly the 3 gemmi-limited CIF tests — ZERO new failures), and the retraction fact evidence (ls-remote, merge-base, redacted line extracts) | Failures 2–3 |
| `outputs/round11-fix/SHA256SUMS` | SHA-256 of every file in this folder plus the register (except itself) | — |

## What is NOT in this delivery, and why (fail-closed record)

The Round-11 audit directs the register to contain **C2 completion, the
G-series rows, the H-9…H-13 rows, owner action A-08 (and A-07), and the
OD-04 amendment**. Their defining text lived in the Round-11 directive /
audit message. That message is **not present in this environment** — the
session that received it did not survive the sandbox reset, and a full
sweep of the workspace, all local git objects, the fetched remote refs,
the sealed r1/r2 archives, and the persisted tool-result records found
no copy. Per EC §1 (fail closed on ambiguity) and EC §13 (no
reconstruction of missing rules), the register (§6.3) records these as
reconstruction gaps by ID with status OPEN and **no invented content**.
They close when the owner re-issues the Round-11 directive/audit text;
no coder-side guessing is permitted.

The register's `outputs/OWNER_DECISIONS.md` (rev 4) and
`outputs/DELIVERY_README.md` from the prior delivery are preserved
byte-identical on this branch — never rewritten (EC §11).

## How to verify

```bash
git fetch origin evidence/round11-fix-2026-10-20
git checkout evidence/round11-fix-2026-10-20
sha256sum -c outputs/round11-fix/SHA256SUMS

# N8 v2 applies cleanly on the science base:
git worktree add --detach /tmp/n8check ebc1fd79
cd /tmp/n8check
git apply --check outputs/round11-fix/patches-n8-per-residue-default-v2/0001-*.patch

# Branch-tree secrets posture (the retracted claim, now stated correctly):
git show origin/evidence/round11-fix-2026-10-20:system_diagnosis.py | sed -n '8p'
```

Test evidence commands (exact): see `test-evidence/*.txt` headers.
Environment: Python 3.12.15, numpy 2.1.3, pytest 9.0.2,
`TOSCANINI_DATA_DIR` set to a writable directory, `gemmi` ABSENT
(PyPI unreachable in the delivery sandbox).

## Standing governance state

Hard stop in force (OD-15/OD-16 OPEN): no merges to `main`, no
deployments, no v25 run, no held-out evaluation, no threshold or cohort
changes, no golden-artifact changes. All 21 ODs remain OPEN. The only
new push is this evidence branch, owner-authorized for the audit
remediation.
