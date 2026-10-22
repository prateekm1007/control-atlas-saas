# OD-scope iteration 3 — constitution-first reading note

Recorded: 2026-10-22T12:41:30Z
Repository: `https://github.com/prateekm1007/control-atlas-saas.git`
Branch: `review/od-scope-v25`
Base: `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
HEAD: `9412e943014df6804e34c670d9cc9496e1f34c9c`
Worktree: clean before the iteration-3 evidence files were created.

## Controlling directions

The owner's 2026-10-22 coder directions control this follow-up iteration. They
authorize a focused follow-up patch on the existing review branch, preserve the
delivery commit and its evidence, and require the engine, calibration-tool, semantic-lock,
and push-log defects to be corrected before any merge or deployment. They do not
authorize a v25 calibration run, held-out evaluation, threshold/cohort change, or
golden-artifact regeneration.

The owner also ruled that `LAW_CANON` is the runtime authority for `LAW-100`: the
canonical metric is the percentage of bonds exceeding 4 sigma, with `<= 5.0%`.
The unregistered experimental RMSZ comparison at `8.0` must be removed. This is
recorded in `owner-ruling-2026-10-22.md` without creating or reusing an OD ID.

## Read record

All six governing documents were read in full in this iteration.

| Document | Bytes | SHA-256 | Git blob |
|---|---:|---|---|
| `AGENTS.md` | 4021 | `5e320939a11135a0d8aed5caa231df6aeaf633d724fcb4dba97b5909e51420b4` | `85767237399ffb30f3e45894d9dfbe6f4240584f` |
| `constitution.md` | 622 | `cfdd5ba3babda31513b55894a49597fbfb18c2e870535a068c8b0ea423c53ce8` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` |
| `ENGINEERING_CONSTITUTION.md` | 20572 | `e0f256f78f2c1256228c4286311cba15800dcd5a6c8005583508ee66d27983f5` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` |
| `docs/coder/BRIEF.md` | 3671 | `8f07edfbba9d3e57de623e242125e942369bc5d60df5805343b26b9048c52e07` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` |
| `docs/decisions/OWNER_DECISIONS.md` | 72329 | `59ccd5e84176a8a610c8c234cbb95fec4734478cc1f35668599ad2e3aef3172a` | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` |
| `docs/P0_STATUS.md` | 35347 | `3f129dd18b9e5cf8888b763158cda42de270a8429a48ed7cc9a3aa675e810e0f` | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` |

- Controlled complete-read output window: `2026-10-22T12:33:57.967Z` to
  `2026-10-22T12:33:59.031Z`.
- Hash-verification timestamp: `2026-10-22T12:33:59.270Z`.
- Remote-ref verification timestamp: `2026-10-22T12:40:43Z`
  (`git ls-remote --heads origin`, 17 heads; output saved at `remote-heads.txt`).
- First edit in this iteration: `2026-10-22T12:41:30Z` (this reading note).

The iteration-2 reading note and preflight in
`docs/audit-evidence/2026-10-21/od-scope/iteration-2/` are preserved unchanged
as historical evidence. This note is a new record, not an edit of them.

## Active decision IDs

The following recorded requirements remain binding for this iteration:

- `OD-01` — canonical governance layout.
- `OD-02` — coverage-first, fail-closed verdict precedence.
- `OD-03` — pLDDT advisory only for LAW-130/LAW-200; retire LAW-155 from
  enforced semantics; no new pLDDT threshold.
- `OD-04` — apply OD-03 to disordered proteins; honest
  `NOT_MEASURED`/`MEASUREMENT_FAILED`/`INDETERMINATE` states.
- `OD-05` — retire unsafe laws with IDs retained; no implicit PASS.
- `OD-09` — freeze unapproved scientific changes, calibration, threshold
  tuning, cohort edits, and golden artifacts.
- `OD-10` — no cohort/threshold criteria from contaminated metadata.
- `OD-13` — pin AFDB model/version, host/redirect policy, format/size,
  coordinate SHA-256, and run provenance; no moving `latestVersion`.
- `OD-14` — V4 and the historical verdict table remain diagnostic only.
- `OD-16` — hard stop: no merge, deploy, v25 run, held-out evaluation,
  threshold/cohort change, or golden-artifact regeneration.
- `OD-17` — preserve v23 dataset/results byte-for-byte; corrected metadata
  in a new version.
- `OD-18` — keep `1I0T` and `1DDL` in a separate audit set.
- `OD-19` — no quota assumptions; prespecify feasibility and
  family/sequence-disjoint splits after reconciliation.
- `OD-20` — keep NMR separate from resolution-based strata; no LAW-182
  policy change under this ID.
- `A-01` — accepted risk; the owner rotates the PAT at the owner's
  pleasure; no credential operations, no token value in any artifact.
- `A-04` — branch-protection configuration remains unverified/open.

## Governance discrepancy: canonical register path

`AGENTS.md` says the controlling register is
`docs/decisions/OWNER_DECISIONS.md` and that its canonical copy lives on the
`docs/owner-decisions` branch until the governance series is merged.

Live verification at `2026-10-22T12:40:43Z`:

- `git ls-remote --heads origin` returned 17 heads.
- It returned **no** `refs/heads/docs/owner-decisions`.
- It returned `refs/heads/main` at
  `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`.
- It returned `refs/heads/review/od-scope-v25` at
  `9412e943014df6804e34c670d9cc9496e1f34c9c`.
- It also returned
  `refs/heads/review/od-scope-v25-2026-10-22` at
  `e2cd17d0f9c5231c489c42b8dd90d74055d1f648`, which is not the delivery
  branch and is not modified by this iteration.

The local register text is available on this branch but its canonical path is
still absent. The owner's written directions are the controlling instructions
for this follow-up; the register is cited where consistent, and the discrepancy
is preserved rather than resolved by assumption.

## Rules this iteration enforces

- `constitution.md` OBS-153 and CRIT-154: no unreviewed change to the domain
  metrics. The LAW-100 correction is an owner ruling returning runtime behavior
  to the published canonical metric, not a threshold experiment.
- EC §§2, 4, 5, 6, 7, 12, 13, 14: preflight before edit; distinguish measured
  from missing; smallest change; evidence bound to the exact candidate; no
  silent semantic change.
- OD-03/04: no pLDDT threshold expansion, no atom-coverage gate substitution,
  keep LAW-155 out of enforced semantics, and preserve honest non-success states.
- OD-13/16: no unpinned data acquisition, no v25/held-out run, and no merge or
  deployment from this branch.
- All new evidence must identify this exact candidate and must not claim that
  the iteration-2 Stability Gate success proves Governance Acceptance.
