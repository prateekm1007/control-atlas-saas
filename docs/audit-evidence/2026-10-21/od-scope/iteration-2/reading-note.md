# OD-scope iteration 2 — constitution-first reading note

Recorded: 2026-10-21T21:09:13Z
Repository: `https://github.com/prateekm1007/control-atlas-saas.git`
Branch: `review/od-scope-v25`
Base and HEAD: `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
Index: the pre-ruling staged candidate; not committed, not the delivery candidate
Worktree: clean relative to the index at the time of reading; new iteration-2 evidence is untracked

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

- Complete-read timestamp: `2026-10-21T21:08:25Z` (file access times; the
  register access is the last of the six).
- Hash-verification timestamp: `2026-10-21T21:09:13Z`.
- Remote-ref verification timestamp: `2026-10-21T21:09:13Z`
  (`git ls-remote --heads origin`, 15 heads; output saved at
  `remote-heads.txt`).
- First edit in this iteration: `2026-10-21T21:13:45Z` (preflight and ruling
  evidence written after this reading note; no engine code changed before the
  reading was complete).

The recorded test count does not waive this protocol; the reading was
completed before any change in this iteration.

## Active decision IDs

The following recorded requirements are treated as binding for this
iteration, per the owner's direction and the available register text:

- `OD-01` — canonical governance layout.
- `OD-02` — coverage-first, fail-closed verdict precedence.
- `OD-03` — pLDDT advisory only for LAW-130/LAW-200; retire LAW-155 from
  enforced semantics; no new pLDDT threshold.
- `OD-04` — apply OD-03 to disordered proteins; honest
  `NOT_MEASURED`/`INDETERMINATE`.
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

`AGENTS.md:21-23` says:

> `docs/decisions/OWNER_DECISIONS.md` — the controlling decision register
> … Until the governance series is merged, the canonical copy lives on the
> `docs/owner-decisions` branch.

Live verification at `2026-10-21T21:09:13Z`:

- `git ls-remote --heads origin` returned 15 heads.
- It returned **no** `refs/heads/docs/owner-decisions`.
- It returned **no** `refs/heads/review/od-scope-v25`.
- The live `main` is `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`.
- The only relevant remote feature refs are
  `engine/od-scope-v25` @ `97d62a89f0efb2917e80a82808a930d34af7f492`
  and `fix/law-130-clash-boundary` @
  `5c216a882972d56f20be5ddae22bf0557465fdbd`.

The local `docs/decisions/OWNER_DECISIONS.md` is therefore the available
register text on this branch, but it is **not** declared canonical here.
The missing remote branch is an unresolved governance discrepancy, not an
assumption that the local copy is authoritative. The owner's written
directions in this iteration are the controlling instructions for the
work; the register text is cited where it agrees, and the conflict is
preserved in `owner-ruling-2026-10-21.md`.

## Pre-ruling state exported

The staged pre-ruling candidate was exported before modification:

- Directory: `docs/audit-evidence/2026-10-21/od-scope/pre-ruling-candidate/`
- `staged.patch` — complete staged binary diff.
- `unstaged.patch` — complete unstaged diff (empty at export time).
- `status-v2.txt` — `git status --porcelain=v2`.
- `changed-files.tsv` — staged name-status manifest.
- `source-patch-object.txt` — `44efa2c397cad9ce86fc5e2f676e3a040c1d0324`
  (`git diff --cached --binary | git hash-object --stdin`).

That staged tree is **not** committed as the delivery candidate. It is
retained as historical evidence of the rejected scope.

## Rules this iteration enforces

- `constitution.md` OBS-153 and CRIT-154: no change to the domain metrics
  in this work.
- EC §2, §4, §5, §6, §7, §12, §13, §14: preflight before edit; distinguish
  measured from missing; smallest change; evidence bound to the exact
  candidate; no silent semantic change.
- OD-03/04: no pLDDT threshold expansion, no atom-coverage gate
  substitution, LAW-155 out of enforced semantics, honest non-success
  states.
- OD-16: isolated feature branch only; no merge/deploy/calibration.
