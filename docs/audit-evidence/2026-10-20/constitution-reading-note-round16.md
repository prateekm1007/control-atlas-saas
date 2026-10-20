# Constitution Reading Note — Round 16 (2026-10-20)

Session context: execution of the Round-16 directions under the amended
freeze ("correction pushes are permitted"): (1) the reading-note
corrective commit on this branch, (2) replacing
`.github/workflows/toscanini-gate.yml` with the science-head version,
(3) the census `.gitleaksignore` leftover cleanup on
`feature/modrs-census-v2`, and (4) the written report correction. This
note is the direction-1 corrective commit for this branch. The
direction-1 subject — the Round-15 reading note — is verified present
on this branch already (section 3 below records the evidence), so this
commit carries the Round-16 note plus that presence record instead of
duplicating a byte-identical file (a git no-op).

## 1. Sources, verification, and ordering

| Step | Timestamp (UTC) | Evidence |
| --- | --- | --- |
| Fresh fetch of all remote refs | 2026-10-20 ~09:5x | all tips match the Round-15 report: governance `0933871`, register `d7a6ea7`, N8 `c78ef15`, census `f2b4fc5`, od07 `c02f774`, hygiene `5149435`, provenance `6ce5ee8`, main `2195a048`; parallel writer's `evidence/round11-fix-2026-10-20` @ `742e36b` and `proposal/ec16-ec26-enforcement` @ `17eb8818` untouched |
| Domain Constitution blob re-verified | 2026-10-20 ~09:5x | `constitution.md` @ `main` → blob `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` (12 lines incl. title) |
| Engineering Constitution blob re-verified | 2026-10-20 ~09:5x | `ENGINEERING_CONSTITUTION.md` @ `policy/engineering-constitution` → blob `b6bb35e84956d84a8c9d9e65734c212baf3ed075` (178 lines) |
| Complete re-read of both documents | 2026-10-20 ~09:5x | content unchanged from Rounds 13/14/15 (blob IDs identical) |
| Proposed amendments tip re-verified | 2026-10-20 ~09:5x | `proposal/ec16-ec26-enforcement` @ `17eb8818` (unchanged from the Round-15 full read; still NOT in force) |
| First edit after the re-read | 2026-10-20 ~10:0x | this reading note (governance and register branches); no repo file was modified before this point in this session |

## 2. Round-16 pre-coding verification (EC §2)

Bases pinned before work: this branch at `0933871` (governance) /
`d7a6ea7` (register), both descending from `main` @ `2195a048`; the
census branch at `f2b4fc5` (science-based, base `ebc1fd79`). The
workflow replacement source is the science head
`science/p0-correctness-2026-10-19` @ `ebc1fd79`
(`.github/workflows/toscanini-gate.yml` → blob `8b86ed22…`). Every
step the replaced workflow performs was locally verified against this
branch's tree BEFORE the commit (section 5, EC §2/§6): the gitleaks
v8.30.1 `--no-git` scan, the canon/matrix hash computations, and the
pip resolution of `backend/requirements-lock.txt`.

## 3. Round-15 reading note — presence record (EC §1 conflict record)

The Round-16 audit states the Round-15 note exists on
`feature/n8-altlocs-per-residue` only. The tree evidence, gathered
this session before any edit, is:

- `docs/audit-evidence/2026-10-20/constitution-reading-note-round15.md`
  is present on all three Round-15 tips as ONE identical blob:
  `f4ff9ecfe800237fb95769ce12b88147e3c2ace6` (8,090 bytes) — on N8
  `c78ef15`, governance `0933871`, register `d7a6ea7`. `git log
  --follow` attributes it to exactly those three commits.
- Verified three independent ways: (a) `git ls-tree` on refs fetched
  from origin this session; (b) a whole-tree path search (no other
  reading-note paths differ); (c) GitHub raw endpoints returned HTTP
  200 with byte-identical content for all three branch names
  (sha256 `204302673ca9c3d8…` each).
- Re-verification commands:
  `git ls-tree origin/<branch> -- docs/audit-evidence/2026-10-20/` and
  `curl -s https://raw.githubusercontent.com/prateekm1007/control-atlas-saas/<branch>/docs/audit-evidence/2026-10-20/constitution-reading-note-round15.md | sha256sum`.
- Observed context, offered without argument: at the ROUND-14 tips
  (`f368bd2a` governance, `0773f0b1` register) the file is absent —
  only the `2026-10-20` and `round14` notes exist there. A tree
  listing taken at those SHAs would match the audit's observation.

Per EC §1 the conflict is recorded with evidence rather than resolved
by editing: the Round-15 claim ("committed to all three pushed
branches") was accurate at the Round-15 tips, and no duplicate file is
committed here (re-adding an identical blob is a git no-op). The
in-writing correction owed for Round-15 is instead the one the audit
sustained: the governance and register CI runs were FAILING
(`No module named pytest`, inherited install defect) at report time,
and the proposed "remediation via governance merge" was circular —
both corrected in the Round-16 report.

## 4. Domain Constitution (`constitution.md`, blob `211ad2c5…`)

- **OBS-153 — Avidity Quantification. ACTIVE.** Not modified this
  round: no verdict, threshold, cohort, or engine change is in scope;
  the census commit touches only the scan-exemption list, and the
  workflow/docs commits touch no science path.
- **CRIT-154 — Structural Admissibility. ACTIVE.** Not modified this
  round; untouched by reading notes, workflow replacement, and
  scan-config cleanup.

## 5. Engineering Constitution in force (`ENGINEERING_CONSTITUTION.md`, blob `b6bb35e8…`)

- **EC §1 — ENFORCED.** Both Round-16 audit rejections were
  re-verified against the trees before acting: the reading-note
  presence conflict is recorded with evidence (section 3); the CI
  state is re-characterized with locally proven root causes (below).
  No source was edited to make either finding easier.
- **EC §2 — ENFORCED.** Bases pinned (section 2); constitutions
  re-read in full from re-verified blobs before any edit; every
  workflow step pre-verified against this branch's tree:
  gitleaks v8.30.1 `--no-git` → 1 finding (`system_diagnosis.py:8`,
  rule `gcp-api-key` — the known H-1 committed credential, removal
  pending on `hygiene/h1-remove-gemini` @ `5149435`); canon hash
  `6a9cd4b4349b81de` MATCH; matrix hash `ed1731ace09a9a3b…` MATCH;
  `pip install --dry-run -r backend/requirements-lock.txt` →
  `ResolutionImpossible` (gradio==4.31.0/gradio_client pins vs
  websockets==15.0.1 — the WP-7 conflict, fixed only on the science
  stack); `backend/requirements-dev.txt`,
  `governance/law_semantic.sha256`, and
  `compute_law_semantic_hash` absent from this main-based tree.
- **EC §3 — ENFORCED.** No secret is added, and none is suppressed:
  the H-1 finding is disclosed, not exempted — no `.gitleaksignore`
  entry is created on this branch. The census exemption (direction 3)
  remains scoped to one sealed file and one rule, 40 verified
  fingerprints, zero inert entries after the cleanup.
- **EC §4 — ENFORCED.** The workflow commit message states the
  predicted per-job outcomes BEFORE the runs exist (fail-closed
  disclosure); no failing state is converted into a pass.
- **EC §5 — ENFORCED.** Smallest complete changes: one note commit
  (this one), one workflow commit carrying the science-head file
  verbatim (blob `8b86ed22…`, no edits), one census-line removal. No
  scope expansion; the verbatim carry specifically avoids inventing a
  weaker "docs-branch gate".
- **EC §6 — ENFORCED.** Exact commands, versions, and results are
  recorded in the commit messages and the Round-16 report; CI run IDs
  will be reported as pending until the push executes, never as
  passed.
- **EC §7 — ENFORCED.** The gate carried onto this branch is the
  exact science-stack gate (secret scan + lock/dev install + three
  hash locks when the tree provides them). It is not weakened,
  edited, or re-ordered to obtain green: the unit-tests job is
  expected to fail at the first install command on this main-based
  tree, and the secret-scan job is expected to fail on the true H-1
  finding. Those red results are truthful evidence about this
  branch's base (`main` @ `2195a048`); the green run required for the
  A-05 governance merge can only exist on a tree that carries the
  science stack (post-merge) — recorded in the Round-16 report.
- **EC §8 — APPLICABLE, NO NEW CLAIMS.** No reliability claims are
  made for docs-only branches.
- **EC §9 — APPLICABLE, NO CHANGE.** No customer-facing behavior is
  touched.
- **EC §10 — NOTED.** The inherited workflow's install line
  hand-installs `google-generativeai`; the science-head replacement
  installs from the lock instead. The complete provider removal
  (lock-level) is H-1's scope on `hygiene/h1-remove-gemini`, pending
  merge.
- **EC §11 — ENFORCED.** Each commit carries problem statement,
  evidence, and rollback (revert the single commit).
- **EC §12 — ENFORCED PER ITEM.** The Round-16 report lists changes,
  evidence, unresolved risks (H-1 key on main-based trees; unrun CI
  until push), and blocked checks.
- **EC §13 — ENFORCED.** No gate is weakened, deleted, or bypassed to
  make CI green; the known-red outcomes are disclosed rather than
  gamed (interlock with EC §7 above).
- **EC §14 — ENFORCED.** Evidence (this note, scan logs, hash runs,
  pip resolution output, mapping JSON) is bound to the proposed
  commits.
- **EC §15 — NOTED, NOT CLAIMED.**
- **Final Rule — ENFORCED.** Behavior proven locally before it is
  claimed; the impossible (a green science gate on this main-based
  tree without code changes) is reported with evidence instead of
  forced.

## 6. Proposed amendments (NOT in force — pending owner approval)

`proposal/ec16-ec26-enforcement` @ `17eb8818` — tip re-verified
unchanged this session; content read in full in Round 15. EC-16
through EC-26 and meta rules R-AS1..R-AS4 / R-DL1..R-DL5 are NOT
treated as authority. Where this round's work happens to align (e.g.
EC-16 remote-verifiable SHAs; EC-22 runs-with-conclusions; EC-21
declared bases in commit messages), the alignment is a coincidence of
intent with the in-force EC §6/§7/§12, not compliance with a rule.

## 7. Standing constraints honored

- Push freeze (Round-16 amendment): only directions 1–3 are pushed —
  the two docs-branch commits and the census cleanup; nothing else
  moves on the remote.
- Hard stop (OD-16): no merges to `main`, no deploys, no v25 run, no
  held-out evaluation, no threshold/cohort changes, no
  golden-artifact regeneration; C3–C11 untouched.
- No merges anywhere this round; the A-05 governance merge waits for
  A-10, A-09, and a green run on a tree that can produce one.
- One ID, one meaning; no `D-xx` in user-visible files.
