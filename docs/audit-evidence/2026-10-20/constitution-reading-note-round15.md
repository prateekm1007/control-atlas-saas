# Constitution Reading Note — Round 15 (2026-10-20)

Session context: execution of the Round-15 directions — the written
report correction (N8 stacking defect; census ignore entry count), the
N8 rebase onto the OD-07-fixed `c02f774`, the CI workflow-trigger fix
on the governance and register branches, the census ignore-entry
verification, and the push freeze (only directions 2 and 3 are pushed).
This note is committed to the first branch pushed this round (the N8
branch) and to the direction-3 branches.

## 1. Sources, verification, and ordering

| Step | Timestamp (UTC) | Evidence |
| --- | --- | --- |
| Fresh fetch of all remote refs | 2026-10-20 ~06:2x | parallel activity observed: `evidence/round11-fix-2026-10-20` advanced to `742e36b` (05:27–05:28 retraction commits); `proposal/ec16-ec26-enforcement` @ `17eb8818` (06:07) read for this note |
| Domain Constitution blob re-verified | 2026-10-20 ~06:2x | `constitution.md` @ `main` → blob `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` (11 lines) |
| Engineering Constitution blob re-verified | 2026-10-20 ~06:2x | `ENGINEERING_CONSTITUTION.md` @ `policy/engineering-constitution` → blob `b6bb35e84956d84a8c9d9e65734c212baf3ed075` (178 lines) |
| Complete re-read of both documents | 2026-10-20 ~06:2x | content unchanged from Rounds 13/14 (blob IDs identical) |
| Proposed amendments read (NOT in force) | 2026-10-20 ~06:3x | `docs/constitution/amendments/2026-10-20-EC16-EC26-PROPOSED.md`, `…META-RULE-AUTHORING-PROPOSED.md`, `…META-DIRECTIVE-LIFECYCLE-PROPOSED.md` on `proposal/ec16-ec26-enforcement` @ `17eb8818` |
| First edit after the re-read | 2026-10-20 ~06:4x | N8 merge preparation; no repo file was modified before this point in this session |

## 2. Domain Constitution (`constitution.md`, blob `211ad2c5…`)

- **OBS-153 — Avidity Quantification. ACTIVE.** Not modified this
  round; the N8 merge carries a re-run verdict diff proving zero
  verdict changes on the pinned cohort, and every hash lock is
  bit-identical.
- **CRIT-154 — Structural Admissibility. ACTIVE.** Not modified this
  round; untouched by identity, altloc-policy, workflow, and
  evidence-tooling changes.

## 3. Engineering Constitution in force (`ENGINEERING_CONSTITUTION.md`, blob `b6bb35e8…`)

- **EC §1 — ENFORCED.** The Round-14 report errors (N8 stack
  certification without checking the base; census count 40 vs 41) are
  corrected in writing rather than defended.
- **EC §2 — ENFORCED.** Bases pinned before work: the N8 rebase merges
  `c02f774` (verified first commit of the fix) into the remote tip
  `b5afd6c`; the workflow fix targets the main-based governance and
  register branches at `f368bd2a` / `0773f0b1`.
- **EC §3 — ENFORCED.** Push freeze honored: only directions 2 and 3
  are pushed; the PAT remains runtime-env-only, single-invocation use;
  every pushed tree scanned for it.
- **EC §4 — ENFORCED.** The N8 merge makes the identity binding
  (coordinate + canon + matrix) hold on the altloc branch too — the
  same structure under different law semantics can no longer share an
  ID anywhere in the P0 stack.
- **EC §5 — ENFORCED.** The merge resolution is minimal and total: the
  two change-sets were disjoint (processor/altloc tests vs identity
  files), so no behavior was silently dropped; the merged tree was
  verified field-by-field (per_residue default; compute_audit_id
  wiring; superseded test absent).
- **EC §6 — ENFORCED.** Full suite on the merged tree (537 passed, only
  the 4 environmental failures, reproduced at base); the cohort verdict
  diff re-run and recorded; CI run IDs reported as **claimed** where
  the auditor cannot verify them (their sandbox gets 403 from the
  Actions API) — per the Round-15 direction 2.
- **EC §7 — ENFORCED.** Append-only pushes (merge commits, no
  force-push); the workflow-trigger fix is itself a CI repair (EC §7:
  required checks must run on the exact candidate).
- **EC §8 — APPLICABLE, NO NEW CLAIMS.**
- **EC §9 — APPLICABLE, NO CHANGE.**
- **EC §10 — ENFORCED (carried).** No narrative provider anywhere in
  the P0 stack after H-1.
- **EC §11 — ENFORCED.** The parallel writer remains active (new
  commits at 05:27–05:28 on `evidence/round11-fix`); recorded, not
  interacted with; attribution stays with owner action A-09.
- **EC §12 — ENFORCED PER ITEM.** The census mapping verification
  records what was checked (41 entries = 40 verified fingerprints + 1
  inert leftover) and what remains (inert-line removal queued behind
  the freeze).
- **EC §13 — ENFORCED.** No gate weakened: the census exemption story
  is now fully documented (each of the 40 fingerprints maps to exactly
  one digest row; the inert line masks nothing).
- **EC §14 — ENFORCED.** Evidence files, exact SHAs, commands, and
  counts bound to the proposed commits.
- **EC §15 — NOTED, NOT CLAIMED.**
- **Final Rule — ENFORCED.** Evidence before claims; corrections
  issued in writing before new work was claimed.

## 4. Proposed amendments read this round (NOT in force — pending owner approval)

Read in full from `proposal/ec16-ec26-enforcement` @ `17eb8818`. They
are cited here for awareness and are NOT treated as binding; nothing in
this round's work invokes them as authority. Where this round's work
happens to align with one, the alignment is noted as coincidence of
intent, not compliance with a rule.

- **EC-16 (remote-verifiable claims).** Result counts only if its
  commit is on a remote with its SHA in the report. *This round:* all
  pushed SHAs are remote-verified post-push; run IDs are labeled
  "claimed" because the auditor cannot independently read them.
- **EC-17 (push on completion).** *This round:* the N8 rebase and the
  workflow fixes are pushed immediately on completion (also required by
  the Round-15 directions themselves).
- **EC-18 (append-only evidence).** *This round:* evidence files are
  additive; nothing previously committed is edited.
- **EC-19 (one active directive).** This report names its directive
  (Round 15) and the items it closes (directions 1–6).
- **EC-20 (rulings copied verbatim).** *This round:* the OD-07 ruling
  text is carried verbatim in the merge commit and the tests assert it
  (same coords + different canon/matrix → different ID).
- **EC-21 (declared base).** *This round:* every branch's base is
  declared in its commit message (N8: b5afd6c + c02f774; governance:
  main @ 2195a048; register: main @ 2195a048).
- **EC-22 (status vocabulary).** *This round:* statuses are stated as
  RUN_PASS-with-run-ID (marked claimed), NOT_RUN, or PUSHED; the
  auditor's decision explicitly invokes "EC-22: no progress without a
  run" as the rationale for the workflow fix.
- **EC-23 (full-tree secret scan).** The census exemption remains
  scoped to one file/rule with the full 41-entry accounting recorded.
- **EC-24 (push log).** No server-side enforcement point exists on a
  personal repository (the proposal itself classifies the
  only-authorized-names clause as a candidate ADVISORY); the push
  history remains reconstructable from commit metadata and this
  worklog. A-09 stays with the owner.
- **EC-25 (no deploy on push).** The `main` deploy-on-push violation is
  recorded on the proposal branch; remediation is owner action A-10
  (disable in the Actions UI). This round changes no deploy workflow.
- **EC-26 (rules map to acceptance).** This note maps every in-force EC
  rule to its check this round, per the standing instruction.
- **Meta R-AS1..R-AS4 (rule authoring)** and **R-DL1..R-DL5 (directive
  lifecycle)**: read; not in force; no action taken under them.

## 5. Standing constraints honored

- Push freeze: only directions 2 (N8 rebase) and 3 (workflow triggers)
  are pushed; the census mapping evidence and the inert-line removal
  stay local until the freeze lifts.
- Hard stop (OD-16): no merges to `main`, no deploys, no v25 run, no
  held-out evaluation, no threshold/cohort changes, no golden-artifact
  regeneration; C3–C11 untouched.
- One ID, one meaning; no `D-xx` in user-visible files.
