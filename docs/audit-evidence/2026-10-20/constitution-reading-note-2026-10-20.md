# Constitution Reading Note — Round 13 (2026-10-20)

Session context: execution of the Round-13 auditor directions (OD-07 revert to
the deterministic-SHA-256 ruling, H-1 provider-credential removal, R-3 feature
branch pushes, A-01 accepted-risk record, register rev 5, N8 default
confirmation). This note is committed to every feature branch pushed in this
round, per the auditor's instruction to cite every EC, OBS, and CRIT ID
enforced.

## 1. Sources, verification, and ordering

| Step | Timestamp (UTC) | Evidence |
| --- | --- | --- |
| Fresh clone + fetch of all remote refs | 2026-10-20T04:05–04:14Z | clone + `git fetch '+refs/heads/*:refs/remotes/origin/*'`; 21 remote branches enumerated |
| Domain Constitution blob verified | 2026-10-20T04:14:44Z | `constitution.md` @ `main` (2195a048) → blob `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` (11 lines) |
| Engineering Constitution blob verified | 2026-10-20T04:14:44Z | `ENGINEERING_CONSTITUTION.md` @ `policy/engineering-constitution` (4993cb3b) → blob `b6bb35e84956d84a8c9d9e65734c212baf3ed075` (178 lines) |
| Complete read of both documents | 2026-10-20T04:15–04:17Z | Full sequential read of both blobs extracted via `git show` and re-hashed with `git hash-object` (both matched the pinned IDs above) |
| First edit after the read | 2026-10-20T04:40Z (approx.) | First branch edit (H-1 worktree); no repo file was modified before this point in the session |

Constitution-first ordering was observed: no write, patch, commit, or
test-driven code change preceded the complete read of both documents.

## 2. Domain Constitution (`constitution.md`, blob `211ad2c5…`)

- **OBS-153 — Avidity Quantification. Status ACTIVE.** Bipodal success is
  measured by Contact Density: success requires ρ(Bipod) ≥ 2.5 × ρ(Monopod),
  counting unique heavy-atom contacts within 4.5 Å. Enforcement in this round:
  **not modified and not modifiable** — no threshold, operator, or measurement
  semantics work is authorized (standing hard stop; OD-02/OD-10/OD-14 keep
  semantics frozen). Every patch in this round carries a verdict-impact
  statement confirming zero change to measurement semantics.
- **CRIT-154 — Structural Admissibility. Status ACTIVE.** Lead validation
  requires Two-Stage Warhead-RMSD < 1.0 Å against the clinical benchmark
  (PDB 5NIU); any design relying on non-physical target deformation (steric
  clash < 1.2 Å) is discarded. Enforcement in this round: **not modified**;
  the veto semantics are untouched by the ID, ingestion-policy, provenance,
  hygiene, and documentation patches delivered here.

## 3. Engineering Constitution (`ENGINEERING_CONSTITUTION.md`, blob `b6bb35e8…`)

Reconciliation of EC §1–§15 and the Final Rule against this round's work:

- **EC §1 (Purpose and authority) — ENFORCED.** Both constitutions read in
  full before any edit; source conflicts are recorded rather than silently
  resolved. The Round-12 `OD-07` register text (UUID form) conflicts with the
  Round-13 ruling (deterministic SHA-256 form); the Round-13 message is the
  controlling direction and the register records the supersession explicitly.
- **EC §2 (Mandatory pre-coding protocol) — ENFORCED.** Exact base SHAs
  pinned per branch before work (`main` 2195a048, `science/p0-correctness-2026-10-19`
  ebc1fd79); affected callers and consumers read before editing
  (`backend/main.py`, `backend/tos/ingestion/processor.py`,
  `backend/tos/storage/audit_store.py`, `backend/tos/security/auth.py`,
  `backend/tos/saas/key_store.py`, `backend/tests/*`, workflows, compose and
  helm wiring); expected behavior, risks, and acceptance criteria are stated
  in each branch's commit message and evidence files.
- **EC §3 (Security is a release requirement) — ENFORCED.** Fail-closed and
  deny-by-default are preserved or strengthened: H-1 deletes the committed
  provider credential and the runtime falls back to the deterministic offline
  narrative (no open-auth path is introduced); OD-08 key scope is enforced on
  audit lookup (cross-tenant and missing-record denials are
  indistinguishable, resisting enumeration); the push credential is used only
  from a runtime environment variable, never written to any file, and every
  pushed tree is scanned for it before the push.
- **EC §4 (Product and scientific integrity) — ENFORCED.** OD-07 audit IDs
  bind provenance (canonical input digest combined with the law canon hash
  and matrix hash) so an authoritative result remains traceable to the exact
  input and law semantics in force; `structure_digest` stays informational
  and is never an identifier; no error, absence, or timeout is converted into
  a favorable verdict by any patch in this round; the N8 default change
  alters only alternate-conformer selection policy, and its verdict impact on
  the pinned cohort is stated explicitly (diff shown even where empty).
- **EC §5 (Architecture and maintainability) — ENFORCED.** Smallest complete
  change per branch; the provider-compiler replacement keeps the existing
  module contract (`get_compiler()`, `synthesize_dossier_content`,
  `model_used`) so callers and response schemas remain stable; import graphs
  and callers were inspected before deleting or replacing any module; no new
  competing source of truth is introduced.
- **EC §6 (Testing and evidence) — ENFORCED.** `test_ids.py` is updated to
  assert the ruling rather than the superseded assumption; every claimed
  result in this round's report includes the exact command, commit SHA, and
  pass/fail counts; pre-existing baseline failures are reproduced on the
  unmodified base and reported as such, never hidden or skipped to obtain
  green CI.
- **EC §7 (CI/CD and software supply chain) — ENFORCED.** All work is
  developed on feature branches; no direct push to `main`, no force-push, no
  merge; CI runs on the exact pushed SHAs (the Stability Gate triggers on
  push to any branch) and its output is reported as-is, including the
  inherited red baseline of `main`-based branches.
- **EC §8 (Reliability, observability, operations) — APPLICABLE, NO NEW
  CLAIMS.** No new production reliability surface is added this round; the
  provenance manifest records git/dependency context for traceability; no
  availability or SLO claim is made.
- **EC §9 (Customer experience, accessibility, cost control) — APPLICABLE,
  NO CHANGE.** No UI, quota, or cost surface is modified; removing the
  external narrative provider eliminates an uncontrolled external dependency
  rather than adding one.
- **EC §10 (LLM and external-provider boundaries) — ENFORCED.** H-1 removes
  the external narrative provider from the runtime entirely; the
  deterministic narrative remains available when no provider exists; the
  narrative layer never overrides the authoritative engine verdict
  (unchanged, re-verified by existing response-schema tests).
- **EC §11 (Change, release, and incident discipline) — ENFORCED.** Each
  patch states problem, acceptance criteria, risk, tests, and rollback; the
  exposed credential is treated as compromised until its issuer confirms
  revocation (A-02 remains an owner action; removal from the tree does not
  revoke it); no history is rewritten on any pushed branch.
- **EC §12 (Definition of done) — ENFORCED PER ITEM.** Final reports list
  changes, evidence, unresolved risks, and unverified items; "not run" is
  never reported as "passed"; no item is called complete while material
  checks remain unrun.
- **EC §13 (Prohibited shortcuts) — ENFORCED.** No gate is bypassed or
  weakened for green CI; no secrets are committed or echoed; no law meaning,
  threshold, or ID is altered; no module is deleted before dependency
  analysis; no merge or deploy of any candidate in this round.
- **EC §14 (Required task artifacts) — ENFORCED.** Each branch carries its
  baseline SHA, changed-file manifest, machine-readable test results bound to
  the proposed commit, and an impact note (see the per-branch evidence files
  and commit messages).
- **EC §15 (External engineering references) — NOTED, NOT CLAIMED.** NIST
  SSDF, OWASP ASVS, and SLSA are used as guidance mappings only; no
  compliance or certification claim is made for this repository.
- **Final Rule — ENFORCED AS THE ROUND'S OPERATING PRINCIPLE.** Trust before
  features; proven behavior before success claims; evidence preserved before
  changing authority; automatable invariants automated (hash locks, secret
  scan, ID tests).

## 4. Standing constraints honored throughout this round

- Hard stop (OD-16): no merges to `main`, no deploys, no v25 calibration run,
  no held-out evaluation, no threshold or cohort changes, no golden-artifact
  regeneration.
- A-01: recorded as an accepted risk per the Round-13 ruling; the PAT remains
  unrevoked at the owner's discretion and is never a precondition for
  authorized work; rotation is recommended in the record only.
- One ID, one meaning: `OD-xx` identifiers are used only in the sense fixed
  by the register; retired internal `D-xx` identifiers do not reappear in
  user-visible files.
