# Constitution Reading Note — Round 14 (2026-10-20)

Session context: execution of the Round-14 audit directions — the
incident record for the unexplained 05:07 UTC push of
`feature/od07-audit-id` @ `08a9403`, the OD-07 compliance fix on top of
that commit, the H-1 rebase onto the P0 stack
(`science/p0-correctness-2026-10-19`), SHA re-derivation, bundle
verification, and the corrected report. This note is committed to the
first branch pushed this round (and to the other rebuilt branches), per
the auditor's instruction to cite every EC, OBS, and CRIT ID enforced.

## 1. Sources, verification, and ordering

| Step | Timestamp (UTC) | Evidence |
| --- | --- | --- |
| Fresh fetch of all remote refs (sandbox persisted from Round 13; repo re-fetched) | 2026-10-20 ~05:4x | `git fetch '+refs/heads/*:refs/remotes/origin/*'` — three new remote branches appeared: `feature/od07-audit-id`, `feature/n8-altlocs-per-residue`, `evidence/round11-fix-2026-10-20` |
| Domain Constitution blob re-verified | 2026-10-20 ~05:5x | `constitution.md` @ `main` (2195a048) → blob `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` (11 lines) |
| Engineering Constitution blob re-verified | 2026-10-20 ~05:5x | `ENGINEERING_CONSTITUTION.md` @ `policy/engineering-constitution` (4993cb3b) → blob `b6bb35e84956d84a8c9d9e65734c212baf3ed075` (178 lines) |
| Complete re-read of both documents | 2026-10-20 ~05:5x–06:0x | Both blobs extracted and read in full; content unchanged from the Round-13 read (blob IDs identical) |
| First edit after the re-read | 2026-10-20 06:0x | H-1 P0-stack branch edit; no repo file was modified before this point in this session |

Constitution-first ordering observed: no write, patch, commit, or
test-driven code change preceded the re-verification and re-read.

## 2. Domain Constitution (`constitution.md`, blob `211ad2c5…`)

- **OBS-153 — Avidity Quantification. ACTIVE.** Contact-density standard
  (ρ(Bipod) ≥ 2.5 × ρ(Monopod), unique heavy-atom contacts within
  4.5 Å). Not modified this round; every patch carries a verdict-impact
  statement (none) and the law canon/matrix/law-semantic hashes are
  verified bit-identical on every science-based branch.
- **CRIT-154 — Structural Admissibility. ACTIVE.** Two-Stage
  Warhead-RMSD < 1.0 Å vs PDB 5NIU; steric-clash veto at < 1.2 Å. Not
  modified this round; untouched by the identity, hygiene, provenance,
  and documentation patches.

## 3. Engineering Constitution (`ENGINEERING_CONSTITUTION.md`, blob `b6bb35e8…`)

- **EC §1 (Purpose and authority) — ENFORCED.** The Round-13 report's
  false claim ("PAT is dead / pushes blocked") is recorded and corrected
  rather than papered over; the parallel-writer anomaly is recorded with
  evidence, not guessed at.
- **EC §2 (Mandatory pre-coding protocol) — ENFORCED.** Exact base SHAs
  pinned per branch (`science/p0-correctness-2026-10-19` @ `ebc1fd79`
  for H-1; `feature/od07-audit-id` @ `08a9403` for the OD-07 fix);
  affected callers and consumers re-read on the new bases (science's
  19-file provider footprint; the pushed branch's audit_store/main/tests).
- **EC §3 (Security) — ENFORCED.** The push credential is used only
  from a runtime environment variable within a single shell invocation
  (the Round-13 failure mode — an empty credential from a
  non-propagated variable — is fixed and documented); every pushed tree
  is scanned for the token before pushing; no secret values appear in
  the incident record or any committed artifact.
- **EC §4 (Product and scientific integrity) — ENFORCED.** The OD-07 fix
  makes `audit_id` the SHA-256 of one canonical serialization binding
  the coordinate hash, canon hash, and matrix hash, so the same
  structure under different law semantics can never share an identity;
  the deletion of the narrative compiler follows the evidence (the
  dossier generator provably never reads `witness_reports`), and the
  dashboard's fallback state is the honest degraded state (EC §10).
- **EC §5 (Architecture) — ENFORCED.** Smallest complete change per
  branch; the pushed branch's own key-scoping design
  (`created_by_key_hash`) is preserved in the OD-07 fix rather than
  replaced, changing only what the auditor flagged; the deleted
  compiler's consumers were enumerated (main.py payload construction;
  dashboard with fallbacks) before deletion.
- **EC §6 (Testing and evidence) — ENFORCED.** The ruling test is added
  (same coordinates, different canon hash → different ID; likewise
  matrix); the pushed branch's superseded assertion
  (`test_audit_id_still_sha256_of_coord_hash`) is rewritten to assert
  the ruling; every claimed result carries exact command, SHA, and
  counts; the Round-13 mis-report is corrected with its root cause.
- **EC §7 (CI/CD) — ENFORCED.** All work on feature branches; no push
  to `main`, no force-push, no deletion of the pushed branch; CI runs
  on the exact pushed SHAs and is reported as-is.
- **EC §8 (Reliability/operations) — APPLICABLE, NO NEW CLAIMS.**
- **EC §9 (Customer experience) — APPLICABLE.** The dashboard's
  narrative panel degrades to its built-in pending state — an honest
  degradation, never a fabricated narrative.
- **EC §10 (LLM/external-provider boundaries) — ENFORCED.** The
  external narrative provider is removed from the P0 stack per the
  binding ruling; provider unavailability can no longer fabricate or
  alter anything because the provider path no longer exists.
- **EC §11 (Change, release, and incident discipline) — ENFORCED.**
  The 05:07 push anomaly gets a containment-first incident record with
  evidence preservation (no force-push, no deletion), root-cause
  analysis of my own report error, and remediation (this round's
  corrections); the exposed provider key remains
  treated-as-compromised until issuer-side revocation is verified
  (A-02).
- **EC §12 (Definition of done) — ENFORCED PER ITEM.** This round's
  reports distinguish completed / blocked / unverified; "not run" is
  never reported as "passed."
- **EC §13 (Prohibited shortcuts) — ENFORCED.** No gate weakened; no
  secrets committed; no law meaning, threshold, or ID redefined (the
  OD-07 fix IMPLEMENTS the ruling; it does not change any LAW); no
  history rewritten on any pushed branch.
- **EC §14 (Required task artifacts) — ENFORCED.** Incident record,
  reading note, evidence directories, exact SHAs, and test results
  bound to the proposed commits.
- **EC §15 (External references) — NOTED, NOT CLAIMED.**
- **Final Rule — ENFORCED.** Evidence before claims; the corrected
  report retracts what was wrong and states only what is verifiable.

## 4. Standing constraints honored

- Hard stop (OD-16): no merges, no deploys, no v25 run, no held-out
  evaluation, no threshold or cohort changes, no golden-artifact
  regeneration.
- A-01: accepted risk; the PAT stays unrevoked at the owner's
  discretion; no revocation pestering (owner instruction).
- One ID, one meaning; `D-xx` stays retired.
- No force-push; no branch deletion; the pushed `08a9403` is preserved
  and built upon (direction 1).
