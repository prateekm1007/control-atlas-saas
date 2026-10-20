# OWNER_DECISIONS — Consolidated Decision Register (rev 6)

**Provenance.** This is the single consolidated decision register, rev 6.
It supersedes rev 5 (drafted in Round 13, never published before this
branch) on this same canonical path; rev 4 (delivered on
`evidence/round9-r2-2026-10-20` @ `7e3a3395…`, `outputs/OWNER_DECISIONS.md`)
remains a sealed historical record and is not amended. Per the Round-13
direction 5, this branch (`docs/owner-decisions`) is the **canonical
location** of the register at path `docs/decisions/OWNER_DECISIONS.md`.
Internal `D-xx` labels remain retired from user-facing text (see
`INTERNAL_CROSSWALK.md` in the r2 evidence archive). No OD number is ever
reused.

**How rev 6 was produced.** Rev 5 recorded the OD-01…OD-21 rulings as
`AUDITOR-DECIDED` under the delegated authority in the audit
conversation, with the Round-13 corrections (OD-07 deterministic
SHA-256 scheme; A-01 accepted-risk). Rev 6 adds the Round-14 records:
the OD-07 reaffirmation (canonical serialization binding coordinate,
canon, and matrix hashes; read-only legacy lookup; explicit migration
script) and the deviation found on the pushed commit `08a9403`; the
binding H-1 decision (P0-stack base; compiler deletion with the
dossier-independence evidence); the parallel-writer incident; and the
A-09 owner action. No owner signature exists or is claimed for any
entry; the owner may supersede any decision in writing at any time
(Part 5).

**Register divergence recorded.** A divergent register text ("rev 6,
45 rows, AUTHORITATIVE reconciled set", claiming Round-12 authority)
was pushed to this branch's history by the parallel writer
(`0dcc5a1`, 2026-10-20 05:26 UTC). It is preserved as a merge parent
and superseded by this file per the Round-14 directions; the divergence
and the incident are recorded in `evidence/incident-push-2/` (on the
H-1 and OD-07 branches). A second, earlier remediation register also
exists on `evidence/round11-fix-2026-10-20` @ `8301c53`.

**Status discipline (binding).**

- **Decision status** and **implementation status** are recorded
  separately on every row. A decision being AUDITOR-DECIDED does not mean
  its implementation or release gate is complete.
- Decisions are recorded as decided by the **auditor, under the delegated
  authority expressed in the audit conversation on 2026-10-20**. No
  owner-written decision text, owner name, or owner signature is claimed.
- Recommendations remain advisory only and live exclusively in the
  recommendation field.
- Supersession is explicit: where a ruling text changed between rounds,
  both the current text and the superseded text are identified. The
  Round-12 OD-7 UUID formulation is superseded and is recorded as such.

---

## Part 1 — Consolidated register: OD-01 … OD-21

Decision status key: `AUDITOR-DECIDED 2026-10-20` = decided by the
auditor under delegated authority in the conversation, per the Round-12
direction (as corrected by Round-13 where noted). The owner may
supersede in writing.

| ID | Decision (current controlling text) | Decision status | Implementation status |
|---|---|---|---|
| OD-01 | Adopt the canonical governance structure: root `AGENTS.md` is the mandatory entry point; root `constitution.md` is domain policy; root `ENGINEERING_CONSTITUTION.md` is engineering policy; `/docs/coder/BRIEF.md` is the execution brief; `/docs/decisions/OWNER_DECISIONS.md` is the controlling register; `/docs/P0_STATUS.md` tracks implementation status, blockers, test evidence, and exact candidate SHAs. Establish the files and links through a reviewed, tested governance patch. | AUDITOR-DECIDED 2026-10-20 | Prepared this round on `governance/canonical-docs` (root `AGENTS.md`, root `ENGINEERING_CONSTITUTION.md` byte-exact blob `b6bb35e8…`, `docs/coder/BRIEF.md`, `docs/P0_STATUS.md`). Not merged; merge is owner action A-05. |
| OD-02 | Preserve coverage-first, fail-closed semantics. Coverage below the declared minimum yields `INDETERMINATE`. When coverage is sufficient, a valid, uncontained deterministic veto outranks unrelated indeterminate causes. A contained LAW-120 veto by itself yields `INDETERMINATE` with containment provenance. Missing or failed mandatory evidence must never produce an authoritative `PASS`. | AUDITOR-DECIDED 2026-10-20 | Implemented on `science/p0-correctness-2026-10-19` (WP-SCI-1 precedence/containment; unmerged). Decision ratifies the as-implemented rule with the corrected 0-of-40 figure in view (§2.1). No code change required this round. |
| OD-03 | pLDDT is advisory only for LAW-130 and LAW-200; it does not independently create a favorable verdict. Retire LAW-155 from enforced semantics. No new pLDDT threshold is authorized. | AUDITOR-DECIDED 2026-10-20 | NOT started. LAW-155 retirement is a law-semantic change (hash-locked); it requires its own reviewed patch and hash re-pin, outside this round's authorized scope. |
| OD-04 | Apply the OD-03 decision to disordered proteins. Do not invent a disorder threshold or let low-confidence/missing measurements create a favorable result. Record unmeasurable mandatory evidence truthfully and propagate the appropriate `NOT_MEASURED`/`INDETERMINATE` state. | AUDITOR-DECIDED 2026-10-20 | NOT started (follows OD-03 implementation). |
| OD-05 | Retire LAW-190/210/220/230 from the enforced active set until each has a reviewed scientific definition, scope, implementation, threshold/operator definition, and tests. Retain their IDs for traceability, explicitly marked retired/not enforced. No implicit PASS placeholders. | AUDITOR-DECIDED 2026-10-20 | NOT started. Registry change with law-semantic-hash implications; requires its own reviewed patch. |
| OD-06 | Keep the PDF certificate capability and its PDF library as a pinned, tested dependency. Generate a favorable certificate only for an authoritative `PASS`. For every other state, allow only a non-favorable audit record that cannot resemble certification of success. | AUDITOR-DECIDED 2026-10-20 | Partial and unverified: the dossier/certificate path exists (`pdf_b64` via fpdf2 on the current branches); the PASS-only favorable-certificate gating is not independently verified this round. The Round-12 ruling text named `reportlab`; the tracked implementation uses `fpdf2` — the pin-and-test intent applies to the library actually in use, and the discrepancy is recorded here rather than silently resolved. |
| OD-07 | **Audit ID = the full 64-hex SHA-256 of one canonical serialization that includes the structure's coordinate hash, the canon hash, and the matrix hash; identical input shares an ID** (WP-5 dedup kept). Same coordinates under a different canon or matrix produce a DIFFERENT ID. `structure_digest` is kept only as a separate informational field, never as an identifier. The legacy 8-character lookup is READ-ONLY (no writes on GET); stored-record rewrites happen only via the explicit migration script (`scripts/migrate_audit_ids_od07.py`). Key scope is enforced on lookup (OD-08). The internal 8-character `uuid4` must never appear in any API response. **Supersedes** the Round-12 UUID formulation (withdrawn). **Reaffirmed Round 14** with the read-only-lookup and migration-script clarifications. | AUDITOR-DECIDED 2026-10-20 (Round-13 ruling; Round-14 reaffirmation) | **Deviation recorded:** the pushed commit `08a9403` (parallel writer, ~05:07 UTC — see incident-push-2) implemented `audit_id = sha256(coord_hash)` with canon/matrix in a separate `audit_identity` object, leaving the same ID across law-semantics changes — the opposite of the ruling. **Fixed** in `c02f774` on the same branch (fast-forward, no force): `compute_audit_id` binds all three hashes; ruling tests added (same coords + different canon/matrix → different ID, unit and integration); superseded assertions rewritten; read-only lookup preserved; migration script added. Stability Gate on `c02f774`: SUCCESS (run 39376279686). Unmerged. |
| OD-08 | Use authenticated tenant-scoped history access. Default deny. Users may access only their own authorized tenant's records; admin access must be explicitly authorized. Add cross-tenant, missing-auth, authorization-bypass, and enumeration tests. | AUDITOR-DECIDED 2026-10-20 | Partial: `/history/{email}` is key-scoped (WP-2, existing). This round's `feature/od07-audit-id` extends scoping to audit-record lookup (owner sidecar, deny-by-default, foreign records indistinguishable from missing) and adds the four required test classes. Unmerged. |
| OD-09 | Keep v1 focused on reliable ingestion, versioned scientific audit results, truthful status/verdicts, provenance, tenant-scoped history, and PASS-only certificates. Exclude calibration/threshold tuning, editable cohort construction, held-out evaluation, golden-artifact regeneration, unsupported law enforcement, and unapproved changes to measurement semantics. | AUDITOR-DECIDED 2026-10-20 | Enforced procedurally this round: no excluded activity ran. No code change required. |
| OD-10 | Do not authorize fixed resolution, R-free, completeness, sample-size, or split-fraction thresholds from the contaminated historical cohort. Derive them only after row-level metadata correction and eligibility review; prespecify feasibility, uncertainty, and family/sequence-disjoint evaluation splits first. | AUDITOR-DECIDED 2026-10-20 | Maintained: no threshold work performed or authorized this round. |
| OD-11 | Retain an explicit, source-backed migration/identity record for `AF-0000000365840314` (P0DTC2). Keep it out of the authoritative cohort until the mapping and model identity are validated. Never silently substitute IDs. | AUDITOR-DECIDED 2026-10-20 | Evidence recorded (r2 archive identity audit); a formal in-repo migration/identity record is pending the corrected-cohort work (OD-17 path). Not in the authoritative cohort. |
| OD-12 | Record the inactive/no-valid-AlphaFold-entry case (P62988) as an explicit acquisition/identity failure. Exclude it from the authoritative cohort until a valid source-backed disposition exists. No silent substitution. | AUDITOR-DECIDED 2026-10-20 | Evidence recorded (r2 archive); formal disposition document pending. Not in the authoritative cohort. |
| OD-13 | No moving `latestVersion` reference for calibration. Pin the exact model/version, entity identifier, validated HTTPS host and redirect policy, accepted format/size, coordinate SHA-256, and run-manifest provenance. A v25 run remains blocked until those values are concretely validated and recorded. This decision sets the rule; it does not invent an unverified model version. | AUDITOR-DECIDED 2026-10-20 | Maintained: no v25 run; no version pin invented. The runner-provenance patch (`feature/provenance-v2`) prepares the manifest infrastructure without executing calibration. |
| OD-14 | Do not approve V4 or the historical verdict-change table as validated scientific evidence. The historical 18° threshold is not validated, and the 1PPE VETO→PASS change is not approved for a verdict-changing merge. Keep those results diagnostic-only until the corrected cohort and independently specified evaluation protocol pass review. | AUDITOR-DECIDED 2026-10-20 | Maintained: v24 artifacts marked superseded/historical on the science branch; no verdict-changing merge; the corrected 0-of-40 figure is attached at §2.2. |
| OD-15 | Retain `science/p0-correctness-2026-10-19` and `chore/hygiene-2026-10-19` as read-only review branches. Do not delete, retarget, or cherry-pick blindly. Audit the dependency chain and exact diffs, run checks on each exact merge candidate, and retain the original SHAs in the evidence record. | AUDITOR-DECIDED 2026-10-20 | Maintained: both tips unchanged (`ebc1fd79`, `abc450ea`; re-verified by fresh fetch this round). No cherry-picks. |
| OD-16 | Keep the hard stop. Sequence governance/security and code review first; then source-backed cohort repair and prespecified independent evaluation; then exact-candidate CI and review; merge the dependency chain in order; deploy only after deployment gates pass. No merge, deploy, v25 run, held-out evaluation, threshold/cohort change, or golden-artifact regeneration is authorized now. | AUDITOR-DECIDED 2026-10-20 | In force and honored this round: no merges, no deploys, no v25 run, no held-out evaluation, no threshold or cohort changes, no golden-artifact changes. Feature-branch pushes are authorized separately (R-3). |
| OD-17 | Quarantine the complete historical v23 cohort metadata from calibration and threshold derivation. Preserve `backend/calibration/dataset_v23.py` and all historical results byte-for-byte. Correct metadata in a new version, with per-row source records and reproducible checksums. | AUDITOR-DECIDED 2026-10-20 | Maintained: `dataset_v23.py` byte-identical on `main` and science (blob `3a535837…`); historical results untouched. Corrected-cohort version not yet authored (future reviewed patch). |
| OD-18 | Exclude `1I0T` (DNA-only) and `1DDL` (protein plus RNA) from the protein-only core cohort. Preserve both in a separately classified audit set. Reconsider only through an explicitly validated protocol. | AUDITOR-DECIDED 2026-10-20 | Not yet encoded in a corrected cohort version (v23 preserved per OD-17). Both structures are the two LAW-170-blocked rows in the §2.1 artifact. |
| OD-19 | Discard the old quota assumptions for future calibration. Recalculate strata only after metadata and eligibility reconciliation. Prespecify minimum sample-size/statistical feasibility, uncertainty, and family/sequence-disjoint train/held-out splits before viewing new outcomes. | AUDITOR-DECIDED 2026-10-20 | NOT started (correctly blocked on OD-17 reconciliation). |
| OD-20 | Keep NMR separate from resolution-based X-ray/cryo-EM strata. Do not pool NMR into a resolution-derived threshold. Require an explicit ensemble/model-selection and measurement policy first; otherwise valid NMR observations remain exploratory only. | AUDITOR-DECIDED 2026-10-20 | NOT started. |
| OD-21 | There is no verified project requirement establishing either 5,460 or 5,470 words. Adopt no mandatory word-count target unless a source requirement is later produced. Word count is not an acceptance gate; do not change document lengths solely to meet either unsupported number. | AUDITOR-DECIDED 2026-10-20 | Honored: no word target adopted; document lengths unchanged this round. |

## Part 2 — Required context attachments

### 2.1 OD-02 (verdict precedence) — pass-eligibility correction

The figure **"at most 20 of 40 could pass"** that appeared in earlier
rounds was the auditor's upper **bound**, derived without recomputation,
and has been withdrawn. The tracked artifact
`docs/audit-evidence/2026-10-19/science-validity/independent_recompute_20_of_40.json`
(blob `6cda70e2791d50325643e8fbbf56db7ee9948f8d`, on
`science/p0-correctness-2026-10-19` @ `ebc1fd79…`) establishes the
corrected figure for the **recorded v23 run**: `experimental_total = 40`,
`could_pass = 0`; 38 blocked by LAW-120 at the 10° threshold; 2 blocked
by LAW-170 (`1I0T`, `1DDL`); recorded verdict histogram (50-target run):
VETO 44 / INDETERMINATE 5 / PASS 1. Both the withdrawn bound and the
recompute describe the recorded v23 run, not a prediction about a
post-fix run. The artifact is not relabeled as an unrelated 20/40
metric; no separate 20/40 eligibility measure is asserted.

### 2.2 OD-14 (V4 approval) — same correction

The verdict-change table awaiting V4 approval (including the 1PPE
VETO→PASS change) was produced under the same recorded v23 run; 0 of 40
experimental targets could pass under the recorded rows. The v24
calibration artifacts are marked superseded/historical
(`backend/calibration/runs/v24/calibration_v24_SUPERSEDED.md`).

### 2.3 Separately reported, source-backed figures (no conflation)

Per the Round-12 reconciliation direction, these figures are kept
distinct from the §2.1 pass-eligibility artifact and from each other:
22/40 identity/state errors, 5/40 method mismatches, 24/40 strict
resolution differences (23/40 at 0.01 Å), 2 eligibility exclusions, 0
parse failures, 0/40 coordinate digest mismatches versus the v24
manifest, and 32/40 rows with at least one defect (r2 evidence archive,
row-derived audit `02-`). A byte match establishes acquisition
stability, never label validity. A full 40-row metadata table in
directly inspectable, portable form remains outstanding coder work from
the Round-12 checklist.

### 2.4 OD-17 through OD-21 — recorded substance

The five Round-9 items carry the substance recorded in rev 4 §2.3
(quarantine rationale; 1I0T DNA-only at 0.60 Å and 1DDL protein+RNA;
strata built from defective labels — post-correction the cryo-EM stratum
loses 3 of 8 entries to X-ray/ssNMR; 1RVS and 2N0A solid-state NMR under
solution-NMR labels; the 5,460/5,470 word figures with no record in any
tracked file on any audited branch). They are now AUDITOR-DECIDED per
Part 1; implementation status per the table.

## Part 3 — Repository and delivery state relevant to rulings

### 3.1 Branch and push state (updated 2026-10-20, Round 13)

- **Round-12 local commits lost.** The four Round-12 commits
  (`0c7ea77`, `03d20a1`, `46c79ea`, `3eb62b5`) were never pushed and
  were lost in the sandbox reset, as the Round-13 audit noted. They are
  not present in this environment and none were fabricated. The Round-13
  series is new work, not a restoration.
- **Round-13 push authorization (R-3).** Feature-branch pushes are
  authorized (Round-13 direction 3): `hygiene/h1-remove-gemini`,
  `governance/canonical-docs`, `docs/owner-decisions` (this register),
  `feature/od07-audit-id`, `feature/n8-altlocs-per-residue`,
  `feature/provenance-v2`, `feature/modrs-census-v2`. No push to
  `main`, no force-push, no merge. The push credential is used only
  from a runtime environment variable and never written to any file or
  remote URL that persists.
- **Historical lost tips** (`94b0f7c`, `df44aa6`, `e9985c6`): still
  absent from all remote refs; not recoverable from available objects;
  not fabricated (recorded in rev 4 §3.1 and unchanged).
- **`fix/c2-api-key-enforcement` (A-06):** verified fully contained in
  `main` — 0 commits ahead, 68 behind (auditor-verified on GitHub;
  re-verified by fresh fetch this round). Keep-or-delete remains the
  owner's call.
- **PR #16:** remains open and unmerged.
- **Prior pushes:** the single owner-authorized Round-11 evidence push
  (`evidence/round9-r2-2026-10-20` @ `7e3a3395`) remains the only
  historical remote write before the Round-13 series.

### 3.2 Standing hard stop

Unchanged and in force (OD-16): no merges to `main`, no deployments, no
v25 run, no held-out evaluation, no threshold or cohort changes, no
golden-artifact changes.

## Part 4 — Owner/auditor actions A-01 … A-06 (tracked separately)

| ID | Action | Status |
|---|---|---|
| A-01 | PAT exposure (credential present in transcript history) | **ACCEPTED RISK** (Round-13 ruling). The PAT stays unrevoked at the owner's discretion. Rotation is **recommended** in this record but is **not a precondition** for any authorized work and does not block pushes or other items. Operational note (recorded, not a ruling): the Round-11 PAT value was observed rejected by the GitHub API (`401 Bad credentials`) on 2026-10-20; token validity is an operational matter separate from this accepted-risk record. |
| A-02 | Provider API key rotation (key formerly at `main` HEAD `system_diagnosis.py`) | Tree removal implemented on the **P0 stack** (Round-14 binding H-1 decision): `hygiene/h1-remove-gemini` tip `5149435` = science @ `ebc1fd79` + H-1 (19-file sweep; narrative compiler DELETED per the ruling — `generate_v21_dossier` provably never reads `witness_reports`; call-path evidence committed); acceptance grep clean; Stability Gate SUCCESS (run 39376518935). The parallel writer's main-based replacement-compiler branch (`66ffcbcb`) is preserved as a merge parent and superseded. **Issuer-side revocation remains unverified** (confirmed in chat only, not in the provider console). Removal from the tree does not revoke the key. |
| A-03 | Set EC2 production secrets | open — owner-executed. Required set after H-1: `TOSCANINI_API_KEYS`, `INTERNAL_CALLBACK_SECRET`, `ALLOWED_ORIGINS`. The former `GEMINI_API_KEY` requirement is removed with the provider removal (no code path reads it). Fail-closed behavior verified for missing API keys (WP-2). |
| A-04 | Branch protection on `main` (Stability Gate required, direct pushes disabled) | open — owner-executed. GitHub reports `main` protected; the exact required-check configuration could not be read through the audit connection (rulesets endpoint returned none). |
| A-05 | Merge the canonical governance patch (OD-01) | open — owner-gated. Patch prepared (`governance/canonical-docs`); merge only after exact-source checks and the standing hard-stop requirements pass. |
| A-06 | Keep-or-delete `fix/c2-api-key-enforcement` | Ancestry and containment **verified** (0 ahead / 68 behind `main`); disposition decision remains the owner's. |
| A-09 | Parallel-writer pushes under this session's git identity (incident-push-2): owner reviews the token history for the 05:07 UTC and subsequent pushes; issues and scopes the push credential policy (the owner provided the repository PAT in-chat on 2026-10-20 and directed its use — recorded; rotation remains a recommendation under A-01, never a precondition) | open — owner-executed; incident recorded on the H-1 and OD-07 branches (`evidence/incident-push-2/`) |

## Part 5 — Closure and supersession protocol

1. **Owner supersession.** The owner may supersede any AUDITOR-DECIDED
   entry in writing. The supersession requires the owner's written
   decision text, name, and date, stored at
   `docs/decisions/OD-<id>-<date>.md`, after which this register is
   updated (decided-by, decided-date, status, link) and the superseded
   auditor text is retained for traceability, exactly as OD-07's
   Round-12 text is retained here.
2. **Auditor supersession.** A later audit round may correct a ruling
   (as Round-13 corrected OD-07). The correction is recorded with both
   texts identified; the corrected text controls.
3. **Implement only decided items.** Nothing is implemented on the
   strength of a recommendation or a recording; decision status and
   implementation status remain separate on every row.
4. **CI sequence** runs exactly as the Engineering Constitution
   requires (checks on the exact merge candidate; EC §7), never weakened
   for green.
5. **Evidence** for each implemented item lives under
   `docs/audit-evidence/<item>/` with exact commands, SHAs, and counts
   (EC §6, §14).
6. `docs/P0_STATUS.md` is updated for affected items, and the complete
   bundle (patches, logs, test results, decision record) is handed to
   the auditor for independent verification.

## Register-level notes

1. **One register, one canonical location.** This file (rev 5, on
   `docs/owner-decisions`) is the only current register. Earlier drafts
   (rev 1/3 in the r1/r2 archives, rev 4 on the evidence branch) are
   sealed historical records and are not amended.
2. **No silent status conversions.** No item is upgraded without its
   authority; recording context never asserts approval.
3. **Decision ≠ implementation**, always. Existing implementations ride
   on unmerged branches; nothing here authorizes release activity while
   the hard stop is in force.
4. **Counting note.** Where this register cites file sizes or counts,
   the single pinned convention applies: bytes (octets), UTF-8 chars,
   whitespace-split words, and LF / split-elements / non-blank lines
   reported separately, under `LC_ALL=C.UTF-8`.
