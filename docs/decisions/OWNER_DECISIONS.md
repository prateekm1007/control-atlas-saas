# OWNER_DECISIONS — Consolidated Decision Register (rev 9)

**Provenance.** This is the single consolidated decision register, rev 9
(Round-20, 2026-10-20). Rev 8 — the Round-19 records (the P-015-CONFIRM
collision and force-update deviation, the repair-forward adjudication, the
live end-to-end enforcement, and the precise activation states per rule) —
was prepared on the `docs/owner-decisions-r19` branch but never entered the
reviewed path (no register PR was open; auditor Round-20 direction 6). Rev 9
carries the FULL rev-8 content into this reviewed path **content-only** (no
merge of the divergent register-branch history; every prior decision and
evidence preserved unchanged) and appends Part 12 — the Round-20 records
(the auditor's "Round-19 is not closed" ruling, the two critical EC-23/EC-24
acceptance-defect fixes, the EC-26 trigger-liveness gap closure, and the
R19 follow-up-push reconciliation). Every row of rev 7 and rev 8 is preserved
unchanged (append-only register discipline; status conversions never happen
silently).
Rev 7 appended the 2026-10-20 auditor adjudications of the
EC-16…EC-26 amendment proposals and the two meta-amendments (Part 6), the
updated owner/auditor action rows (Part 7), the push-authorization
reconciliation (Part 8), and the preserved-blockers statement (Part 9).
Rev 7 superseded rev 6 (Round 14, published as `0773f0b`) on this same
canonical path; rev 6 in turn superseded rev 5 (drafted in Round 13, never published before that
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

1. **One register, one canonical location.** This file (rev 7, on
   `docs/owner-decisions`) is the only current register. (Correction,
   2026-10-20: this note previously said "rev 5" — a stale reference left
   from the rev-6 edit; the register it describes was rev 6 at the time and
   is rev 7 now.) Earlier drafts
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

## Part 6 — EC-16…EC-26 adjudications + meta-amendments (AUDITOR-ADJUDICATED, 2026-10-20)

The user's instruction delegated these policy decisions to the auditor. They are recorded
below as **`AUDITOR-ADJUDICATED` under delegated authority, dated 2026-10-20**. No owner
signature is claimed or implied. These are decisions on the INTENDED POLICY CONTENT: the
amendments are **not active** in the repository Constitution, no workflow is activated, and
nothing is merged until the implementation defects are fixed and the exact merge candidate
passes the required validation gates (auditor direction 2026-10-20, sections 2–3).

| Amendment | Decision | Binding direction (auditor, 2026-10-20) |
|---|---|---|
| EC-16 — Remote-verifiable claims | **ACCEPT WITH REVISIONS** | Require every claimed SHA to be a full 40-character lowercase hexadecimal commit ID. Verify it against the intended repository using GitHub's commit API or another verified remote-object route. Check that the returned commit SHA exactly matches the request. Reject mixed valid/invalid lists; never silently ignore short or malformed IDs. API permission, rate-limit, and network errors must fail closed as verification failures, not be misreported as nonexistent commits. |
| EC-17 — Push on completion | **ACCEPT** | Each completed work item must be durably stored on an authorized remote feature/evidence branch or delivered in a bundle with verified `SHA256SUMS` before subsequent work begins. A proposed rule cannot authorize its own push. Enforce the existing scope of push authorization. |
| EC-18 — Append-only evidence | **ACCEPT WITH REVISIONS** | Retain UTC timestamps, unique entry IDs, monotonic sequence checking, and explicit `RECONSTRUCTED` labels. Add a check that the current log preserves the exact previous version as a prefix, with no edits or deletions to prior entries. Bind recorded actions to source evidence; timestamp comparison alone cannot prove that an action occurred at the claimed time. |
| EC-19 — One active directive | **ACCEPT** | Maintain one active directive series in the single index. Each report and evidence package must identify its directive and the items it closes. Historical backfills must be clearly marked as reconstructed. |
| EC-20 — Rulings copied verbatim | **ACCEPT WITH REVISIONS** | Copy the controlling ruling into the work item and name machine-verifiable acceptance tests. Require the named tests to exist and actually pass on the exact candidate. A file-path presence check is not a test run. Any deviation requires a recorded written waiver before implementation. |
| EC-21 — Declared base | **ACCEPT WITH REVISIONS** | Require a full base SHA and a validated `Stack-parent` value. Check ancestry against fetched, complete history and validate the declared stack class — `main-line`, `science-line`, or `policy-line` — against the expected base. Do not accept a local Git object's existence as proof that it is on the remote. |
| EC-22 — Status vocabulary | **ACCEPT WITH REVISIONS** | Only `RUN_PASS` accompanied by an actual CI run ID counts as CI progress. Retain `NOT_RUN`, `RUN_FAIL`, `PUSHED`, and `MERGED`. Reject unknown status tokens. Progress percentages must be computed from verifiable run records, not trusted because a coder typed a percentage beside a run ID. |
| EC-23 — Full-tree secret scan | **ACCEPT WITH REVISIONS** | Scan the intended repository history and all required refs, including evidence branches. Pin the Gitleaks action to a reviewed immutable commit SHA, declare least-privilege workflow permissions, and test detection and failure behavior using safe fixtures. A scan of only the current push range must not be represented as a whole-repository historical scan. |
| EC-24 — Push log | **ACCEPT WITH REVISIONS** | Require a checkable push log with branch, full SHA, timestamp, and credential identifier where technically observable. Classify "only authorized credential names can push" as `ADVISORY` until a real server-side enforcement mechanism exists. A post-hoc log does not prevent an unauthorized credential from pushing. |
| EC-25 — No ungated deploy on push | **ACCEPT WITH REVISIONS** | Production deployment must be dispatch-only, restricted to the approved deployment ref, gated on exact-SHA validation and successful required checks, with environment protection. No deployment on ordinary `push` or `workflow_run` events. Remediation of the existing workflow must be part of the same reviewed activation plan as the lint. |
| EC-26 — Rules map to acceptance | **ACCEPT WITH REVISIONS** | Derive the in-force rule inventory from the canonical directive index — not an arbitrary list supplied by the pull-request author. For each applicable rule, require a mapping to an actual enforcement point and named acceptance test. Validate that those paths and checks exist. |
| R-AS1–R-AS4 — Rule-authoring standard | **ACCEPT WITH REVISIONS** | Every in-force rule needs a real enforcement point. Mark unenforceable clauses `ADVISORY`; advisory rules do not count toward compliance. Maintain the single classification index. Retrofit existing rules without silently changing their meaning. |
| R-DL1–R-DL5 — Directive lifecycle | **ACCEPT WITH REVISIONS** | Maintain one active series and clear supersession records. Read-only audits, local drafts, and verification may proceed within the active directive. The default must not authorize unlisted remote writes, merges, deployments, credential operations, spending, or scientific-semantic changes. Existing scoped push authorization remains scoped. |

**Adoption order (auditor-decided):** Batch 1 = EC-17 and EC-25, with the rule-authoring
standard and directive-lifecycle safeguards; Batch 2 = EC-22 and EC-23; Batch 3 = EC-16,
EC-18, EC-19, EC-20, EC-21, EC-24, and EC-26.

**Implementation status (separate from decision status, per register discipline):** the
correction cycle implementing auditor items 3.1–3.9 lives on `proposal/ec16-ec26-enforcement`
(inert under `proposals/` and `scripts/ci/`; local acceptance suite AT-01…AT-16 recorded in
`docs/evidence/ec16-ec26/LOCAL_TEST_LOG.md` on that branch; **GitHub Actions CI for those
commits: NOT RUN** — no run ID exists and no CI-pass claim is made). Activation requires:
defect fixes verified, exact merge candidate passing the required validation gates, and the
reviewed activation plan (`proposals/remediation/` for EC-25) applied in order. Nothing here
authorizes activation or merge.

## Part 7 — Owner/auditor actions (rev 7 updates to Part 4; statuses unchanged unless stated)

| Action | Rev 7 direction (2026-10-20) | Implementation status |
|---|---|---|
| A-01 — PAT risk | Retain the existing `ACCEPTED RISK` status. Do not issue another rotation reminder. Never write or reproduce the token value. | UNCHANGED — ACCEPTED RISK; no new rotation reminder issued in this cycle; the token value appears in no file, log, or register row. |
| A-02 — Former Gemini API key | Require issuer-side revocation/rotation evidence. Removal from a feature branch is not key revocation. Keep deployment blocked while the deployed `main` code and secret requirements remain unresolved. | OPEN — tree removal remains implemented on `hygiene/h1-remove-gemini` @ `5149435` (unmerged); issuer-side evidence still not presented; deployment stays blocked. |
| A-03 — EC2 secrets | Require the mandatory production secrets to be configured and verified through the approved secure mechanism; fail closed if any required secret is absent. Do not claim configuration has been completed without evidence. | OPEN — no configuration evidence exists; fail-closed behavior is verified in code (WP-2) but the production secrets remain unconfigured as far as recorded evidence shows. |
| A-04 — Branch protection | Require the Stability Gate and secret scan as required checks and disable direct pushes to `main`. The `main` protected flag alone does not prove the exact required-check configuration. | OPEN — protection flag observed; the exact required-check set remains unread/unverified and is NOT claimed as configured. |
| A-05 — Canonical governance merge | Keep the governance patch unmerged until its exact-source validation and required checks pass. No merge is authorized by this audit. | OPEN — `governance/canonical-docs` remains unmerged; no merge authorized. |
| A-06 — Legacy `fix/c2-api-key-enforcement` branch | Keep it read-only for this task. Its reported zero-ahead/68-behind state supports cleanup, but preserve the ref evidence and confirm no open PR or unique data before any later deletion. | OPEN (read-only) — no writes to that branch this cycle; disposition remains the owner's. |
| A-09 — Parallel-writer incident | Preserve the incident evidence append-only. Do not infer credential identity from Git author/committer fields alone. No further unscoped remote writes by parallel actors; every claim must be supported by exact remote-ref evidence. | OPEN — incident records preserved; no attribution added beyond recorded evidence; no unscoped writes performed in this cycle. |

## Part 8 — Push-authorization reconciliation (honest record)

- Rev 6's authorization list (Round-13 direction 3) named seven branches:
  `hygiene/h1-remove-gemini`, `governance/canonical-docs`, `docs/owner-decisions`,
  `feature/od07-audit-id`, `feature/n8-altlocs-per-residue`, `feature/provenance-v2`,
  `feature/modrs-census-v2`. It did **not** list `proposal/ec16-ec26-enforcement`, yet the
  proposal package was pushed to that branch (push-log P-011, tip `17eb881…`).
- **Discrepancy recorded, not smoothed over:** that push preceded this authorization text and
  was NOT covered by rev 6's list. It was not self-authorized by EC-17 (which is a proposal,
  not an in-force rule — a proposed rule cannot authorize its own push).
- **Authorization now in force (auditor direction 2026-10-20, section 3):** continued
  focused, auditable corrections on `proposal/ec16-ec26-enforcement` — appended commits only,
  no force-push, no merge to `main`, no activation of workflows — are authorized, and the
  branch is hereby added to the standing push-authorization list for that scoped purpose.
  The two-phase push log (INTENT + CONFIRM records, `PUSH_LOG_SCHEMA.md` on the proposal
  branch) documents each subsequent push of this cycle.
- Still in force: no push to `main`; no force-push of any published branch; no merge; no
  deployment; the push credential is used only from a runtime environment variable / protected
  file outside the repository and its value is never written to any file, log, or URL.

## Part 9 — Preserved blockers (explicitly NOT closed by this cycle)

1. **Portable 40-row cohort metadata table** — OUTSTANDING (unchanged from rev 6). The
   directly inspectable, portable, row-level cohort metadata table has not been produced;
   nothing in the amendment work marks it complete.
2. **r2 archive acceptance** — BLOCKED on the auditor side: the delivered archive bytes must
   be independently downloaded, hashed, extracted, and verified in a clean environment before
   acceptance; the coder-side verifications previously delivered do not substitute for that.
3. **v23 dataset and historic calibration results** — preserved unchanged; metadata
   quarantine, the OD-17 corrected-cohort process, and the scientific hard stop remain in
   force; no v25 calibration, no held-out evaluation, no threshold/cohort changes, no
   golden-artifact regeneration.
4. **Hard stop (OD-16)** — unchanged: no merges to `main`, no deployments, no unlisted
   remote writes, no credential operations.

## Part 10 — Concurrent-writer reconciliation during this cycle (recorded 2026-10-20)

While rev 7 was being prepared, the remote `docs/owner-decisions` branch advanced
`0773f0b` → `a99348f` (three commits, committer date 2026-10-20T10:03:38Z, identity
`Z User <z@container>` — the same sandbox identity documented in the incident records)
carrying a `toscanini-gate.yml` trigger fix and two constitution reading notes labeled
Round-15/Round-16. The concurrent work cites "Round-16 direction 2"; that directive text
is not present in this environment and is NOT evaluated here. Reconciliation = a merge
(`b421237`), not a rewrite: file sets are disjoint, both lineages preserved, no
force-push, per the A-09 discipline. Attribution of the concurrent pushes follows A-09:
git identity alone does not establish credential identity; the owner's token-history
review remains the authority.

## Part 11 — Round-19 mainline repair records (AUDITOR-ADJUDICATED, 2026-10-20)

Recorded from the auditor's Round-19 directions 1–8 (full delegation, same
discipline as Part 6: decision status and implementation status are separate;
nothing here claims an owner signature).

### 11.1 Merge-collision adjudication (direction 2) — APPLIED

Current `main` @ `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` is the baseline and is
repaired FORWARD through a reviewed, narrowly scoped PR
(`repair/round19-mainline-acceptance`); it is not replaced by or merged with
`c4a0fcc` (preserved as evidence; no further force-push). The
`a0e321c -> c4a0fcc` force-with-lease update is recorded as a deviation from
rev 7 Part 8 (appended commits only) in
`docs/evidence/push-log/INCIDENT-2026-10-20-P015-COLLISION.md`, together with
the `P-015-CONFIRM` identifier collision between the two competing histories.
No push is attributed to a person from git metadata alone (A-09 discipline).

### 11.2 Mainline acceptance defects (direction 3) — REPAIRED on the branch

- 3A: the two remaining non-documentation provider-SDK references removed
  (reviewed content carried from `c4a0fcc` content-only).
- 3B: the remediation README corrected to the current live state; patches
  labeled superseded/historical with exact blob links; policy status and
  active enforcement recorded separately.
- 3C: package manifest scope defined and machine-enforced
  (`scripts/ci/manifest_check.py`); `SHA256SUMS` regenerated over the full
  38-file scope; verified in a clean checkout (all 38 digests OK via
  `sha256sum -c`).

### 11.3 EC-20 activation failures (direction 5) — ROOT CAUSE FOUND AND FIXED

PRs #18/#19/#20 all fail their EC-20 job (runs 39417485854 / 39417488957 /
39417492355) because `verify_claimed_shas.py --check-tests` collected tests
repo-wide; legacy script-style modules under `gpu_worker/tests/` execute at
import time and `sys.exit()` at module level, aborting pytest collection with
INTERNALERROR so the named acceptance tests registered as "not collectable".
Reproduced locally on the exact head `748544f` in a pytest-only environment
before fixing; the fix (per-file scoped collection, AT-18) was verified to
execute the named tests on that same head. The three duplicate workflow
additions are consolidated into one activation candidate
(`activation/ec16-ec19-ec20-r19`); PRs #18/#19/#20 are superseded (branches
preserved as historical evidence).

### 11.4 Activation states (direction 6) — precise, as of this register revision

| Rule | Policy status | Implementation status |
|---|---|---|
| EC-17 / EC-22 / EC-25 | ACTIVE (in force, prior ruling) | Live automated enforcement begins with the Round-19 repair's `governance-acceptance.yml` (before it: policy in force, checks inert — recorded, not smoothed over) |
| EC-16 / EC-19 / EC-20 | PROPOSED — NOT active | Fixes landed in the repair; consolidated activation candidate prepared (`activation/ec16-ec19-ec20-r19`); activation requires the complete gate green on the exact candidate AND the auditor's register recording |
| EC-21 | PROPOSED — NOT active (accepted in principle) | Original push confirmed against the live remote ref (`7f5d50c`, run 39417494195 on that head); rebased candidate `activation/ec21-r19` prepared with the reconciliation record; same activation condition |
| EC-26 | PROPOSED — NOT active | Checker gap fixed in the repair (ACTIVE rules must cite live workflows; AT-17 negatives); candidate `activation/ec26-r19` prepared; same activation condition |
| EC-18 / EC-23 / EC-24 | PROPOSED — NOT active | Candidates prepared (`activation/ec18-r19`, `activation/ec23-r19`, `activation/ec24-r19`); EC-23's live job runs the checksum-pinned scanner with positive/negative controls and the documented H-1 history baseline (A-02 stays open); same activation condition |

### 11.5 Live end-to-end enforcement (direction 7) — IMPLEMENTED on the branch

`.github/workflows/governance-acceptance.yml` (push + pull_request, all
branches; the Stability Gate is untouched) runs on every relevant PR and push:
the full acceptance suite (AT-01..AT-20), reading-note/inventory validation
with the ACTIVE-liveness rule, status-vocabulary validation against an
in-job runs manifest built from Actions records, push-log schema + unique-ID +
pairing + append-only parent-prefix validation, the deploy-trigger lint
against the actual `.github/workflows/` tree, the package manifest +
forbidden-reference checks, PR-body claim verification with EXECUTED named
acceptance tests, declared-base verification, and the full-history all-refs
secret scan with controls. A green base gate remains necessary, not
sufficient.

### 11.6 Standing items (direction 8) — UNCHANGED

A-01 remains ACCEPTED RISK (no rotation reminder issued; the token value
appears nowhere). The scientific hard stop is preserved in full: no law
semantics, calibration, thresholds, cohorts, held-out evaluation, golden
artifacts, or deployment. The portable 40-row cohort metadata table and the
independent r2 archive acceptance remain OPEN and are not closed by the
amendment or repair work. Merges of the repair and the activation candidates
go through reviewed PRs on exact-candidate green checks; no force-push; no
administrative merges.

## Part 12 — Round-20 records (AUDITOR-ADJUDICATED, 2026-10-20)

Recorded from the auditor's Round-20 ruling ("Round-19 is not closed") and
directions 1–8, under the same delegated authority and register discipline as
Part 6/Part 11: decision status and implementation status remain separate; no
owner signature is claimed; nothing here converts a status silently.

### 12.1 The ruling and the two critical deficiencies

The auditor independently verified the live refs and ruled Round-19 NOT
closed despite green checks on PR #30's then-head `f44bf4d…` (gate run
`39443515998`, Governance Acceptance run `39443516042`): green is necessary,
not sufficient. Two critical deficiencies were found in the new controls,
plus one EC-26 gap and one register gap:

1. **EC-23 false-positive-classification defect.** The all-refs scan
   tolerated EVERY `generic-api-key` finding under `docs/` or `evidence/` as
   an evidence-digest class — a category exemption under which a newly
   introduced credential in those directories would also pass. Separately,
   the positive control treated any nonzero scanner exit as proof of
   detection, so a scanner crash or invocation error could masquerade as a
   successful control.
2. **EC-24 verification defect.** The schema checker validated field syntax
   but never verified that a recorded `remote_sha` equals the live remote ref
   for the stated branch, nor that the `evidence=` field supports the recorded
   SHA. Materially: the P-017…P-024 confirmations record the earlier
   first-push tips while the current PR heads are different tips
   (`f44bf4d`, `bd8d86a`, `32d272d`, `151b8c1`, `97067a0`, `fa854f8`,
   `28ee983`, `df98f5e`); those confirmations are valid for the pushes they
   describe but do not confirm the newer tips, and the pushes that carried
   them were not themselves reconciled. The Round-19 report's claim that the
   push protocol was fully closed was therefore not supported.
3. **EC-26 gap.** The liveness test checked the workflow's path and job
   existence but not that its `on:` triggers actually fire on the events the
   claimed enforcement requires (a manual-dispatch-only workflow under
   `.github/workflows/` with the named job would have satisfied it).
4. **Register gap.** Register rev 8 sat on `docs/owner-decisions-r19` outside
   the reviewed path; PR #30 still carried rev 7 and no register PR was open.

### 12.2 EC-23 fix — reviewed fingerprint allowlist + verified controls

- The category exemption is REPLACED by an explicit, reviewed allowlist:
  `docs/evidence/ec16-ec26/secret_scan_baseline.json` (47 documented
  fingerprints: source, classification, and reason per entry; documented in
  `SECRET_SCAN_BASELINE.md` in the same directory) enforced by
  `scripts/ci/secret_baseline_check.py`. Unknown fingerprints FAIL wherever
  the file lives — including `docs/` and `evidence/`. A fingerprint that
  matches while its rule/file contradict the reviewed entry FAILS. The
  baseline is validated itself (required fields, unique fingerprints, known
  classification tokens, H-1 entries require the A-02 note).
- The two historical H-1 findings (the former provider key at pre-release
  commit `b66f0e3b…`) remain a SEPARATE `DOCUMENTED-HISTORICAL-A-02` class,
  disclosed on every run and tolerated only while register A-02 stays open
  (issuer-side revocation unverified; deployment blocked; history never
  rewritten as a substitute for revocation).
- The positive control now runs with `--exit-code 101` (findings-found is
  distinct from execution errors, which exit 1) and is verified by
  `scripts/ci/secret_control_check.py`: the report must exist and parse; the
  expected finding for the seeded synthetic key must be present; the expected
  rule (`aws-access-token`) must actually have fired; no unexpected findings;
  and the assembled key must be redacted everywhere (absent from report bytes
  and the captured scanner log; the key travels via an environment variable,
  never argv or file content). The clean-repository negative control must
  exit 0 with an empty report. Acceptance tests: AT-21 (allowlist), AT-22
  (controls).
- The all-refs scan itself runs with `--exit-code 0`, so any nonzero exit is
  a scanner execution failure and fails closed (the `|| true` masking is
  removed).

### 12.3 EC-24 fix — remote-ref verification and truthful reconciliation

- `scripts/ci/evidence_sequence_check.py --remote-verify` now: (a) fetches
  the live ref (`git ls-remote`) for every CONFIRM/RECONCILE's claimed branch
  and compares the exact SHA — equal = CURRENT, different = HISTORICAL
  (valid for the push it describes, explicitly NOT a confirmation of the
  newer tip); (b) requires the `evidence=` field to support the recorded SHA
  (the SHA and the claimed branch must appear in the evidence or entry body;
  a cited Actions run ID whose manifest `head_sha` contradicts the recorded
  SHA FAILS — fabricated support); (c) requires the branch under test's live
  tip to be DECLARED — by an unpaired INTENT (PENDING, disclosed, never
  implied complete) or an exact-tip CONFIRM/RECONCILE — with an undeclared
  tip (the Round-19 follow-up-push class), an ls-remote failure, or an absent
  remote ref failing closed; (d) classifies missing confirmation as PENDING,
  never backdated or fabricated.
- The eight Round-19 follow-up pushes are reconciled TRUTHFULLY by
  `push-reconcile` records **P-050…P-057** (schema v3, defined in
  `PUSH_LOG_SCHEMA.md`; narrative in
  `RECONCILIATION-2026-10-20-R19-FOLLOWUP-PUSHES.md`): dated, live-ref
  verified at the reconciliation window `2026-10-20T15:58:48Z`, with
  check-run citations for each exact head, and an explicit disclosure that no
  INTENT existed for those pushes (recorded, not smoothed over). No prior
  entry was edited; the append-only parent-prefix check stays in force.
- The Governance Acceptance `push-log-integrity` job now runs the full
  verification on every push/PR event with the event's exact branch and tip.
  Acceptance tests: AT-24.

### 12.4 EC-26 fix — trigger liveness

ACTIVE rules' `workflow:` references must now declare the events their
enforcement claims (`~events=push+pull_request` suffix in the INDEX
inventory), and the live workflow's `on:` block must actually contain each
declared event with branch filters covering the enforcement branch (`main`).
Manual-only (`workflow_dispatch`-only), branch-disabled, and otherwise
untriggered workflows FAIL; the previous AT-17 path/job checks are retained.
Acceptance tests: AT-23. EC-26 remains PROPOSED until the revised checks
pass on the exact candidate (auditor direction 5).

### 12.5 Register discipline for this revision (status distinctions)

Per auditor direction 6, this register distinguishes, for every rule touched
this round: the auditor's decision status (recorded here), the policy status
(ACTIVE/PROPOSED in the INDEX inventory), the implementation status (what is
built and where), the exact candidate (the repair-branch tip carrying this
revision — recorded in the push log P-058 records and PR #30's body, and
verified by the claim-verification and push-log-integrity jobs on that exact
SHA), the CI results on that exact SHA (delivered by those jobs' runs; run
IDs recorded in the push-log confirmations and the round report), and whether
a rule's automation is LIVE ON `main` versus merely present on a PR branch
(nothing in this Part is live on `main`; `main` remains at `f535a71…` until
PR #30 merges through the reviewed path).

### 12.6 Activation adjudications applied (auditor direction 7)

| Rule | Decision status | Policy status | Implementation status |
|---|---|---|---|
| EC-16 / EC-19 / EC-20 | PROPOSED until the consolidated PR has every required check green on the exact final candidate and its activation is recorded in the register | PROPOSED | fixes live on the repair branch; consolidated candidate `activation/ec16-ec19-ec20-r19` awaits the post-merge rebase/revalidation |
| EC-18 | PROPOSED until its evidence selection and sequence validation pass on the exact candidate | PROPOSED | parent-prefix + remote-verification layers implemented (AT-13/AT-24) |
| EC-21 | Accepted in principle only; its EC-21 failure on the latest PR events is expected while the declared repair base is not on `main`; do not weaken the check; revalidate after the repair merges | PROPOSED | declared-base job live on the repair branch |
| EC-23 | PROPOSED; the positive-control and broad-exemption blockers are fixed on the repair candidate; activation recorded only when the revised checks pass on the exact candidate | PROPOSED | fingerprint allowlist + verified controls (AT-21/AT-22); all-refs scan rerun on the exact candidate |
| EC-24 | PROPOSED; push-confirm verification now verifies live remote refs | PROPOSED | remote-ref verification layer + RECONCILE records (AT-24) |
| EC-26 | PROPOSED until trigger liveness actually passes on the exact candidate | PROPOSED | `~events=` grammar + live-on trigger validation (AT-23) |
| EC-17 / EC-22 / EC-25 | Prior decision remains in force; no claim that automated enforcement is live on `main` until the reviewed workflow is merged and verified there | ACTIVE (policy) | live jobs run on the repair branch's own pushes/PRs only — NOT on `main` |

Merge path (binding): PR #30 merges only after the critical defects are
fixed, ALL checks are green on the exact final candidate, and the required
independent review is satisfied (GitHub reported no approving review; do not
bypass branch protection or repeat the administrative-merge path). After
that merge, PRs #31–#36 are serially rebased/revalidated against the new
`main` and their complete PR-event workflows rerun (their current EC-21
failures are not passes and must not be waived). PRs #23–#29 close as
superseded only after a file-by-file comparison proves their unique changes
are covered by the reviewed R19 series; their branch refs and audit evidence
are preserved.

### 12.7 Standing items (direction 8) — UNCHANGED

A-01 remains ACCEPTED RISK: the owner rotates the PAT at the owner's
pleasure; no rotation reminder is issued and none is implied here; the token
value appears nowhere. The scientific hard stop is preserved in full: no
deployment, no v25 calibration, no held-out evaluation, no threshold or
cohort changes, no law-semantic changes, no golden-artifact regeneration.
The independent r2 archive verification and the portable 40-row cohort
metadata table remain OPEN and are not closed by this governance work.

## Part 13 — Round-21 records (AUDITOR-ADJUDICATED, 2026-10-20; register rev 10)

### 13.1 The ruling

The auditor ruled: **PR #30 is not approved for merge yet.** Round-20's green
CI results are real but necessary, not sufficient — the independent review
found defects not covered by the acceptance tests: (a) the EC-23 control
checker could disclose the very secret it detects in its failure path
(critical); (b) EC-24 did not bind an unpaired intent to the current
candidate, did not pair confirmations to their own intent identifiers, and
did not validate `expected_remote_tip_before_push` (critical); (c) `HISTORICAL`
did not prove a historical push was real and absent manifest runs were
tolerated as unverified (critical); (d) the protected-main merge transition
was unproven; (e) EC-26 lacked path-filter and invalid-filter-type
validation; (f) status records required reconciliation. This Part records the
adjudication and the Round-21 repair; per the discipline of Part 12.5, the
decision status, policy status, implementation status, exact candidate, CI
status, and main-activation status are recorded separately below.

### 13.2 Repairs implemented on the Round-21 candidate (branch
`repair/round19-mainline-acceptance`, fix-forward from `62fab8a…`)

1. **EC-23 disclosure fix (direction 3):** the control checker prints no
   `Secret`/`Match` values or raw report/log bytes on any path; a requested
   capture log must exist and be readable (fail closed); the HEAD-tree clean
   path requires the report to exist, parse, and carry zero findings — scanner
   exit status alone is never trusted. AT-25 proves the synthetic key never
   reaches stdout/stderr on success or failure paths.
2. **EC-24 binding fix (direction 4):** identifier pairing; the single
   accepted terminal intent must be introduced by the exact candidate under
   test, declare the verified pre-push tip (matching the newest recorded
   branch state; no intervening ref update), and the live remote tip must be
   the exact expected tip; excess unmatched intents FAIL unless voided by an
   `intent-reconcile` record; the terminal state is labeled PENDING, never
   COMPLETE, and no push is manufactured merely to eliminate a correctly
   disclosed terminal PENDING (auditor adjudication). AT-26.
3. **EC-24 historical proof (direction 5):** a HISTORICAL confirmation is
   proven by ancestry from the live tip, or — off the lineage — by a reviewed
   `evidence-back` record whose SHA is a real remote commit (object store or
   the commits API); every run citation anywhere in an entry verifies via the
   manifest, else the Actions API; failed/cancelled runs are disclosed and
   cannot prove a push; a manifest omission fails the gate. The four
   off-lineage legacy confirmations are backed by reviewed records
   P-062..P-065. AT-27.
4. **Protected-main merge transition (direction 6):** the `merge-intent` verb
   pre-declares a reviewed merge (P-061: PR #30, base `f535a71…`, head
   DEFERRED-resolved to the exact final candidate). A push to `main` is
   DECLARED only by a VERIFIED merge-intent: merge-commit parents = declared
   base + approved head; the intent rides the PR branch; the PR is merged
   with an approving non-author review. Raced, unrelated, direct, or
   unapproved updates FAIL. Deterministic AT-28 covers the legitimate merge
   and every failure class. The merge itself is NOT performed by the writer.
5. **EC-26 path coverage (direction 7):** ACTIVE workflow references declare
   `~governs=`; restrictive or invalid path/branch filters FAIL. AT-29.
   Policy status and enforcement status remain recorded separately; a checker
   running on a feature branch activates nothing.
6. **Status and evidence alignment (direction 8):** P0_STATUS top metadata
   corrected and section 10 added; INDEX R19 row reads CLOSE-OUT REJECTED and
   R20 SUPERSEDED-BY-R21 with nothing live on `main`; this Part 13 appended;
   push-log schema v4; the Round-21 reading note replaces the prior one; the
   Round-21 evidence table is generated from live evidence (generator +
   committed output). The PAT-hygiene sweep remains a LOCAL coder claim —
   GitHub CI does not independently verify it; no token value appears in any
   repository file, worklog, or evidence artifact.

### 13.3 Adjudication table (direction 10)

| Rule | Decision status | Policy status | Implementation status | Live on main |
|---|---|---|---|---|
| EC-23 | PROPOSED until the secret-disclosure regression is fixed and the exact-candidate scan controls pass | PROPOSED | disclosure fix + report validation on the Round-21 candidate (AT-25); the exact-candidate CI run is the completing evidence | NO |
| EC-24 | PROPOSED until intent binding, identifier pairing, historical verification, reconciliation validation, and the protected-main transition pass | PROPOSED | all five layers implemented (AT-26/AT-27/AT-28); exact-candidate CI run is the completing evidence | NO |
| EC-26 | PROPOSED until trigger, branch, and relevant path-filter checks pass on the exact candidate | PROPOSED | trigger + branch + path-filter + invalid-type validation implemented (AT-23/AT-29) | NO |
| EC-16 / EC-19 / EC-20 / EC-21 | PROPOSED until their reviewed and rebased activation candidates pass all required checks | PROPOSED | live jobs run on the repair branch's own events; activation candidates await the post-#30 serial revalidation | NO |
| EC-17 / EC-22 / EC-25 | In force as policy under earlier auditor decisions | ACTIVE | live jobs run on the repair branch's own pushes/PRs only | NO (until PR #30 merges and the post-merge verification succeeds) |

### 13.4 Merge path (binding)

PR #30 merges only after: the Round-21 blockers are fixed on the exact final
candidate; the full Stability Gate and Governance Acceptance pass on that
SHA's push and PR events; the post-merge `main` transition is demonstrated
(AT-28 plus the pre-declared P-061); and the required independent approval is
recorded on the PR. The merge is executed only through GitHub's protected PR
mechanism with a MERGE COMMIT. After that merge: PRs #31–#36 revalidate
serially against the resulting `main` (declared bases updated; fresh
exact-head results; their own merge-intents declared before each merges;
EC-21 is not waived); PRs #23–#29 close only after a path-by-path diff proves
their unique changes are covered and their historical refs/evidence are
preserved.

### 13.5 Standing items — UNCHANGED

A-01 remains ACCEPTED RISK: the owner rotates the PAT at the owner's
pleasure; no rotation reminder, no rotation initiation, no credential
operations; the token value is never written into the repository or any
evidence. The scientific hard stop is preserved in full: no production
deployment, no v25 calibration, no held-out evaluation, no threshold or
cohort changes, no law-semantic changes, no golden-artifact regeneration
without their approved gates. The independent r2 evidence-archive
verification, the portable 40-row cohort metadata table, A-02 issuer-side
revocation, and the EC2 deployment secrets remain OPEN and are not closed by
this governance work.

## Part 14 — Round-22 records (AUDITOR-ADJUDICATED, 2026-10-20; register rev 11)

**Status of this Part:** auditor-adjudicated per the binding Round-22 decision;
NOT owner-signed. It explicitly supersedes Round-21's closure CLAIM (not the
historical record): Round-21 is NOT CLOSED, PR #30 is NOT APPROVED FOR MERGE,
and EC-23/EC-24/EC-26 remain PROPOSED and are not live on `main`.

### 14.1 The ruling (transcribed from the auditor's Round-22 decision)

Green CI is necessary but not sufficient; the independent audit found further
fail-open cases in EC-23, EC-24 and EC-26 that the existing acceptance tests
did not cover. `main` remains at `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`;
no merge, no activation of governance rules on `main`, no gate weakening. Work
proceeds on `repair/round19-mainline-acceptance` through ordinary fast-forward
pushes and the reviewed PR path. The mandated repairs: (2) EC-23 — no
report-derived value may reach stdout/stderr; per-element schema validation
failing closed; redaction checks not weakened; AT extension covering
RuleID/File/StartLine/unexpected fields/Secret/Match/non-object entries/
malformed reports with captured stdout+stderr assertions. (3) EC-24 — commit
existence, branch-push provenance, and CI success as separate claims; run
event/head_branch validation; re-audit of P-062..P-065; P-011/P-013 proof
search with honest NOT PROVEN classification when insufficient; machine-
verifiable review/decision references for historical backing records. (4)
EC-24 — branch-scoped, unique, valid intent reconciliation; global
unmatched-intent invariant; expected_remote_tip_before_push bound to trusted
push-event before/after state (all-zero before for new branches; no PR-event
substitution); ABSENT backed by the ref's absent prior state. (5) merge
transition — `pull-requests: read` + live canary; merge_commit_sha, source
ref, merged_at, fresh non-author approval on the exact head (dismissals and
supersessions handled; pagination); pre-merge PR-event verification; A-04
stays an open blocker; DEFERRED supersession fixed; sequential-merge test.
(6) EC-26 — sound conservative glob coverage; ordered negative patterns;
empty/malformed filters fail closed. (7) document/evidence reconciliation
with the append-only push log preserved. (8) two-phase push mechanics with
truthful reclassification, then stop for review. A-01 remains ACCEPTED RISK
(no rotation reminders, no rotation initiation, never reveal or store the
token). No merge of PR #30 until all Round-22 blockers are fixed and
independently reviewed on the exact new head; then live Stability Gate +
Governance Acceptance, functioning pre-merge approval verification, and the
protected-PR merge-commit strategy; after a verified merge, #31–#36
revalidate serially and #23–#29 close only on proven unique-content coverage.
The scientific hard stop is preserved in full.

### 14.2 Repairs delivered on the Round-22 candidate

EC-23: `scripts/ci/secret_control_check.py` type-validates every report
element before examination and never interpolates report-derived values into
diagnostics (counts, booleans, argv-supplied expectations, validated integers
only); malformed reports/rows fail closed with fixed messages; the raw-key
checks cover the report bytes even when the schema fails; broad-exception
paths fail closed with a fixed message (AT-30). EC-24: claim separation
(commit existence / push provenance / CI success) with event + head_branch
validation on every cited run; off-lineage confirmations require push-event
runs, live Events-API PushEvents, adjudicated `push_event=` citations, or
`push-not-proven` reclassifications with machine-verifiable `review_ref=`
references into THIS register (the reference must resolve to a Part whose
section mentions the backed entry — "reviewed" is never merely a string
written by the same record's author); intent reconciliations enforce
preceding eligible targets, branch scope (merge-intent: source_branch),
one-void-per-target, no voiding of consumed intents, and the GLOBAL
unmatched-intent invariant; the terminal intent's declared pre-push state is
bound to the PushEvent before/after state when the events source is
configured (the CI gate always configures it; local tool mode discloses).
Merge transition: `pull-requests: read` + live canaries; merge_commit_sha /
source head.ref / merged_at / fresh non-author approval on the exact head
(commit_id, pre-merge timing, latest-state-per-reviewer) / merge PushEvent
before; the new pre-merge-verification job is the functioning pre-merge gate
(expected red until the independent approval exists — by design, disclosed).
EC-26: conservative witness-based glob coverage with ordered negative
patterns and fail-closed empty/malformed filters (AT-34). Acceptance suite
AT-01..AT-34 green (no check weakened or deleted; AT-27/AT-28 fixtures
strengthened to the new required PR fields). Schema v5 documents every new
rule.

### 14.3 Adjudication table (direction 8; supersedes the Round-21 table's rows for these rules)

| Rule | Decision status | Policy status | Implementation status | Live on main |
|---|---|---|---|---|
| EC-23 | PROPOSED until the report-derived disclosure regression is fixed and the exact-candidate scan controls pass | PROPOSED | disclosure fix + type validation + AT-30 on the Round-22 candidate; the exact-candidate CI run is the completing evidence | NO |
| EC-24 | PROPOSED until claim separation, provenance, reconciliation invariants, pre-push-state binding, and BOTH merge-transition verifications pass | PROPOSED | all layers implemented (AT-31/AT-32/AT-33; schema v5); exact-candidate CI run is the completing evidence | NO |
| EC-26 | PROPOSED until conservative coverage, ordered negative patterns, and malformed-filter rejection pass on the exact candidate | PROPOSED | witness-based coverage + ordered semantics + AT-34; exact-candidate CI run is the completing evidence | NO |
| EC-16 / EC-19 / EC-20 / EC-21 | PROPOSED until their reviewed and rebased activation candidates pass all required checks | PROPOSED | live jobs run on the repair branch's own events; activation candidates await the post-#30 serial revalidation | NO |
| EC-17 / EC-22 / EC-25 | In force as policy under earlier auditor decisions | ACTIVE | live jobs run on the repair branch's own pushes/PRs only | NO (until PR #30 merges and the post-merge verification succeeds) |

### 14.4 Historical evidence-back re-adjudication (direction 3; P-062..P-065)

Search performed 2026-10-20T20:17:51Z with the owner PAT (A-01; the value is never written
anywhere), covering every machine-verifiable source the auditor named:

- **Actions runs on the exact heads** (GET /actions/runs?head_sha=…):
  `17eb88184c9f…` total_count=0; `dab33d946fe3…` total_count=0 — no run of
  any event type exists on either head; the pre-gate-trigger era is
  confirmed by live query, not by trusting the prior record. P-012/P-014's
  cited runs (39408223554 on docs/owner-decisions, 39408297643 on
  governance/canonical-docs) DO exist as push-event runs with matching
  head_branch and head_sha; their conclusion=failure is disclosed separately
  and never counted as successful gate evidence.
- **Events API** (full retained window: 300 events,
  2026-10-18T16:00:59Z..2026-10-20T19:13:33Z, 20 PushEvents): no PushEvent
  with after == `17eb88184c9f…` and none with after == `dab33d946fe3…`;
  exactly ONE retained PushEvent on refs/heads/proposal/ec16-ec26-enforcement: id 24612279632, created
  2026-10-20T10:32:21Z, before = `dab33d946fe3…` (after/size elided by the
  API) — a later push's own record that the branch's prior tip was
  `dab33d946fe3…`, consistent with P-013-CONFIRM's ls-remote at 10:28:58Z.
- **Incident and retained-Git evidence**: the P-013-CONFIRM body's
  "fast-forward 17eb881..dab33d9" narrative and the durability bundles prove
  commit existence and describe the transition, but a log narrative is not
  independent machine-verifiable ref-transition evidence (transcript-only
  proof is rejected per direction 3.6).

Adjudication (this Part is the canonical record the review_ref= fields of
P-074/P-075 resolve to):

- **P-011-CONFIRM (proposal/ec16-ec26-enforcement @ `17eb88184c9f…`) →
  NOT PROVEN.** Reclassified by P-074 (push-not-proven): the push's own
  event (~09:34Z per the log) rolled off the retained window, no Actions
  run exists on the head, and no other ref-transition evidence survives.
  The commit exists (P-062's commits-API HTTP 200); commit existence is not
  push provenance. The original entries stay unchanged; this push is never
  counted as evidence for any declaration.
- **P-013-CONFIRM (proposal/ec16-ec26-enforcement @ `dab33d946fe3…`) →
  PUSH-PROVEN by prior-tip citation.** Backed by P-075 (evidence-back,
  review_ref=this Part): PushEvent 24612279632 (2026-10-20T10:32:21Z, same
  ref) records its before as `dab33d946fe3…` — the ref provably pointed
  there (refs move only by pushes), establishing the branch update that
  P-013-CONFIRM records. The event is re-verifiable live while retained;
  the citation is recorded in P-075 in the schema-v5 prior_tip form.
- **P-012-CONFIRM and P-014-CONFIRM → PUSH-PROVEN by push-event runs**
  (direction 3.3, verified live again this cycle): 39408223554 (event=push,
  head_branch=docs/owner-decisions, head_sha=`78b147d5…`,
  conclusion=failure) and 39408297643 (event=push,
  head_branch=governance/canonical-docs, head_sha=`b12f30e9…`,
  conclusion=failure). The push event on the exact branch and head proves
  the ref pointed there; the FAILED conclusions are disclosed separately
  and are not counted as successful gate evidence.

### 14.5 Merge path (binding; unchanged in substance from Part 13.4, with the Round-22 additions)

PR #30 merges only after: all Round-22 blockers are fixed on the exact final
candidate; the full Stability Gate and Governance Acceptance pass on that
SHA's push and PR events; the pre-merge-verification job passes WITH the
required independent, non-author, current-head approval recorded on the PR
(until then that job is expected red — the gate functioning); the merge
executes only through GitHub's protected PR mechanism with a MERGE COMMIT;
the post-merge `main` run verifies the merge-intent (parents,
merge_commit_sha, source ref, merged_at, fresh approval, PushEvent before).
After that merge: PRs #31–#36 revalidate serially against the resulting
`main` (declared bases updated; fresh exact-head results; their own
merge-intents declared before each merges; EC-21 is not waived); PRs #23–#29
close only after a path-by-path diff proves their unique changes are covered
and their historical refs/evidence are preserved.

### 14.6 Standing items — UNCHANGED

A-01 remains ACCEPTED RISK: the owner rotates the PAT at the owner's
pleasure; no rotation reminder, no rotation initiation, no credential
operations; the token value is never written into the repository or any
evidence. A-04 (exact branch-protection required-check configuration) is an
OPEN blocker until independently verified; protection is never bypassed. The
scientific hard stop is preserved in full: no production deployment, no v25
calibration, no held-out evaluation, no threshold or cohort changes, no
law-semantic changes, no golden-artifact regeneration without their approved
gates. The independent r2 evidence-archive verification, the portable 40-row
cohort metadata table, A-02 issuer-side revocation, and the EC2 deployment
secrets remain OPEN and are not closed by this governance work.

## Part 18 - Round-24 records (AUDITOR-ADJUDICATED, 2026-10-21; register rev 13)

### 18.1 The decision

The governance-only workstream is STOPPED and the unfinished LAW-160 engine
correction resumes immediately. The next round is not spent on the P-122
reconciliation, archive tags, pull-request pruning or the O-1 approval process:
those activities do not fix the outstanding engine or calibration defects.
PR #37 stays open and unmerged; the engine patch does not go into it. A new,
isolated feature branch is cut from freshly verified `main`
(`63da0e347002ed9b77071ff1ed0cefda7c0ce88a`); the LAW-160 patch goes there.

### 18.2 The defect and the correction (LAW-160, Chain Integrity)

`get_sequential_residues()` returns the residues that are PRESENT (sorted,
de-duplicated); it does not guarantee that two successive entries are
consecutive in the biological sequence. LAW-160 measured the C-alpha distance
between successive entries, so a chain with missing intervening residues (a
modelling gap) was reported as a physical chain break - a FALSE VETO.

The correction factors the adjacency predicate into one source of truth
(`Tier1Measurements._adjacent`), keeps the semantics of the existing
`_is_sequential()` VERBATIM, adds the tuple-form
`_is_sequential_residues()` for the `(res_seq, insertion_code, res_name)`
tuples, and measures only sequence-adjacent pairs. `get_sequential_residues()`
is unchanged (other laws use it). The canonical 4.5 A threshold, the law ID,
the operator, the units and the classification are unchanged; the canon and the
modality matrix are untouched.

**Verdict impact (explicit):** correcting pair eligibility can remove false
LAW-160 vetoes and therefore change the verdict of affected structures. A
structure that previously vetoed on a modelling gap now passes, or - when no
pair is eligible at all - is reported as unmeasured (18.3). Pre- and
post-correction results are not comparable for chains with gaps.

### 18.3 OPEN: the zero-eligible-pair representation (release blocker)

With no sequence-adjacent C-alpha pair, LAW-160 is not measured. The previous
code emitted `observed=0.0`, `PASS`, `sample=0` - a fabricated value and a PASS
for a law that was never evaluated. The patch removes that: the engine emits
no measurement, and the adjudicator's EXISTING missing-measurement default
applies (status FAIL, sample 0 - the representation every other unmeasured law
already receives). Nothing is fabricated, no status is invented, no PASS is
claimed.

That is not the final answer: "FAIL" for an unmeasurable law still reads as
a scientific failure. A dedicated unmeasured state needs an owner decision on
the approved status vocabulary (the engine's vocabulary today is exactly
PASS / VETO / FAIL) and on the response schema
(`backend/tos/schemas/response_v1.py`: `PhysicsMeasurement.observed` is a
REQUIRED float, so a null observation is not available without a schema
change). The proposal, the inspected engine/adjudicator/schema/consumers and
the two options are recorded in
`docs/audit-evidence/2026-10-21/LAW-160-R24-EVIDENCE.md`, section 5. The task is
reported as PARTIALLY COMPLETE with this item explicitly outstanding.

### 18.4 Evidence

`docs/audit-evidence/2026-10-21/LAW-160-R24-EVIDENCE.md` records the base and
candidate SHAs, the runtime versions, every command with its exit code and
output, the test counts against the unmodified baseline in the same
environment, the canon/matrix hash results, the secret scan, the changed-file
digests, the verdict impact, the zero-sample proposal and the unresolved
limitations. The LAW-160 target tests and the full backend suite were run on
the exact candidate.

### 18.5 Standing items

The scientific/production hard stop is preserved in full: no merge to `main`,
no deployment, no v25 run, no held-out evaluation, no threshold or cohort
change, no golden-artifact regeneration. PR #37 is not merged and the new
feature branch is not merged. A-01 remains ACCEPTED RISK. LAW-130 (the
fail-open clash measurement and the `<` where the canon specifies `<=` at the
20.0 boundary) and the calibration runner (the missing `X-API-Key` header,
continuing after failures, and writing partial results) are the NEXT engine
and calibration corrections, each in its own separate branch.

## Part 19 - Round-24 NOT_MEASURED contract and its enforcement (AUDITOR-ADJUDICATED, 2026-10-21; register rev 14)

### 19.1 The approved contract (recorded, not re-decided)

The auditor approved the explicit unmeasurable state for LAW-160 (the Round-24 Option A
contract, as adjudicated and recorded in
`docs/audit-evidence/2026-10-21/LAW-160-R24-EVIDENCE.md` and in the pull request body of the
canonical branch). The contract is:

- An unmeasured law is reported with `status = "NOT_MEASURED"`, `measured = false`,
  `observed = null`, `deviation = null` and `sample_size = 0`. A measurement that could not be
  taken is never presented as a numerical zero, and never as a scientific failure.
- A required deterministic law that is NOT_MEASURED yields the verdict `INDETERMINATE` - it
  neither passes the structure nor fails it scientifically. A known deterministic veto takes
  precedence and remains a VETO.
- The response model enforces the cross-field consistency of a measurement row: a measured row
  carries a numeric observation, an unmeasured row carries the NOT_MEASURED status with null
  observation, null deviation and a zero sample. A contradictory row is rejected.
- The PDF report states the unmeasured laws explicitly, in their own section, independently of
  whether any other law was violated; an unmeasured law is not counted as a violation, is never
  named as a VETO cause, is never a refinement target, and is never rendered as a bare `None`.

Part 18 is preserved unchanged; this part records the contract and the state of its enforcement
on the canonical branch `feature/law-160-sequence-adjacency` (pull request #40).

### 19.2 State of the implementation on the canonical branch

- The engine emits the NOT_MEASURED row for LAW-160 when no eligible C-alpha pair exists
  (`backend/tos/engine/tier1_measurements.py`).
- The adjudicator maps a required deterministic NOT_MEASURED law to INDETERMINATE and keeps a
  real deterministic veto as a VETO (`backend/tos/engine/adjudicator.py`).
- The response model validates measurement consistency (`backend/tos/schemas/response_v1.py`).
- The PDF reporting paths render the unmeasured section, exclude unmeasured laws from
  violations, VETO causes and refinement recommendations, and use the shared observation
  rendering helper (`backend/tos/forensic_artifacts/pdf_generator.py`).
- The end-to-end tests are rendering-level: the generated PDF is decoded and its content is
  asserted (`backend/tests/test_not_measured_pdf_contract.py`), and the regressions of the
  superseded R24 candidate are preserved in `backend/tests/test_law_160.py`.

### 19.3 Standing items - UNCHANGED

The scientific/production hard stop is preserved in full: no merge to `main`, no deployment, no
v25 calibration run, no held-out evaluation, no threshold or cohort change, no golden-artifact
regeneration. A-01 remains ACCEPTED RISK. LAW-130 and the calibration runner remain the NEXT
separate corrections, each in its own branch, and are not touched by this record.
