# OWNER_DECISIONS — Consolidated Decision Register (rev 4)

**Provenance.** This is the single consolidated owner-decision register
adopted per the Round-11 direction (U5). It supersedes every earlier
register draft (rev 1 in the r1 archive, rev 3 inside the r2 evidence
archive), which remain sealed, byte-identical historical records. The
OD-01…OD-21 list below is the auditor's consolidated table, reproduced
verbatim in its Decision and Blocks columns. Internal D-xx labels are
retired from all user-facing text and live only in
`INTERNAL_CROSSWALK.md` (marked internal). No OD number is ever reused.

**Status discipline (binding, U6).**

- Every decision below is **OPEN** as of 2026-10-20. A decision may be
  recorded as **DECIDED** only with (a) the owner's written decision text,
  (b) the owner's name, and (c) a date. Until all three exist, the status
  stays OPEN and nothing is treated as ruled.
- **Recommendations are advisory only** and live exclusively in the
  "recommendation" field — never in the status field, never phrased as a
  ruling. A recommendation by the coder confers no authority.
- Decision status and implementation status are separate. Where
  implementation already exists on a branch (e.g. OD-02's precedence rule
  since WP-SCI-1), that fact is recorded as context; it is not an owner
  approval and does not close the decision.

---

## Part 1 — Consolidated register: OD-01 … OD-21 (all OPEN)

| ID | Decision (auditor's consolidated table, verbatim) | Blocks | Status | Recommendation (advisory only) |
|---|---|---|---|---|
| OD-01 | Canonical governance file: AGENTS.md pointer, ENGINEERING_CONSTITUTION.md, or both | Docs merge | **OPEN** | none — layout choice is the owner's |
| OD-02 | Verdict precedence: measured VETO over INDETERMINATE. Implemented since WP-SCI-1; gate open | Verdict-table approval | **OPEN** | ratify the as-implemented rule (see context §2.1; rule with the corrected 0-of-40 figure, not the withdrawn 20-of-40 bound) |
| OD-03 | pLDDT scope for LAW-130, 155, 200 (LAW-200 flips VETO→PASS: 477 vs 68 Å³/atom) | Threshold work | **OPEN** | none — scientific threshold call |
| OD-04 | Disordered-protein policy for LAW-200 | Merged with OD-03 unless owner splits it | **OPEN** | decide together with OD-03 unless a split is explicitly wanted |
| OD-05 | LAW-190, 210, 220, 230: implement or retire | Registry | **OPEN** | retire-or-implement per entry; do not leave declared-but-unenforced entries visible as enforced (EC §4 registry discipline) |
| OD-06 | Certificate: keep with reportlab, or delete | Dependency cleanup | **OPEN** | tie to the OD-09 scope cut: if certificates survive v1, keep reportlab pinned and tested; if not, delete both together |
| OD-07 | Audit-ID policy: identical structures share an ID | ID schema | **OPEN** | none — ID semantics call |
| OD-08 | History access: key-scoped or admin-only | /history routes | **OPEN** | none — security posture call |
| OD-09 | v1 scope cut | Feature removal | **OPEN** | cut before threshold work so calibration scope is stable |
| OD-10 | Cohort numerics: resolution, R-free, completeness, split fraction | Calibration, thresholds | **OPEN** | derive only from the post-reconciliation cohort (see OD-17–OD-19) |
| OD-11 | P0DTC2 disposition: migrated ID `AF-0000000365840314` | Cohort membership | **OPEN** | none — cohort membership call |
| OD-12 | P62988 disposition: UniProt Inactive, no AF entry | Cohort membership | **OPEN** | none — cohort membership call |
| OD-13 | AlphaFold version pin for calibration | v25 run | **OPEN** | pin before any v25 run; the v25 run cannot start without this ruling |
| OD-14 | V4 approval of the verdict-change table, including 1PPE VETO→PASS | Any verdict-changing merge | **OPEN** | none — scientific approval call; see context §2.2 for the corrected pass-eligibility figure |
| OD-15 | Disposition of the two pushes (`science/p0-correctness`, `chore/hygiene`) | Branch retention | **OPEN** | none — see record §3.1 |
| OD-16 | Merge and deploy sequence | Any merge to main | **OPEN** | none — sequencing call; hard stop remains in force until ruled |
| OD-17 | Metadata integrity (Round-9 item) | Metadata work | **OPEN** | approve the quarantine-then-rebuild path (corrected-cohort proposal v2, classes A–D) rather than in-place silent correction |
| OD-18 | Cohort eligibility (Round-9 item) | Cohort | **OPEN** | exclude 1I0T (DNA-only) and 1DDL (protein+RNA) from the protein-only core cohort; keep both in an explicitly classified audit set |
| OD-19 | Stratified population (Round-9 item) | Cohort, thresholds | **OPEN** | prespecify strata and sample-size feasibility only after identity/modality/resolution reconciliation |
| OD-20 | NMR subclass (Round-9 item) | Modality rules | **OPEN** | keep NMR out of pooled resolution strata until model selection and ensemble treatment are prespecified |
| OD-21 | Word target: confirm the number and its source requirement | Document length | **OPEN** | no number chosen; the figures 5,460 / 5,470 have no record in any tracked file on any audited branch — keep the target "unconfirmed" and document lengths unchanged until the owner supplies the number and its source |

## Part 2 — Required context attachments

### 2.1 OD-02 (verdict precedence) — pass-eligibility correction (U7)

The figure **"at most 20 of 40 could pass"** that appeared in earlier
rounds was the auditor's upper **bound**, derived without recomputation,
and has been withdrawn. The tracked artifact
`docs/audit-evidence/2026-10-19/science-validity/independent_recompute_20_of_40.json`
(blob `6cda70e2791d50325643e8fbbf56db7ee9948f8d`, on
`science/p0-correctness-2026-10-19` @ `ebc1fd79…`) establishes the
corrected figure for the **recorded v23 run**:

- `experimental_total = 40`, `could_pass = 0`
- 38 blocked by LAW-120 at the 10° threshold
- 2 blocked by LAW-170: `1I0T`, `1DDL`
- recorded verdict histogram (50-target run): VETO 44 / INDETERMINATE 5 /
  PASS 1

Both numbers describe the **recorded v23 run** — the 20-of-40 bound and
the 0-of-40 recompute are historical statements about that run, not
predictions about a post-fix run. The owner should rule on verdict
precedence with the 0-of-40 figure in view.

### 2.2 OD-14 (V4 approval of the verdict-change table) — same correction (U7)

The verdict-change table awaiting V4 approval (including the 1PPE
VETO→PASS change) was produced under the same recorded v23 run. The same
corrected figure applies: under the recorded rows, 0 of 40 experimental
targets could pass (38 LAW-120@10°, 2 LAW-170). Approval of the table and
of any verdict-changing merge should proceed with the corrected figure,
not the withdrawn bound. The v24 calibration artifacts are marked
superseded/historical (`backend/calibration/runs/v24/calibration_v24_SUPERSEDED.md`).

### 2.3 OD-17 through OD-21 (Round-9 items) — recorded substance

These five are **directives recorded for owner decision, not decisions**.
Decided-by and decided-date are blank for each; none is implemented.

- **OD-17 (metadata integrity).** Quarantine the historical v23 cohort
  metadata from calibration and threshold derivation until a corrected,
  source-backed cohort proposal is reviewed. Preserve
  `backend/calibration/dataset_v23.py`
  (blob `3a5358374b2037b9475bd2d494754944485da5ac`) and historical
  calibration outputs unchanged. Row-derived audit (r2 evidence archive,
  `02-`): 22/40 descriptive errors (21 identity + 1 state), 5/40 method,
  24/40 strict resolution (23/40 at 0.01 Å), 2/40 non-protein-only,
  0 parse failures, 0/40 coordinate-byte mismatches vs the v24 manifest —
  32/40 rows carry at least one metadata defect. Byte-match establishes
  acquisition stability only, never label validity.
- **OD-18 (cohort eligibility).** 1I0T is DNA-only (Z-DNA, 0.60 Å);
  1DDL is protein+RNA. Both are excluded from the protein-only core
  cohort absent explicit protocol support, and both are exactly the two
  LAW-170-blocked structures in the artifact cited in §2.1.
- **OD-19 (stratified population).** Recorded strata are built from
  defective labels (post-correction, the cryo-EM stratum loses 3 of 8
  entries to X-ray/ssNMR). Recompute strata only after reconciliation;
  prespecify sample-size feasibility before outcomes are reviewed.
- **OD-20 (NMR subclass).** 1RVS and 2N0A are solid-state NMR recorded
  under a solution-NMR label. Keep NMR separate from resolution-based
  strata until model selection, ensemble treatment, and applicable laws
  are prespecified.
- **OD-21 (word target).** The transcript figures 5,460 / 5,470 have no
  record in any tracked file on any audited branch (auditor search and
  independent coder `git grep` across `origin/main`,
  `origin/science/p0-correctness-2026-10-19`,
  `origin/integration/p0-hardening-2026-10-19`,
  `origin/policy/engineering-constitution` agree). No number is chosen;
  the target stays "unconfirmed"; **document lengths are unchanged**
  (U8). See `WORD_TARGET_OD-21_UNCONFIRMED.md` in the r2 evidence
  archive.

## Part 3 — Repository and delivery state relevant to rulings

### 3.1 Not-pushed and unrecoverable branches (updated this round)

- **Docs branch, tip `94b0f7c`** (and the earlier `df44aa6`): never
  pushed; the commits are **not present in the current environment**
  (lost with the prior sandbox reset — recorded in the lost-work
  register, r2 evidence archive `08-`). Bundles or patch series for
  these exact tips **cannot be produced from available objects**, and
  none were fabricated. The register content itself survives and is
  delivered as this file. GitHub verification: both tips absent from all
  remote refs.
- **N8 work branch, tip `e9985c6`** (per-residue altloc fix, first of a
  4-commit series reported as 506/5 → 545/5): never pushed; commits
  lost with the sandbox reset. Same non-fabrication rule applies. A
  **fresh** per-residue altloc feature patch is delivered this round
  under U9 (`patches-n8-per-residue-altlocs/`) — it is new work on the
  science branch tip, explicitly NOT a restoration of `e9985c6`, and its
  commit SHAs differ.
- **Pushes performed this round (owner-authorized):** exactly one — the
  evidence-delivery branch `evidence/round9-r2-2026-10-20` carrying the
  `outputs/` folder (this register, the r2 evidence archive, verifier,
  SHA256SUMS, and the U3/U9 materials). No existing branch, tag, or ref
  was moved; nothing was merged or deployed; PR #16 remains
  open/unmerged.
- **Stale legacy branch** `fix/c2-api-key-enforcement` (February 2026,
  tip `2d4de2f`, 0 ahead / 68 behind `origin/main`): untouched;
  keep-or-delete is owner action A-06.

### 3.2 Standing hard stop

Until the owner rules (OD-15, OD-16): no merges to `main`, no
deployments, no v25 run, no held-out evaluation, no threshold or cohort
changes, no golden-artifact changes. The only authorized activities are
the U9 unblocked items (per-residue altlocs as a patch, the MODRES
census as documentation, runner provenance as a patch with a
demonstration manifest) and the evidence delivery itself.

## Part 4 — Owner actions A-01 … A-06 (tracked separately)

| ID | Action | Status |
|---|---|---|
| A-01 | Revoke and reissue the PAT (credential exposed in transcript history) | open — owner-executed; no revocation evidence to date |
| A-02 | Rotate the Gemini API key still present in `main` HEAD tree | open — owner-executed |
| A-03 | Set EC2 production secrets (`TOSCANINI_API_KEYS`, `INTERNAL_CALLBACK_SECRET`, `ALLOWED_ORIGINS`, `GEMINI_API_KEY`) | open — owner-executed |
| A-04 | Branch protection on `main` requiring the Toscanini Stability Gate; disable direct pushes | open — owner-executed |
| A-05 | Merge the constitution (canonical governance per OD-01) | open — owner-gated |
| A-06 | Decide keep-or-delete for `fix/c2-api-key-enforcement` (February 2026 branch, diverged; coder-measured 0 ahead / 68 behind `origin/main`) | open — recorded |

## Part 5 — Closure protocol

When the owner answers an OD in writing:

1. Create `docs/decisions/OD-<id>-<date>.md` carrying the owner's
   written decision text, name, and date, and its exact scope.
2. Update this register: set decided-by, decided-date, decision status
   (DECIDED / superseded), and link the decision document.
3. Implement **only DECIDED items**; nothing is implemented on the
   strength of a directive, a recording, or a coder recommendation.
4. Run the CI sequence exactly as the Engineering Constitution requires
   (checks on the exact merge candidate; EC §7), never weakened for
   green.
5. Write evidence under `evidence/<item>/` with exact commands, SHAs,
   and counts (EC §6, §14).
6. Update `docs/P0_STATUS.md` for the affected items.
7. Hand the complete bundle (patches, logs, test results, decision
   record) to the auditor for independent verification.

## Register-level notes

1. **One register.** This file is the only current register. Earlier
   drafts inside the sealed r1/r2 evidence archives are historical
   records and are not amended.
2. **No silent status conversions.** No owner-blocked item is upgraded;
   recording context or a recommendation never asserts owner approval.
3. **Decision ≠ implementation**, always. Existing implementation (where
   noted) rides on unmerged branches; nothing here authorizes release
   activity while the hard stop is in force.
4. **Counting note.** Where this register cites file sizes or counts,
   the single pinned convention applies: bytes (octets), UTF-8 chars,
   whitespace-split words, and LF / split-elements / non-blank lines
   reported separately, under `LC_ALL=C.UTF-8`.
