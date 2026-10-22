# Port Rationale: Review Candidate `review/od-scope-v25-2026-10-22`

Candidate head: `f962c6a` (base: `63da0e3` = verified `main`).
Commits: `818717a` (duplicate-tree removal + `tier1_laws.py` cleanup),
`ff22b72` (ported reviewed round-25 changes + OD-04 fixes), `f962c6a`
(LAW-182 per-chain correction found during C5 regeneration).
Delivery record compared against: `engine/od-scope-v25` @ `97d62a8` (frozen; merge-base `0eae6d5`).
Written under audit round-26 directives 1-5. Read-only references: brief
`052a84ca…`, OWNER_DECISIONS `12906bf5…`, constitution snippets verified earlier
(constitution `211ad2c5…`, ENGINEERING_CONSTITUTION `b6bb35e8…`, P0_STATUS `a31517ab…`).

## 1. Why a new isolated re-based branch

The round-25 delivery branch shares only its first `23` commits with `main`
(`23` ahead / `45` behind). It cannot be merged: it carries `control-atlas-saas/`
and `pipeline/` trees, modularized engine files (`constitution_checks.py`,
`station.py`, etc.), and evidence hashes that do not exist on `main` and were not
reviewed as reviewable here. Rebuilding a candidate directly on `main` makes the
diff auditable file-by-file and keeps the frozen branch untouched.

## 2. Branch push correction (disclosed)

The first push to the new remote branch name was sent with HEAD (`97d62a8`) by
mistake (single refspec push with no explicit ref), moments after branch creation.
It was corrected seconds later by force-push of the initiating review branch to
`63da0e3` from the same workspace. The frozen delivery branch `engine/od-scope-v25`
and `main` were never force-pushed. No other users consumed the branch in that
window. Recorded on `+x` force-push transparency for any session.

## 3. What was ported (and why)

| Area | On the delivered branch | Candidate port | Rationale |
| --- | --- | --- | --- |
| `governance/law_registry.yaml` (4 hunks: LAW-105/130/155/200) | present | ported | Direct reviewed registry changes: scope→core, status_enforced measures text corrections |
| `governance/law_registry.sha256` | re-pinned (`7e3b42d5`) | re-pinned (`069a6c99…`) | bytes differ due to documented deviations below; pin follows registry |
| LAW-105 dual-basis row (residue vs atom coverage) | present | ported | the delivered hunk was verified line-by-line |
| LAW-110 measures sentence | NOT on the delivered branch (commit message overclaimed a doc update) | candidate adds it explicitly | registry↔code agreement: LAW-110 is shipped core-scoped + pair counting; the deviation is this candidate's own, disclosed in §5 |
| semantic LAW-155/182/200 relabeling (FAIL (Advisory), no deterministic VETO) | present | ported + completed | delivered branch message claimed a "LAW-182 advisory-only" registry edit that does not exist in its diff; candidate closes that gap only for code-side status text, no registry edit for LAW-182 (see §4) |
| `adjudicator.py` 1-line advisory normalization (`advisory_experimental` → `heuristic`) | present | ported | single reviewed line |
| INDETERMINATE cause blocking in `_compute_verdict` | delivered branch structure re-cast | candidate implements via `INDETERMINATE_CAUSE_STATUSES` + `indeterminate_causes` field | implements the same behavior without the branch's refactor surface |
| honest empty-set statuses (NOT_MEASURED / NOT_APPLICABLE / MEASUREMENT_FAILED) for changed laws | partial | ported, in the same spirit (`tier1_measurements.py` changes only) | enables fail-closed behavior with main's single-file structure |
| scope text in `station_sop.py` (3 lines LAW-130/155/200) | present | ported | threshold/operator values unchanged; scope text matches the registry |
| 50-70 band removal in forensic artifacts (pdf/viz/3D/heatmap) | partial (PDF hunks) | ported + completed | candidate also removes the fabricated proportion bars (§3b) with a measured replacement |
| LAW-130 `<`→`<=` boundary | NOT on the delivered branch (`clash_score < threshold` at `97d62a8`) | candidate fixes with tests | audit defect list; canonical operator `<=` + exact-boundary semantics |
| self-comparison to the reported value | n/a | candidate compares the ROUNDED observed value | keeps a single deterministic number visible to auditors; SEE change detail in commit ff22b72 |

## 3b. Confidence bar: measured bins vs coverage-derived fractions

The delivered branch still rendered bin proportions from coverage attributes
(`confidence_viz.py` derived `very_high`/`high`/`low` fractions from
`coverage_pct` and its multipliers, and the PDF page duplicated the old video
even before band removal). On the candidate:

- `pdf_generator.py` extracts per-residue pLDDT values once from
  `payload["pdb_b64"]` (`_extract_ca_atoms`) and passes measured
  `bin_counts = {very_high: n>=90, high: 70<=n<90, low: n<70}` to
  `generate_confidence_bar`.
- `confidence_viz.py` renders proportions from the measured counts; when no
  observations exist it renders an explicitly-labelled atom-coverage bar via
  `_atom_coverage_indicator(coverage_pct)` (no synthetic bins).
- `residue_heatmap.py` keeps 4 bands → 3 bands; `structure_viz.py`
  legend and 3D coloring likewise; the PDF text labels drop the 50-70 language.

## 4. What was deliberately NOT ported

| Omitted item | Reason |
| --- | --- |
| `control-atlas-saas/` nested duplicate tree (84 files) | exists only on the delivery branch; duplication of the repo itself causes confusion in context and metadata overlap. Deleted on the candidate in commit `818717a`. |
| `pipeline/` (parallel pipeline/ v2 integration) | not reviewed here; the delivery branch's parallel pipeline stack is the OD consumer, out of scope for the LAW-level candidate |
| modular engine files (`constitution_checks.py` etc.) | not reviewed as part of the round-25 package; main retains the single-file structure |
| `governance/law_semantic.sha256` + `compute_law_semantic_hash` | not reviewed here; the hash binding artifact is out of scope for this candidate |
| WP-SCI-1 measurement-definition package (P0 items 1-13 + field-order residue pool) | not reviewed as part of the round-25 delta (`c78ef15..97d62a8` carries none of it; it pre-dates the round-25 review window) and stays on main's pre-P0 pool semantics; the C5 cross-lineage section documents exactly which symptoms remain (1PPE, 5T4O) |
| `detect_confidence_source` / decompose methods | main does not need them (static confidence profile parsing); not part of the reviewed defect list |
| LAW-182 registry row edit (advisory-only statement) | OD-20 owner decision remains the canonical NMR-separation source; the clarification file is NOT a register amendment (see conflict record) |
| LAW-110 `sample` = `pairs_checked` | delivered LAw-110 keeps `sample: peptide_pairs_checked` (WP-SCI-1, already at c78ef15) | candidate matches that field name/semantics; the remaining deviations are the two listed in §5-2 | the deviation is this candidate's own, disclosed in §5-2 |

## 5. Deviations from the delivered diff (documented, no hidden divergence)

1. **Registry content.** The candidate's registry rows keep the language
   matching main's existing format, avoid the delivered branch's
   `status_semantics` YAML block, and re-pin the flat sha256 file. Hash value:
   `069a6c993d27ffd8450960e734445f6810c2a0d66cf7042357399bb589a9600b`.
2. **LAW-110 pairing basis.** The delivered branch keeps the WP-SCI-1
   `sample = peptide_pairs_checked` field, but walks the FULL file-order
   residue pool and skips non-core endpoints by identity (`core_ids`) —
   preserving file-order adjacency. The candidate instead iterates the
   core-filtered (sorted) pool directly and pairs consecutive core
   members, so for PREDICTED structures the pair set coincides only when
   core residues are also pool-adjacent. The registry `measures` text
   documents the candidate's semantics; disclosed here as a
   candidate-specific deviation (exposed by no verdict change: C5
   attributed comparison stays 0/48).
3. **Adjudicator.** Same reviewed behavior (heuristic advisory normalization +
   INDETERMINATE causes for missing failed evidence), implemented as a
   two-pass `_compute_verdict` (VETO first, then missing-evidence causes) in
   main's single-file structure instead of the delivery branch's separate
   conditions module. Behavior per adjudicated arbitration is pinned in
   `backend/tests/test_adjudicator.py` (added tests).
4. **Atom coverage.** The candidate exposes the atom basis of LAW-105 coverage
   in the same return tuple (no CLI/protocol surface changed); delivered
   branch exposed a `computed` field. Consumers: `backend/main.py` (`_run_physics_sync`)
   and the PDF layer (measured bins). The contract stays a dict with the same
   keys, `coverage` unchanged, `atom_coverage` ADDITIVE (new key).
5. **No offer of LAW-182 in the registry.** Canon stays authoritative from
   OWNER_DECISIONS OD-20. The candidate documents the conflict record separately
   (see `conflict-od-20-clarification.md`).

## 6. Verification

- Full backend suite at committed candidate head `f962c6a`: **399 passed, 5 skipped** (0 failed).
- New pin suites: `backend/tests/test_plddt_scope_pin.py` (5) and
  `backend/tests/test_law_130_boundary.py` (9, incl. exact-20.0 boundary and
  tooling-failure fail-closed adjudication).
- Governance pin updated coherently: `governance/law_registry.sha256`
  = `069a6c99…` == test pin in `backend/tests/test_law_registry.py`.
- Canon/matrix literals in `toscanini-gate.yml` untouched (`6a9cd4b4349b81de…`,
  `ed1731ac…`).
- OD-16 hard stop honored: no merges/deploys/calibration/threshold/cohort/golden
  changes; candidate is a review branch only.

## 7. C5 evidence regeneration (round-26 directive)

`docs/audit/2026-10-22_od_port_candidate/scripts/regenerate_c5_evidence.py`
materializes each engine byte-exactly with `git archive`, re-fetches the v24
pinned cohort (SHA-256 + byte-count verified against the frozen manifests
`9722c3e45d31…` / `f602351f4844…`), runs isolated per-(engine,target)
evaluations, and writes portable artifacts (no machine paths) to
`evidence/engine/` mirroring the delivered `C5` layout:

- Three engines: `c78ef15` (od-scope lineage, reproduces the delivered
  old-engine records: e.g. 1AKE atoms 3312, det 9/10, identical LAW-182
  observed values), `63da0e3` (candidate's own base, no port), current head.
- **Primary attributed comparison (base vs candidate): 0/40 v24-pinned and
  0/9 supplementary verdict changes; 23 law rows change status, all direct
  outcomes of the reviewed delta** — 17 LAW-182 VETO→FAIL (Advisory)
  relabels (observed + sample values IDENTICAL in every row; the
  three-state/per-chain semantics this depends on are re-certified by the
  reviewed commit's own `test_law_182.py`), 4 LAW-200 + 1 LAW-155-style
  core-scope/advisory relabels, and 1I0T's LAW-110 NOT_MEASURED /
  LAW-182 NOT_APPLICABLE for an empty pLDDT core (honest empty states,
  scope pin: LAw-110/130/155/200).
- A first regeneration attempt exposed a translational defect: the ported
  LAw-182 had collapsed c78's per-chain loop into one global centroid list
  (min-across-chains lost). `f962c6a` restores the delivered per-chain /
  worst-chain semantics; the recorded narrative is from the fixed head.
- Cross-lineage documentary comparison (c78ef15 vs candidate): 2 verdict
  changes (1PPE, 5T4O) + 72 law-row changes, all attributable to the P0/
  WP-SCI-1 measurement-definition package that the delivery branch carries
  BEFORE the round-25 window (`c78ef15` already contains it) and the
  candidate does not (explicit §4 omission). 1PPE is the documented P0
  item-4 symptom (sorted residue pool re-pairs chain E around
  183/184/184A/185 and fabricates LAW-110 backbone gaps); 5T4O stems from
  main-lineage LAW-150 empty-sample / LAW-195 not-applicable semantics.
  Neither is an effect of the reviewed round-25 port.
- Supplementary AFDB acquisition reproduces the delivered record parity:
  9/10 acquired; `AF-P62988-F1` unavailable (404, as delivered);
  `AF-P0DTC2-F1` resolves the same migrated entity id
  (`AF-0000000365840314`, model_v1) with the same SHA-256 as the delivered
  artifact. Supplementary cohort: 5 PASS / 0 VETO / 4 INDETERMINATE, 0
  verdict changes on both comparisons.
