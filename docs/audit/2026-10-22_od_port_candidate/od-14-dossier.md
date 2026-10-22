# OD-14 Reconciliation Dossier — Read-Only Supporting Record

Candidate branch: `review/od-scope-v25-2026-10-22` @ `f962c6a`
(base `63da0e3` = verified `main`). Frozen delivery branch:
`engine/od-scope-v25` @ `97d62a8`.

This is a **read-only reconciliation document**. It records how the
candidate's work under audit round-26 stays inside the auditor decisions
OD-10 … OD-21 as they are recorded in the canonical register
(`docs/decisions/OWNER_DECISIONS.md`, branch `docs/owner-decisions`,
verified at `main` = `63da0e3`). No register row is amended, no calibrated
artifact is modified, and no data file is edited by this dossier.

## Read sources (byte references)

| Source | Location | Content used |
| --- | --- | --- |
| Canonical register | `63da0e3:docs/decisions/OWNER_DECISIONS.md` (branch `docs/owner-decisions` is canonical per Round-13 directive 5) | quoted substance of OD-10 … OD-21 (Part 1 consolidated register, `AUDITOR-DECIDED 2026-10-20`), plus §2.2 (OD-14 V4 correction) and §2.4 (OD-17…21 substance) |
| v24 coordinate manifest | `97d62a8:backend/calibration/runs/v24/coordinate_manifest.json` — SHA-256 `9722c3e45d316ff58ed326df074eec138b81ceb520191a492088cc0ddef63762` | 50 target pins (40 SHA-256-verified re-fetches, 10 historical `model_v4` failures re-probed) |
| v24 run manifest | `97d62a8:backend/calibration/runs/v24/run_manifest.json` — SHA-256 `f602351f4844531d409cfc256939508103e3c2d4f57fe21d7f77c34494aa4036` | run provenance pins referenced by the regeneration script |
| Delivered C5 record | `97d62a8:evidence/engine/C5-verdict-change.md`, `c5/*` (protocol shape and old-engine record parity target) | comparison protocol; the regeneration's cross-lineage section reproduces its old-engine rows (1AKE atoms 3312/det 9/10; identical LAW-182 observed values; LAW-182 advisory relabels) |
| Delivered clarification note | `97d62a8:evidence/engine/owner-clarification-od-03-04-20.md` | interpretation/scope text (NOT treated as a register amendment — see conflict record) |
| Superseded v24 approval note | register §2.2 (V4 correction; `backend/calibration/runs/v24/calibration_v24_SUPERSEDED.md` marked historical) | basis for keeping all verdict tables diagnostic-only |

## Per-decision reconciliation

Quoted substance lines are the register's own text (Part 1). "Candidate
state" records what round-26 work did or did not do.

### OD-10 — no thresholds from the contaminated cohort

> "Do not authorize fixed resolution, R-free, completeness, sample-size,
> or split-fraction thresholds from the contaminated historical cohort."

Candidate state: honored. The candidate touches no calibration input,
no threshold constant, and no split. The only governance edits are the
reviewed registry hunks (LAW-105/130/155/200 wording + scope) and their
re-pin; `LAW_CANON` thresholds/operators are byte-identical to `main`
(verified by `test_law_registry` pin `069a6c99…`, canon `6a9cd4b4349b81de…`
and matrix `ed1731ac…` literals untouched in `toscanini-gate.yml`).

### OD-11 — AF-0000000365840314 (P0DTC2) migration record

> "Retain an explicit, source-backed migration/identity record for
> `AF-0000000365840314` (P0DTC2). Keep it out of the authoritative cohort
> until the mapping and model identity are validated. Never silently
> substitute IDs."

Candidate state: honored. The round-26 regeneration writes
`evidence/engine/c5/AF-P0DTC2-F1.json` — a portable source-backed record
naming the served entity id `AF-0000000365840314`, the exact
`https://alphafold.ebi.ac.uk/files/AF-0000000365840314-model_v1.cif`
fetch, and the SHA-256 — the same identity pointer the delivered record
preserves (`97d62a8:evidence/engine/c5/AF-P0DTC2-F1.json`). The target
never silently reuses another accession's bytes, and the supplementary
cohort it belongs to is labeled supplementary (AlphaFold-DB-current),
never authoritative.

### OD-12 — P62988 inactive/no-valid-AlphaFold-entry

> "Record the inactive/no-valid-AlphaFold-entry case (P62988) as an
> explicit acquisition/identity failure. Exclude it from the authoritative
> cohort until a valid source-backed disposition exists."

Candidate state: honored. `evidence/engine/c5/AF-P62988-F1.json` records
the probe failure verbatim (`HTTP Error 404: Not Found`), measured on the
same `alphafold.ebi.ac.uk/api/prediction/P62988` endpoint; 9/10
supplementary acquisitions proceed, this one explicitly fails and is
excluded from every evaluation and from both the pinned and supplementary
cohort tables.

### OD-13 — pin exact model/version; no `latestVersion` movement for calibration

> "No moving `latestVersion` reference for calibration. Pin the exact
> model/version, entity identifier, validated HTTPS host and redirect
> policy, accepted format/size, coordinate SHA-256, and run-manifest
> provenance."

Candidate state: honored. No v25 calibration run exists or is prepared;
the supplementary records (`AF-*.json`, `af-v6-supp-*.json`) are
descriptive diagnostics of AFDB's current served state (e.g. entity
`AF-P68871-F1` version 6 as served on regeneration day) and are never
used to derive, authorize, or move any calibration pin. The v24 pins at
the top of this dossier are the only pins referenced for verification.

### OD-14 — no V4 / historical verdict-table validation (core of this dossier)

> "Do not approve V4 or the historical verdict-change table as validated
> scientific evidence. The historical 18° threshold is not validated, and
> the 1PPE VETO→PASS change is not approved for a verdict-changing merge.
> Keep those results diagnostic-only until the corrected cohort and
> independently specified evaluation protocol pass review."

Candidate state: honored. Everything round-26 produced around verdicts is
**diagnostic**: `evidence/engine/C5-verdict-change-f962c6abfe1d.md` opens
with "**Status:** diagnostic evidence for owner review. Not a threshold
validation." and repeats the constraint. The attributed comparison
(candidate's own base `63da0e3` vs head `f962c6a`) yields **0/40 v24-pinned
and 0/9 supplementary verdict changes**, so the reviewed port produces no
new verdict-changing outcome that could be mistaken for validation. The
historical 1PPE VETO→PASS item is untouched: no merge, no calibration
change, no threshold change. The cross-lineage section documents two
verdict changes (1PPE, 5T4O) that trace to the P0/WP-SCI-1
measurement-definition package pre-dating the round-25 window — recorded
as lineage documentation, equally not a validation claim. The SPR
(schema) and code paths that COULD reintroduce the 1PPE sort-repair are
explicitly listed as not-ported (port rationale §4).

### OD-15 — read-only review branches preserved

> "Retain `science/p0-correctness-2026-10-19` and
> `chore/hygiene-2026-10-19` as read-only review branches."

Candidate state: honored and out of scope that those two tips
(`ebc1fd79`, `abc450ea`) were never referenced for mutation. The
round-26 work reads two other frozen references only:
`engine/od-scope-v25` @ `97d62a8` and lineage rev `c78ef15` (the
delivered C5 old-engine baseline); both stay untouched on the remote,
and `main` is never force-pushed (the candidate branch push correction
is disclosed in the rationale §2).

### OD-16 — hard stop sequencing

> "Keep the hard stop. Sequence governance/security and code review
> first; then source-backed cohort repair and prespecified independent
> evaluation; then exact-candidate CI and review; … No merge, deploy, v25
> run, held-out evaluation, threshold/cohort change, or golden-artifact
> regeneration is authorized now."

Candidate state: honored. The round-26 work is branch-local review
material only: new review branch pushes, CI on the candidate, and audit
documents. No merge to `main`, no deploy, no v25 run, no held-out or
threshold/cohort/golden work. Evidence re-generation reads frozen pins
and writes additive `evidence/engine/` diagnostics.

### OD-17 — quarantine historical v23 metadata and results

> "Quarantine the complete historical v23 cohort metadata from
> calibration and threshold derivation. Preserve
> `backend/calibration/dataset_v23.py` and all historical results
> byte-for-byte."

Candidate state: honored. `backend/calibration/dataset_v23.py` and the
v23 result tree are untouched by the candidate diff (commits `818717a`,
`ff22b72`, `f962c6a` touch no calibration file); the v24 manifest under
`backend/calibration/runs/v24/` is only READ (pinned SHA verification).

### OD-18 — 1I0T / 1DDL exclusions

> "Exclude `1I0T` (DNA-only) and `1DDL` (protein plus RNA) from the
> protein-only core cohort. Preserve both in a separately classified
> audit set."

Candidate state: honored. The diagnostic v24-pinned evaluation table
still LISTS 1I0T/1DDL rows (as the delivered record does) but nothing in
the candidate derives any cohort, stratum, or threshold from them; no
corrected cohort is authored (that is OD-17/OD-19 future work, formally
blocked). The 1I0T row's law-status changes (LAw-110 NOT_MEASURED /
LAw-182 NOT_APPLICABLE for its empty predicted-pLDDT core) are recorded
as diagnostic outcomes of the reviewed scope change, not as cohort edits.

### OD-19 — quota/strata recalcs stay blocked

> "Discard the old quota assumptions for future calibration. Recalculate
> strata only after metadata and eligibility reconciliation. Prespecify…"

Candidate state: honored by absence — no strata, quota, sample-size,
uncertainty, or split work exists anywhere in the round-26 artifacts.

### OD-20 — NMR separation + no heuristic veto (conflict record)

> "Keep NMR separate from resolution-based X-ray/cryo-EM strata. …"

Candidate state: honored. The canonical register remains the only source
of calibration policy; no calibration change of any kind is made. On the
code side, the reviewed advisory outcome (heuristic LAw-155/182/200
statuses `FAIL (Advisory)`, adjudicator method `heuristic`, non-vetoing
row contribution) is implemented exactly as the reviewed commit
certified in its own tests (`backend/tests/test_law_182.py` at
`360373d`), without fabricating a register row — full reasoning in
`conflict-od-20-clarification.md`.

### OD-21 — no word-count gate

> "There is no verified project requirement establishing either 5,460 or
> 5,470 words."

Candidate state: honored. Round-26 documents are sized by content; no
length target is applied anywhere.

## Summary

Every register row OD-10 … OD-21 is either untouched-by-design
(OD-10/13/16/17/18/19/20/21) or satisfied by explicit, source-backed
records the round-26 work adds (OD-11/12/14/15). The candidate branch
contains: the reviewed round-25 port, its C5 diagnostic regeneration with
the base-vs-head attribution evidence, the port rationale, the OD-20
conflict disclosure, and this dossier — no merge, no deploy, no
calibration, no register amendment.
