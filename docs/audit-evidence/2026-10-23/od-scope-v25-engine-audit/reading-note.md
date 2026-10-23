# OD-scope v25 engine audit — required reading note

- Date: 2026-10-23
- Branch: `repair/od-scope-v25-engine-audit`
- Base commit: `ae8a53547a0433f9220ff656fc53a706591c3bd4` (`review/od-scope-v25`, verified remote tip)
- Candidate: under construction; the terminal commit SHA will be recorded in the exact-candidate verification report and in the P-092 confirmation after the final commit exists.
- Governing directive: owner directions of 2026-10-23 for the next legitimate candidate, including preservation of `ae8a535`, P-091 confirmation and P-092 intent, P0 engine/ingestion repair, calibration provenance hardening, and append-only historical evidence repair.

## Complete-read record

The following files were read in full in this session before any edit on this branch. Hashes were independently verified against the working tree at `2026-10-23T04:49:44Z`; the complete reads precede this note.

| File | Bytes | SHA-256 | Git blob |
|---|---:|---|---|
| `AGENTS.md` | 4021 | `5e320939a11135a0d8aed5caa231df6aeaf633d724fcb4dba97b5909e51420b4` | `85767237399ffb30f3e45894d9dfbe6f4240584f` |
| `constitution.md` | 622 | `cfdd5ba3babda31513b55894a49597fbfb18c2e870535a068c8b0ea423c53ce8` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` |
| `ENGINEERING_CONSTITUTION.md` | 20572 | `e0f256f78f2c1256228c4286311cba15800dcd5a6c8005583508ee66d27983f5` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` |
| `docs/coder/BRIEF.md` | 3671 | `8f07edfbba9d3e57de623e242125e942369bc5d60df5805343b26b9048c52e07` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` |
| `docs/decisions/OWNER_DECISIONS.md` | 72329 | `59ccd5e84176a8a610c8c234cbb95fec4734478cc1f35668599ad2e3aef3172a` | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` |
| `docs/P0_STATUS.md` | 35347 | `3f129dd18b9e5cf8888b763158cda42de270a8429a48ed7cc9a3aa675e810e0f` | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` |

The governing directive requires the constitution to be re-read at a new phase or when direction conflicts. No conflict has been identified; the controlling constraints below are applied directly.

## Exact-candidate baseline

- `ae8a535` is published; it must not be pushed again. This branch starts from it.
- Stability Gate run `39897260982` succeeded: Secret Scan job `120230456596`; Unit Tests + Hash Lock job `120230456813`.
- Governance Acceptance run `39897260952` failed only the push-log remote-verification job `120230457210`. EC-22 `120230456874`, AT-01..AT-34 `120230457021`, manifest/forbidden `120230457056`, EC-23 `120230457093`, EC-26 `120230457104`, and EC-25 `120230457140` succeeded. EC-21 `120230457980` and EC-16/19/20 `120230458312` were skipped.
- The failure is 22 historical PUSH_LOG records across 12 deleted refs. Base `e71cec3` run `39825159704`, job `119986401764`, has the same 22 missing-ref failures. This is pre-existing, not a candidate regression.
- Existing P-074/P-084/P-085 `push-not-proven` records do not suppress the failure: `scripts/ci/evidence_sequence_check.py` checks live branch availability before honoring reclassification, and the schema fails closed on missing remote verification.
- P-091 remains `PENDING`; no P-091-CONFIRM has been appended. P-092 has not been introduced.

## Policy clauses that bind this work

- `OD-16` remains in force: no merge, deployment, v25 calibration execution, held-out evaluation, threshold/cohort change, or golden-artifact regeneration.
- Roadmap D2's OD-14 claim does not supersede the controlling OD-14 decision; C5 is diagnostic only.
- Governance Acceptance and evidence-integrity checks must not be removed or weakened absent a demonstrated replacement and explicit controlling decision.
- The append-only push log is evidence, not editable narrative. Prior entries must be preserved byte-for-byte; new state is appended only at the end. No global waiver, invented citation, or rewritten historical record.
- Historical-ref repair must distinguish dead refs, unavailable run citations, and independently provable commits. Only machine-verifiable live evidence may be reclassified.
- `P-091-CONFIRM` must cite the actual remote ref and push-event workflow evidence for `ae8a535`: workflow/run IDs, event type, branch, head SHA, conclusions, and required job results. Do not manufacture evidence.
- `P-092-INTENT` must be introduced only by this candidate's final commit and declare `expected_remote_tip_before_push=ae8a53547a0433f9220ff656fc53a706591c3bd4`.
- OD-07 audit identity must bind full SHA-256 coordinate identity, canonical-law hash, and modality-matrix hash through canonical serialization. Atom-order invariance must be explicit and tested.
- Canonical ingestion must use an explicit residue identity, model policy, altloc policy, deterministic occupancy tie-break, and declared polymer-HETATM treatment. PDB and mmCIF must produce equivalent measurements.
- Coverage must be source-backed. Malformed, discarded, or unparsed required records must not inflate coverage or produce a favorable experimental shortcut.
- Geometry must not treat sequence adjacency as peptide-bond proof. Chain breaks and missing insertion-code intermediates must break continuity. Missing/invalid measurements must not become favorable zeros.
- pLDDT must be observation-backed. B-factors must never be relabeled as pLDDT.
- Calibration labels, family assignments, exclusions, and run artifacts require independent source/byte verification. Hash-shaped strings are not proof. Unpinned or diagnostic-only model inputs must never become authoritative.
- No force-push, no rewrite of published history, and no unrelated commit after the P-092 terminal intent.
- The owner PAT is rotated by the owner; it must not be exposed, committed, logged, or persisted.
- A local pass is not remote CI. The exact pushed SHA and every required push-event job must be verified after the push.

## Affected callers and contracts

- `backend/tos/ingestion/processor.py` -> `backend/main.py::_run_physics_sync` entrypoint and `backend/calibration/runner.py::_audit_structure`.
- `backend/tos/engine/tier1_measurements.py::run_full_audit` -> `main.py`, `calibration/runner.py`, and engine regression tests; the `(results, coverage, is_experimental)` tuple shape must remain compatible.
- `backend/tos/engine/coord_math.py::get_atom_pos` -> Tier1 measurements and LAW tests; canonical residue/atom selection must not change call sites.
- `backend/main.py` audit identity -> API `governance.audit_id`, persisted audit/refinement records, and `backend/tests/test_ids.py`.
- `backend/calibration/evaluate.py` provenance checks -> calibration readiness and evaluator reports; no threshold validation is authorized.
- `backend/tos/engine/batch_processor.py` accepts ZIPs of PDB/ENT bytes; PDB-first behavior remains supported while canonical ingestion adds mmCIF.
- `backend/calibration/protocol/*` is not modified unless required for truthful provenance; `calibration_status=not_approved` and `supports_threshold_validation=false` remain.

## Planned change groups

1. Canonical ingestion/atom selection: PDB and mmCIF through pinned gemmi; HETATM, waters/ligands, altloc, occupancy, model policy, residue identity, and deterministic selection.
2. Source-backed coverage and required-observation accounting, including discarded/malformed source records.
3. OD-07 canonical coordinate hash and full audit ID, with atom-order invariance and component-hash tests.
4. LAW-120/125/130/135/160 continuity and bond-validity corrections, finite-input propagation, and regression tests.
5. P1 LAW-100/145/150/155/182/195 and confidence/B-factor audit; no favorable missing-evidence shortcuts.
6. Calibration provenance verification: resolved source bytes/digests, mandatory integrity manifest, pinned input policy, defect-injection tests.
7. Append-only historical evidence repair and truthful P-091-CONFIRM/P-092-INTENT.
8. Machine-readable audit report, exact-SHA full test/check results, and post-push workflow verification.

## Acceptance criteria

- The reading note and preflight exist before any source/test/workflow/push-log edit.
- PDB and mmCIF equivalents produce identical canonical structures/measurements for the same coordinates.
- HETATM/altloc/occupancy/model policy is declared, deterministic, and covered by tests; waters/ligands are explicitly excluded or classified.
- Coverage uses a source-backed denominator; missing/malformed/discarded records cannot increase it; experimental structures do not receive a fabricated 100%.
- Audit IDs are full SHA-256 values binding coordinate identity, canon hash, and matrix hash; component changes change the ID; atom order does not.
- LAW-120/125/130/135/160 do not measure across chain breaks, missing insertion codes, or non-bonded adjacency; invalid measurements fail closed.
- Confidence/pLDDT is observation-backed and never sourced from experimental B-factors.
- Calibration source references resolve and match actual bytes; missing integrity manifests fail closed; unpinned inputs cannot be authoritative.
- Historical push-log repair is append-only and machine-verifiable; prior P-087/P-088 bytes are preserved.
- P-091-CONFIRM cites exact post-push evidence; P-092-INTENT is the sole terminal intent and is introduced by the final commit.
- Targeted and full suites, acceptance tests, manifest, secret scan, workflow lint, hash locks, and append-only checks pass on the exact candidate SHA.
- The live remote is `ae8a535` immediately before a fast-forward push; post-push live SHA and required jobs are verified.

## Pre-change attestation

No source, test, workflow, push-log, manifest, or governance content has been modified on this branch before this note and the accompanying preflight were written. The only prior actions were reading the required sources, verifying their hashes, and creating this branch at the unmodified base.
