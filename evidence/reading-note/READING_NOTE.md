# Constitution Reading Note (P0 program, items 3–18)

**Coder:** P0 correctness program, branch `science/p0-correctness-2026-10-19`
**Base SHA:** `296f10c` (science head at program start)
**Date:** 2026-10-19

## Sources read in full before any edit

| Document | Location | Version / binding |
|---|---|---|
| SaaS Engineering Constitution | `ENGINEERING_CONSTITUTION.md` on `policy/engineering-constitution` | 1.0.0 (Proposed canonical engineering policy) |
| Toscanini MDI Constitution (domain) | `constitution.md` on science branch | v6.1-Hardened (11 lines, 2 invariants) |

`AGENTS.md` is absent on every branch; it is not relied on. Per the
constitution's own change policy (§ Preamble), amendments are never smuggled
through implementation PRs — none of my commits amend either document.

## Engineering-constitution invariants I will enforce on this program

The Engineering Constitution states its rules per section rather than per
numbered ID. I assign stable reading identifiers `EC-<section>` below and
enforce every MUST/MUST NOT in the document; the list is exhaustive by
section, not paraphrased.

- **EC-2 (Pre-coding protocol):** record repo/branch/tree state and exact
  base SHA before editing; read both constitutions in full; read callers and
  consumers of every changed file; record expected behavior, risks, tests,
  and acceptance criteria before coding.
- **EC-3.1 (Fail closed):** no open-auth defaults, no dev secrets; missing
  security config must fail loudly, never silently permit.
- **EC-3.2 (Deny by default):** every route has an explicit access policy;
  callbacks are signed/audience-bound, never trust caller-supplied identity.
- **EC-3.3 (Authorize objects):** every read/write/artifact scoped to the
  authorized owner or service role.
- **EC-3.4 (Secrets):** no credential in source, fixtures, images, reports,
  logs, or remote URLs; runtime env only; history-exposed credentials are
  treated as compromised until rotated by the issuer.
- **EC-3.5 (Input validation):** type/size/format/content validated at
  boundaries; timeouts and resource bounds on costly work.
- **EC-3.7 (Sensitive data):** user structures and reports are proprietary;
  minimize exposure; keep payloads out of logs.
- **EC-4.1 (Provenance):** every authoritative result binds input digest,
  engine/version, law-registry hash, measurements, config, run status.
- **EC-4.2 (Honest non-success):** acquisition failure, empty input, and
  failed measurement are explicit non-success states — never a PASS, VETO,
  zero, or plausible score.
- **EC-4.3 (Test-only synthetic):** synthetic data never receives an
  authoritative verdict or enters production evidence stores.
- **EC-4.4 (Measurement honesty):** `not measured` is distinct from
  `measured and passed/failed` and `measurement failed`.
- **EC-4.5 (Law classification):** deterministic vs heuristic vs advisory is
  declared and visible; narrative output never overrides the engine verdict.
- **EC-4.7 (Registry integrity):** law IDs, thresholds, and meanings are
  never silently renumbered or redefined; conflicts stay explicit.
- **EC-5.2 (No duplicate sources of truth):** before deleting/moving code,
  inspect import graphs and runtime entry points; identical content does
  not prove interchangeability (cited re: dead `_law_155` and the nested
  `control-atlas-saas/` copy).
- **EC-5.3 (Contracts):** API schemas and persistent records are versioned
  or compatibly migrated.
- **EC-5.6 (Dependencies):** pinned from the correct manifest; backend and
  dashboard deps never confused (cited re: adding the mmCIF parser dep).
- **EC-6.1 (Regression first):** reproduce bugs with failing tests before
  fixing wherever feasible.
- **EC-6.4 (Test evidence):** claimed results include exact command, commit
  SHA, and pass/fail/skip counts; counts are never invented; "not run" is
  never "passed."
- **EC-6.3 (No test weakening):** no test is deleted, skipped, or weakened
  to obtain green CI.
- **EC-7.1/7.3 (Branch + gate discipline):** feature branches and PRs only;
  required checks run on the exact merge candidate; no deployment bypasses
  a failed gate.
- **EC-13 (Prohibited shortcuts):** no committing secrets; no
  errors-into-verdicts conversion; no law-meaning changes without explicit
  decision; no merge of unverified commits; no completion claims while
  material checks are unrun.
- **EC-14 (Task artifacts):** per-item `evidence/<item>/` directories with
  commands, outputs, SHAs, digests, and verdict-change tables.

## Domain-constitution invariants (v6.1-Hardened)

- **OBS-153 (Avidity Quantification):** ACTIVE. Bipodal success is contact
  density ρ(Bipod) ≥ 2.5 × ρ(Monopod), unique heavy-atom contacts within
  4.5 Å. Not touched by this program; no change made.
- **CRIT-154 (Structural Admissibility):** ACTIVE. Two-Stage Warhead-RMSD
  < 1.0 Å against clinical benchmarks (PDB: 5NIU); non-physical deformation
  (steric clash < 1.2 Å) is discarded. Not touched by this program; no
  change made.

## Conflicts recorded (not silently resolved)

1. `constitution.md` on science is a 2-invariant MDI domain policy; the
   Engineering Constitution on `policy/engineering-constitution` is not
   merged into any branch that deploys. Governance-file location is an
   owner gate (checklist item under Governance).
2. The Engineering Constitution §4 requires absent measurements to be
   explicit non-success states; the live adjudicator still carries
   malformed-row defaults of `FAIL` (present-but-statusless rows). Item 8
   aligns these to NOT_MEASURED.
3. `AGENTS.md` is referenced in audit direction but absent from every
   branch — recorded here so the directions' file references stay honest.

## Attestation

No automated invariant was weakened to make this program land. Every test
count reported in `docs/P0_STATUS.md` comes from an executed command whose
exact invocation is recorded in the item's evidence directory.
