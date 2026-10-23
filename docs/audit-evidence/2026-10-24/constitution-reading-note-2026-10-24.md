# Constitution Reading Note — 2026-10-24

**Candidate:** `integration/od-scope-v25-2026-10-23` @ `8a9a9c9b156f8d05a76cf88550b9a36cfb609786`
**Branch:** `feature/evaluator-trust-boundary` (created from the candidate)
**Reading anchors (UTC):** full reads completed before the first edit; attested by the
workspace-preflight record at `2026-10-23T19:05:02Z` and the first push intent
`P-093-INTENT` at `2026-10-23T19:09:59Z`. No read is claimed after a later edit.

## Documents read in full (blob hashes at the candidate)

| Document | Blob SHA-1 |
|---|---|
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` |
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` |
| `docs/decisions/OWNER_DECISIONS.md` | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` |
| `docs/P0_STATUS.md` | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` |
| `backend/calibration/protocol/evaluation_protocol_v1.md` | read in full (protocol scope) |
| `backend/calibration/protocol/evaluation_protocol_v1.1_addendum.md` | read in full (protocol addendum) |
| `scripts/ci/manifest_check.py`, `scripts/ci/evidence_sequence_check.py` | read in full (CI evidence-sequence scope) |

## Enforced clauses relevant to this task

- **OD-16 hard stop:** no merges to `main`, no deploys, no v25 calibration run, no held-out evaluation, no threshold/cohort changes, no golden-artifact regeneration. This task performs none of these.
- **EC §2 (pre-coding protocol):** repository/branch/base SHA confirmed; constitutions read in full; preflight recorded (`docs/audit-evidence/2026-10-24/preflight.json`).
- **EC §4 (integrity):** provenance must bind authoritative results to input/content digest, engine/version, law-registry version/hash, configuration, and run status. `not measured` must be distinguished; synthetic data must not receive authoritative verdicts; a heuristic must not be described as a deterministic law.
- **EC §6/§12 (testing/evidence):** exact command, commit SHA, counts recorded; "not run" never reported as passed.
- **EC §7:** feature branches only; no direct push to `main`; no force-push; no merge.
- **EC §13 (prohibited shortcuts):** no gate weakening; no converting errors/absent metrics into normal verdicts; no silent law-semantics changes; no fabricated evidence.
- **BRIEF §3:** verdict impact must be stated per patch (or "none" with argument).

## Scope of this change

- **Modified:** `backend/calibration/evaluate.py` (trust-boundary hardening) and `backend/tests/test_calibration_evaluate.py` (adversarial tests).
- **Read for context, not modified:** `backend/calibration/runner.py` (diagnostic-only; output run_type not changed).
- **Verdict/semantics impact:** none — no law thresholds, IDs, operators, semantics, cohorts, or hash locks change.

## Intended behavior, risks, file scope, test plan, acceptance, rollback

See the regenerated workspace preflight (`docs/audit-evidence/2026-10-24/preflight.json`) for
the full machine-readable plan: intended behavior, security risks, file scope, regression-test
plan, acceptance criteria, and rollback strategy. This note attests the constitution-first read
and the enforced clauses above.
