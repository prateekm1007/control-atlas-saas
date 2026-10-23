# Constitution Reading Note — 2026-10-24

**Candidate:** `integration/od-scope-v25-2026-10-23` @ `8a9a9c9b156f8d05a76cf88550b9a36cfb609786`
**Branch:** `feature/evaluator-trust-boundary` (created from the candidate)
**Read window (UTC):** 2026-10-23T18:5–Z (fetch/identity) → 2026-10-23T19:05:02Z (first edit)

## Documents read in full (blob hashes at the candidate)

| Document | Blob SHA-1 |
|---|---|
| `AGENTS.md` | `85767237399ffb30f3e45894d9dfbe6f4240584f` |
| `constitution.md` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` |
| `ENGINEERING_CONSTITUTION.md` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` |
| `docs/coder/BRIEF.md` | `052a84ca829d321dc0841a4b1c8702145dcbf6b5` |
| `docs/decisions/OWNER_DECISIONS.md` | `12906bf5997a3e2c6b5c41e5d5c60fee14613274` |
| `docs/P0_STATUS.md` | `a31517abbc046c495930ba19ef7aeaf8bb72f2d0` |

## Enforced clauses relevant to this task

- **OD-16 hard stop:** no merges to `main`, no deploys, no v25 calibration run, no held-out evaluation, no threshold/cohort changes, no golden-artifact regeneration. This task does not perform any of these.
- **EC §2 (pre-coding protocol):** repository/branch/base SHA confirmed; constitutions read in full; preflight recorded (`docs/audit-evidence/2026-10-24/preflight.json`).
- **EC §4 (integrity):** provenance must bind authoritative results to input/content digest, engine/version, law-registry version/hash, configuration, and run status. `not measured` must be distinguished; synthetic data must not receive authoritative verdicts; a heuristic must not be described as a deterministic law.
- **EC §6/§12 (testing/evidence):** exact command, commit SHA, counts recorded; "not run" never reported as passed.
- **EC §7:** feature branches only; no direct push to `main`; no force-push; no merge.
- **EC §13 (prohibited shortcuts):** no gate weakening; no converting errors/absent metrics into normal verdicts; no silent law-semantics changes.
- **BRIEF §3:** verdict impact must be stated per patch (or "none" with argument).

## Scope of this branch

Evaluator trust-boundary hardening and evidence provenance in `backend/calibration/evaluate.py` and `backend/calibration/runner.py`: recursive `MANIFEST.sha256` verification (rejecting traversal/symlinks), source-byte verification of evidence references, cross-manifest target consistency, and exact `run_inputs` schema validation. Verdict/semantics: none (no law thresholds, IDs, or meanings change).
