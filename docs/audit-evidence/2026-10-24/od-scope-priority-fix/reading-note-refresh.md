# Constitution-first reading note — refresh for the P2-API-fix cycle (2026-10-24)

**Read window:** 2026-10-24T01:05Z–01:15Z (first edit of this cycle: 01:15Z)
**Branch:** `fix/od-scope-residue-coverage-advisory` (candidate = the introducing commit of this cycle, resolved on the next push)
**Base:** `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`

## Working-tree blob hashes (measured 2026-10-24T01:10Z; unchanged from prior refresh)

| Document | SHA-256 |
|---|---|
| `AGENTS.md` | `5e320939...` |
| `constitution.md` | `cfdd5ba3...` |
| `ENGINEERING_CONSTITUTION.md` | `e0f256f7...` |
| `docs/coder/BRIEF.md` | `8f07edfb...` |
| `docs/decisions/OWNER_DECISIONS.md` | `b0c64b08...` |
| `docs/P0_STATUS.md` | `3f129dd1...` |

## Recorded conflict (preserved, not resolved here)
Working tree `constitution.md`=`cfdd5ba3...`, `ENGINEERING_CONSTITUTION.md`=`e0f256f7...`; register/auditor-cited canonical blobs are `211ad2c5...` / `b6bb35e8...`. Both preserved; canonical-document authority to be resolved via the reviewed governance process.

## This cycle
- P2 API-path fix: `main.py` serializes unavailable coverage as JSON `null` (no `round(None)`); `tier1_measurements.py` early-return paths return `coverage=None`; response schema `verdict.coverage_pct` is `Optional[float]`; `pdf_generator.py` renders `None` coverage without arithmetic/crash; added API regression `TestSourceUnavailableApiContract` (no server error, null coverage_pct, LAW-105 MEASUREMENT_FAILED, verdict INDETERMINATE).
- Calibration protocol preregistration rule recorded (one-sided 95% exact Clopper-Pearson upper bound ≤5% false vetoes, overall + per stratum; 30/stratum = feasibility floor).

## Rules enforced
EC §2; EC §4 (record conflicts); EC §6/§12 truthful counts; EC §13 no law/gate changes; OD-16 hard stop.
