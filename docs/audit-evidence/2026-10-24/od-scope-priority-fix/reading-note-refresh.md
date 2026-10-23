# Constitution-first reading note — refresh for candidate f4485d8 (audit 2026-10-24)

**Date/read window:** 2026-10-24T00:45Z–00:55Z (first edit of this cycle: 00:56Z)
**Candidate:** branch `fix/od-scope-residue-coverage-advisory` @ `f4485d8badca1756bf010f90b01af479dcb52fe5`
**Base:** `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`

## Actual blob hashes in this working tree (measured 2026-10-24T00:50Z)

| Document | SHA-256 (tree) |
|---|---|
| `AGENTS.md` | `5e320939...` |
| `constitution.md` | `cfdd5ba3...` |
| `ENGINEERING_CONSTITUTION.md` | `e0f256f7...` |
| `docs/coder/BRIEF.md` | `8f07edfb...` |
| `docs/decisions/OWNER_DECISIONS.md` | `b0c64b08...` (includes Part-16 D6) |
| `docs/P0_STATUS.md` | `3f129dd1...` |

## Recorded conflict (no silent resolution)

The auditor's 2026-10-24 direction cites the canonical constitution blobs `constitution.md` = `211ad2c5...` and `ENGINEERING_CONSTITUTION.md` = `b6bb35e8...` (these match the register/P0_STATUS Round-14 references). The working tree actually contains `cfdd5ba3...` / `e0f256f7...`. This divergence is recorded here; both sets are preserved. Per EC §4 the effective runtime path (this tree's bytes) is what was read for this work; the register-cited blobs remain the canonical reference pending reconciliation on the canonical governance branch.

## Documents read in full (this session)
AGENTS.md, constitution.md, ENGINEERING_CONSTITUTION.md, BRIEF.md, OWNER_DECISIONS.md (incl. Parts 14/16), P0_STATUS.md, README.md, SAAS_HANDOFF.md, evaluation_protocol_v1.md, plus the affected implementation and tests for P1/P2 (processor.py, tier1_measurements.py, coord_math.py, main.py, adjudicator.py, identity.py, coverage/parity/advisory tests).

## Rules enforced
EC §2 pre-coding protocol; EC §4 (record conflicts, distinguish measured vs inferred); EC §6/§12 truthful counts; EC §13 no gate weakening / no law changes; OD-16 hard stop.
