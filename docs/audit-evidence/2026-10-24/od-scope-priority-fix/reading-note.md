# Constitution-first reading note — OD-scope priority-fix phase

**Date:** 2026-10-23T20:48:52Z (read window ~2026-10-23T19:52Z–20:49Z)
**Branch:** `fix/od-scope-residue-coverage-advisory` (derived from `integration/od-scope-v25-2026-10-23` @ `8a9a9c9b156f8d05a76cf88550b9a36cfb609786`)
**Base:** `main` @ `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`

## Governing documents read in full

| Document | Blob SHA-256 |
|---|---|
| `AGENTS.md` | `5e320939a11135a0d8aed5caa231df6aeaf633d724fcb4dba97b5909e51420b4` |
| `constitution.md` | `cfdd5ba3babda31513b55894a49597fbfb18c2e870535a068c8b0ea423c53ce8` |
| `ENGINEERING_CONSTITUTION.md` | `e0f256f78f2c1256228c4286311cba15800dcd5a6c8005583508ee66d27983f5` |
| `docs/coder/BRIEF.md` | `8f07edfbba9d3e57de623e242125e942369bc5d60df5805343b26b9048c52e07` |
| `docs/decisions/OWNER_DECISIONS.md` | `59ccd5e84176a8a610c8c234cbb95fec4734478cc1f35668599ad2e3aef3172a` |
| `docs/P0_STATUS.md` | `3f129dd18b9e5cf8888b763158cda42de270a8429a48ed7cc9a3aa675e810e0f` |
| `backend/calibration/protocol/evaluation_protocol_v1.md` | `c880086f04f59b577a34ead13eafa2de6cc5ae337059bc16b1ed87434a724806` |

## Rules enforced in this phase

- EC §2 pre-coding protocol: read constitutions in full, record hashes, pin base,
  map callers, record a machine-readable preflight.
- EC §4: distinguish measured vs inferred/missing/failed; never convert an absent
  measurement into a favorable zero; heuristics are not deterministic laws.
- EC §6/§12: reproduce with failing regression first; record exact command/SHA/counts;
  never invent counts.
- EC §13: no weakening/deleting gates; no law threshold/meaning changes without a
  decision; no force-push; no merge/deploy.
- OD-16: hard stop in force (no merge, deploy, v25, held-out, threshold/cohort,
  golden-artifact regeneration).
- D1–D5 (register rev 12): one integration line from `main` `63da0e3` via content-port;
  D3 LAW-155 advisory (not retired); D4 Governance Acceptance diagnostic-only; D5
  calibration targets preregistered.

## First edit timestamp
2026-10-23T20:49Z (first source edit: `backend/tos/ingestion/processor.py`).
