# H-1 on the P0 stack — evidence (Round 14, direction 3)

Base: `science/p0-correctness-2026-10-19` @ `ebc1fd79f20b9a6b6a262a70ce6942ba069d19e9`
(the P0 stack, per the auditor's binding H-1 decision). The Round-13
H-1 commit (`a761d72`, based on `main`) is superseded by this branch;
it remains preserved in the Round-13 durability bundle and patch series
(`download/round13-branches/`, SHA-256 in `SHA256SUMS`).

## Ruling and how it was applied

Auditor decision (Round 14, binding): "Gemini removal happens on the P0
stack, based on `science/p0-correctness-2026-10-19`. The narrative
compiler stays only if `generate_v21_dossier` depends on it, with
evidence. Otherwise it is deleted."

**The dossier does not depend on the narrative compiler — the compiler
is deleted.** Evidence (call-path analysis on the base `ebc1fd79`):

1. `backend/tos/forensic_artifacts/pdf_generator.py::generate_v21_dossier`
   (line 132) reads, via `payload.get(...)`: `verdict`, `tier1`, `tier3`,
   `characterization`, `provenance`, `confidence_meta`, `governance`,
   `strategic_math`. A full-file search for `witness`, `narrative`,
   `ai_model`, `executive`, `deep_dive`, `recommendation` over
   `pdf_generator.py` returns no payload-field access (only an unrelated
   badge comment and a fixed refinement note). `witness_reports` is
   never read by the dossier generator.
2. The only consumers of the compiler output were
   `backend/main.py:206` (`witness_reports`) and `:209`
   (`ai_model_used`) — both removed with the compiler. Both schema
   fields remain `Optional[...] = None` in `tos/schemas/response_v1.py`
   (backward compatible; absent keys validate as null).
3. The dashboard (`dashboard/dashboard.py:421-431`) renders
   `witness_reports` with explicit fallbacks (`"Analysis pending."`) —
   the honest degraded state (EC §10), not a fabricated narrative.
4. No test on the base imports the enrichment compiler (verified by
   grep over `backend/tests`); `tos/witness/contract.py`
   (WitnessContractEnforcer) is an unrelated module and is untouched.

## What changed (19 files touched at the base, full sweep)

- `backend/tos/enrichment/gemini_compiler.py` **DELETED** (ruling,
  Case: dossier-independent). No replacement module.
- `backend/main.py`: compiler import and the two payload lines removed
  (comment records the ruling and the evidence pointer).
- `backend/requirements.txt` / `requirements-lock.txt`:
  `google-generativeai>=0.3.0` and `google-genai==2.31.0` /
  `google-auth==2.62.0` pins removed (the provider SDK and its only
  auth dependency; no other lock package requires them — verified).
- `scripts/research/system_diagnosis.py`, `scripts/research/internal_audit.py`:
  provider SDK imports and credential probes removed (the pre-existing
  corrupted `models =` line in system_diagnosis.py went with its block;
  AFDB/dependency checks kept).
- `docker-compose.override.yml`, `helm/toscanini/values.yaml`,
  `helm/toscanini/templates/secret.yaml`: provider key wiring removed.
- `README.md`, `docs/DEPLOY.md`: provider env-var rows removed.
- `docs/P0_STATUS.md`: A-02 owner action rephrased (issuer-side
  revocation still pending); provider key removed from the required EC2
  secret set (no code path reads it).
- `backend/tos/engine/adjudicator.py` header comment and
  `backend/tos/glossary/epistemic_definitions.py` PIL-NAR-11 entry:
  provider mentions neutralized.
- Legacy vendored copy `control-atlas-saas/`: compiler DELETED;
  `backend/main.py` rewired to a deterministic inline summary; its
  requirements, `docker-compose.yml`, and `SAAS_HANDOFF.md` cleaned.
- Evidence, the Round-14 constitution reading note, and the
  incident record (`evidence/incident-push-2/`) added.

## Acceptance (auditor's exact criterion, run on this branch)

    git grep -n "AIzaSy\|GEMINI\|generativeai" -- ':!docs/audit-evidence' ':!docs/security'
    -> no matches (CLEAN)

Stricter sweep (case-insensitive
`gemini|generativeai|aizasy|genai|google.generative`, same excludes):
also CLEAN. PAT scan of the tree: CLEAN (one pre-existing
`docs/audit-evidence/2026-10-19/` file contains the pattern-text
`ghp_...` — dots, not a token; present on the base and passing the
science tip's gitleaks CI).

## Commands run and output (Python 3.12.15; gemmi 0.7.5; provider SDK absent)

1. `cd backend && python3 -m pytest tests -q --tb=line`
   -> **4 failed, 502 passed, 5 skipped in 10.34s**
   IDENTICAL to the pristine science base in this environment (4
   environmental TestUsageWidget failures — streamlit absent locally,
   green in CI which installs requirements-dev.txt). Zero regressions.
2. Hash locks (Stability Gate check): canon `6a9cd4b4349b81de`, matrix
   `ed1731ace09a9a3b…`, law-semantic vs
   `governance/law_semantic.sha256` — **all MATCH** (no law-semantics
   change; verdict impact: NONE).

## Verdict impact

NONE. No law, threshold, operator, scope, or adjudication change; the
governance hash locks are bit-identical. The removed layer was advisory
narrative only, which by construction could not and cannot alter the
authoritative deterministic verdict.
