# Security record — Gemini API key revocation and removal (H-1 / A-02)

**Date:** 2026-10-20
**Decision authority:** H-1 (Gemini removed completely) and A-02 (Gemini key
revocation), auditor acting under the owner's delegation, 2026-10-20.
**Classification:** credential incident record. Contains NO secret values.

## 1. What was found

A hardcoded Google (Gemini) API key literal was present in the repository at:

- `system_diagnosis.py`, line 8 (branch `main`, commit `2195a048…`, and
  inherited by earlier evidence branches).

The literal is REDACTED here by policy (EC §3: secret values are never
printed in reports, artifacts, or logs). Its first six characters were
`AIzaSy`, the standard Google API key prefix; the full value exists in prior
audit correspondence and in git history. Known-since record: the literal was
identified in the Round-9/Round-11 audit cycle (2026-10-20 or earlier) and
its presence contradicted a delivery README's "scan clean" claim; that claim
was retracted in writing (see the evidence branch's written retraction).

## 2. Revocation status (A-02)

- The owner confirmed IN CHAT (2026-10-20) that the key was revoked in the
  Google Cloud console.
- The revocation has NOT been independently verified in the Google console by
  the auditor or the coder; verification requires owner-side console access.
  Until the issuer confirms revocation/rotation, the credential is treated as
  COMPROMISED per EC §3 ("Treat a credential found in Git history as
  compromised until its issuer confirms revocation/rotation").
- Per H-1, no git history rewrite was performed; removal from the working
  tree does not revoke anything — revocation is the owner action recorded
  above (A-02).

## 3. Removal from code (H-1) — this branch

Branch `hygiene/h1-remove-gemini` (based on `main` @ `2195a048…`) removes
every Gemini reference outside `docs/audit-evidence/` and `docs/security/`:

| Surface | Removal |
|---|---|
| `system_diagnosis.py`, `internal_audit.py` | research scripts deleted (Gemini-dependent diagnostics) |
| `backend/tos/enrichment/gemini_compiler.py` | deleted; replaced by `narrative_compiler.py` (deterministic, provider-less, same interface) |
| `backend/main.py` | import switched to the deterministic narrative compiler |
| `backend/requirements.txt` | `google-generativeai` dependency removed |
| `backend/requirements-lock.txt` | `google-generativeai==0.8.6` and `google-genai==1.57.0` removed (lock hand-edited; full lock regeneration is a documented follow-up) |
| `.github/workflows/toscanini-gate.yml` | `google-generativeai` removed from the CI pip install set |
| `docker-compose.yml`, `docker-compose.override.yml` | `GEMINI_API_KEY` environment entries removed |
| `helm/toscanini/templates/secret.yaml`, `helm/toscanini/values.yaml` | `GEMINI_API_KEY` secret/value entries removed |
| vendored `control-atlas-saas/` copy | same treatment (compiler replaced, main.py import/usage, requirements, compose, SAAS_HANDOFF architecture section) |
| `backend/tos/glossary/epistemic_definitions.py` | PIL-NAR-11 description updated (deterministic narrative) |

**Acceptance (H-1):** `git grep -n "AIzaSy\|GEMINI\|generativeai"` returns
NOTHING outside `docs/audit-evidence/` and `docs/security/` on this branch —
verified at commit time. Historical/audit references remain unchanged inside
`docs/audit-evidence/` (evidence preservation, EC §11: never erase audit
evidence).

## 4. Residual risk and follow-ups

- The literal remains in git HISTORY (no rewrite per H-1). This is accepted
  under the no-rewrite decision; the compensating control is revocation
  (A-02) plus this record.
- Follow-up (owner-verified): confirm the key shows "revoked" in the Google
  Cloud console and record the confirmation here.
- Follow-up (hygiene, non-blocking): regenerate `requirements-lock.txt`
  through the normal toolchain so the removal is tool-verified rather than
  hand-edited; `google-auth*` transitive entries were deliberately left for
  that regeneration.
