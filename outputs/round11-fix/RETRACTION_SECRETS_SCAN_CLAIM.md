# RETRACTION — Round-11 delivery secrets-scan claim (audit failure 2)

**Date:** 2026-10-20 (fix delivery, session-local UTC+05:30)
**Retracting party:** coder (this session)
**Delivered on:** `evidence/round11-fix-2026-10-20` (branch based on
`origin/main` @ `2195a048e87dcec4aa63e8dc7981b3bb6d237ca0`, carrying the
prior Round-11 delivery tip `7e3a3395be2d0719facf5c61ac532a0d347a49d3`
plus this correction)

## 1. What is retracted, verbatim in substance

The Round-11 delivery README (`outputs/DELIVERY_README.md` on
`evidence/round9-r2-2026-10-20` @ `7e3a339`, commit message "evidence:
Round-11 (U1-U9) outputs delivery") presented the delivery as secrets-
clean: its SHA-256 tables and "Contents" section describe a verified,
self-contained delivery, and the README's closing sections state that
nothing outside the `outputs/` folder was pushed. The implied and
intended claim — **that the delivered branch was checked for
credentials and none are present** — is FALSE, and it is hereby
RETRACTED in writing.

The Round-11 audit's finding is confirmed against the actual bytes:

- The branch `evidence/round9-r2-2026-10-20` is based on `origin/main`
  @ `2195a048…`, whose tree contains `system_diagnosis.py` at the
  repository ROOT.
- `system_diagnosis.py`, line 8, contains the hardcoded literal
  `key = "AIzaSyD1d8…"` (a Gemini API key; the full value is NOT
  reproduced in this retraction — EC §3/§13: no secrets in reports,
  artifacts, logs, prompts, or issue text). The literal is therefore
  PRESENT in every commit of that branch, including the delivery tip
  `7e3a339`, because the branch inherits main's tree.
- The delivery README claimed a clean state for the delivery without
  disclosing this. That claim was wrong.

## 2. Factual posture (verified this session, GET-only + local objects)

| Fact | Value |
|---|---|
| Key literal location on the delivered branch | `system_diagnosis.py` line 8 (repo root), inherited from `main` @ `2195a048…` |
| Is the key revoked? | The Round-11 audit states the key is already revoked; revocation is an issuer-side fact the coder cannot verify and relies on the auditor's statement for |
| Known-since | The Worklog record of the original P0 baseline (Task 1-2, 2026-10-19) already listed "system_diagnosis.py:8 hardcoded Gemini key" as the single secret-grep hit at HEAD |
| Existing fix, not merged | `hardening/WP-1-remove-hardcoded-credentials` @ `ab5bc321…` and the science-branch patch `0001-WP-1-Remove-hardcoded-Gemini-key-add-pinned-gitleaks…` remove the literal; they remain owner-gated and UNMERGED |
| Fail-closed pattern already exists | `scripts/research/system_diagnosis.py` on the science branch (`296f10c…`) reads `GEMINI_API_KEY` from the environment and fails closed; the repo-root copy on `main` was never given that fix |
| New pushes this session | `evidence/round11-fix-2026-10-20` only — a NEW branch; no existing ref moved, nothing merged, PR #16 untouched |

## 3. Why the claim was wrong (honest accounting)

The Round-11 delivery verified the **`outputs/` folder contents** (the
delivered evidence artifacts) and presented that verification as if it
characterized the **branch**. It did not. The correct statement would
have been: "the outputs/ artifacts contain no credentials; the branch
tree, however, inherits main's `system_diagnosis.py:8` hardcoded Gemini
literal (revocation per the auditor's statement), which remains until
the WP-1 fix is merged." That distinction — artifact-level clean vs
branch-tree state — is exactly the conflation EC §4 forbids (measured
vs inferred never conflated), and it is recorded here as a process
error, not smoothed over.

## 4. Corrective posture going forward

1. This retraction is delivered on the fix branch, linked from the fix
   README, and cross-referenced in the register (OD-15 incident-ledger
   context).
2. No re-push, history rewrite, or branch deletion of the affected
   delivery is performed (EC §11: never erase audit evidence or rewrite
   history). The literal remains in `main`'s tree/history until the
   owner approves the WP-1 merge (owner action). Git-history exposure
   is treated as compromise-until-revocation per EC §3; per the
   auditor's statement the key is already revoked, which closes that
   exposure.
3. Every future delivery README must carry an explicit BRANCH-TREE
   secrets statement (scan scope, exact scan command, and the known
   inherited literals) separate from the artifact-level statement.
