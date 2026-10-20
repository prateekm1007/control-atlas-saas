# Incident record — push/merge governance deviations (append-only, EC-18)

- **Grammar:** `## <ISO-8601 UTC> | <ENTRY-ID> | <verb> | <subject> [| RECONSTRUCTED]` — the
  same evidence-log grammar as `docs/evidence/push-log/PUSH_LOG.md`. Entry IDs use the `INC-`
  namespace, which is disjoint from push-log `P-<n>` IDs and is never reused.
- **Append-only:** entries are never edited or deleted; corrections are new entries that name
  the entry they correct. The parent-prefix check (`scripts/ci/evidence_sequence_check.py
  --parent-prefix`) applies to this file from its second tracked version onward.

## 2026-10-20T13:32:54Z | INC-1 | record-incident | PR #17 merge to main bypassed the required review path (administrator push; merge endpoint was blocked)

- **What happened (verified facts):** PR #17 (`proposal/ec16-ec26-enforcement-rebased` →
  `main`) was merged as `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd`. Branch protection on
  `main` required one approving review, and the merge endpoint was blocked for the session's
  credential; the merge was completed by an administrator push that bypassed the required
  review path. The ordinary gate passed on the exact merge candidate (run `39415713844`,
  both jobs success), but the merged tree retained mainline defects (provider-SDK references
  in two `proposals/remediation/` files, stale remediation documentation, incomplete manifest
  coverage) — recorded by the Round-18 audit.
- **Classification:** procedural deviation (review-path bypass). NOT a rewrite: no history
  was rewritten or deleted; the merge is a normal merge commit and stays in force. The
  defect is that the required human review did not happen before the merge.
- **Why it is an incident:** the review requirement is the control that catches exactly the
  class of defects that survived in the merged tree; an administrator push removes that
  control regardless of gate color. The deviation is recorded so the bypass cannot be
  normalized into practice.
- **Remediation (binding on future work):**
  1. No administrator push may bypass branch protection again. Future merges MUST use the
     required review-and-check path (approving review + required checks on the exact merge
     candidate).
  2. No existing historical evidence is rewritten, deleted, or "corrected in place" to
     reflect this incident; this append-only record is the correction mechanism.
  3. The Round-18 repair proceeds as a fix-forward PR from the exact `main` tip
     `f535a71b85d90ab8ff92ac8613a05b3fb77f5dfd` (no reset, no force-push, no re-merge of
     divergent history).
  4. No current activation PR is approved for merging by the Round-18 audit; activation PRs
     merge only after their exact final candidates are green AND the required review path is
     available and used.
- **Credential note:** the credential used is recorded by NAME only (`OWNER-PAT-1`); token
  values are never written anywhere (standing credential-containment rule). The PAT
  disposition remains owner-controlled under A-01 (`ACCEPTED RISK`; no rotation reminders).

## 2026-10-20T13:32:54Z | INC-2 | record-collision | two competing P-015-CONFIRM records exist in two preserved push-log histories (documented, not merged)

- **Facts:** the published `main` history (`f535a71`) carries `P-015-CONFIRM` for
  `proposal/ec16-ec26-enforcement-rebased` with `remote_sha=0f1e9a8856bafec592c657270e5799734ed61d6a`
  (push accepted 11:27:51Z, 2026-10-20). The separate preserved history
  `proposal/ec16-ec26-enforcement` @ `c4a0fcc3357930b7eba61e056489c804b16d2277` carries its own
  `P-015-CONFIRM` for `proposal/ec16-ec26-enforcement` with
  `remote_sha=a0e321c0ef08cd3430e7910936075994791fd8df` (re-verified at that CONFIRM's write time).
  The same entry ID therefore denotes two different confirmations in two different lineages.
- **Disposition:** NEITHER log is copied onto the other history; the two `P-015-CONFIRM`
  records are NOT combined; both histories are preserved untouched. The collision is
  documented here and in `docs/evidence/push-log/PUSH_LOG_RECONCILIATION.md` +
  `PUSH_LOG_RECONCILIATION.json`, which also classify the confirm SHAs that no longer sit on
  their branch's current remote lineage. The push-log protocol checker
  (`scripts/ci/evidence_sequence_check.py --pushlog-schema`) validates pairing per branch
  against the log on the CURRENT history only and consumes the reconciliation record for
  superseded SHAs; it never merges records across histories.
- **ID discipline:** all NEW push-log events from Round-18 onward use globally unused
  `P-<n>` IDs continuing after the highest ID used in ANY preserved history (the highest
  existing ID anywhere is P-016 on the `c4a0fcc` history; Round-18 events therefore start at
  P-017). Incident IDs use the `INC-` namespace.
