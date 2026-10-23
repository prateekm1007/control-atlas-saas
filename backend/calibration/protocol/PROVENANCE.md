# Calibration Protocol Provenance

## Status

- `calibration_status`: `not_approved`
- `supports_threshold_validation`: `false`
- This protocol is recovered as a **document of record only**. It is not a
  live gate in this candidate and does not authorize any held-out evaluation,
  threshold validation, cohort split, or golden-artifact regeneration. OD-16
  remains in force.

## Exact source

- Ref: `engine/od-scope-v25` @ `97d62a89f0efb2917e80a82808a930d34af7f492`
- Introducing commit: `863e6dd56156375c8bf48c9e6a802a0daf4cec14`
- Blob: `ca51ceb794320cf33dc7a83760d299cd18d1e69b`
- Size: 8,364 bytes
- SHA-256: `c880086f04f59b577a34ead13eafa2de6cc5ae337059bc16b1ed87434a724806`
- The blob is unchanged on `engine/od-scope-v25` through `97d62a8`.
- The WIP preservation commit `631f2bf` carries the same byte-identical blob.

The recovery copies the exact bytes above into
`backend/calibration/protocol/evaluation_protocol_v1.md`. Byte equality is
pinned by `backend/tests/test_calibration_protocol_recovery.py`.

## Scope and known gaps

- The review lineage retains its independently hardened
  `backend/calibration/evaluate.py` and `backend/calibration/runner.py`.
  This protocol is not wired into those entry points.
- The protocol refers to `station_sop.LAW_CONTAINMENT`. That symbol is not
  defined in this candidate, so the protocol's containment language is a
  document-of-record requirement, not an enforced runtime gate here.
- The sibling `split_manifest_template.json` and `cohort_v25_proposal.md`
  are owner-gated proposal documents and were not copied.
- No calibration run artifacts were recovered with this protocol. Any run
  would require a separately pinned exclusion inventory and explicit owner
  authorization.

## Hard stop

No merge, deployment, v25 run, held-out evaluation, threshold or cohort
change, or golden-artifact regeneration is authorized by this recovery. The
C5 published report remains diagnostic-only.
