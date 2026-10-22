# Owner ruling recorded for LAW-100 in iteration 3

Recorded: 2026-10-22T12:41:30Z
Branch: `review/od-scope-v25`
Base: `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
Candidate at recording: `9412e943014df6804e34c670d9cc9496e1f34c9c`
Canonical register reviewed: `docs/decisions/OWNER_DECISIONS.md` blob
`12906bf5997a3e2c6b5c41e5d5c60fee14613274`, SHA-256
`59ccd5e84176a8a610c8c234cbb95fec4734478cc1f35668599ad2e3aef3172a`.

## Controlling ruling

`LAW_CANON` in `backend/tos/governance/station_sop.py` is the runtime authority
for LAW-100. The canonical metric is:

- the percentage of covalent bonds whose length deviates more than 4 sigma from
  the Engh-Huber ideal, across the declared core-residue scope;
- threshold `5.0`, unit `%`, operator `<=`, type `Percentage`.

The candidate's experimental-branch RMSZ value and its unregistered `8.0`
comparison are not a canonical LAW-100 metric. They must be removed. The
measurement result must use the same percentage metric, `observed` field, and
`<= 5.0%` decision rule for all source regimes, subject to the existing
modality method classification.

No new threshold is created or approved by this ruling. No cohort or data
selection is changed. The existing modality classification remains unchanged.
No new OD identifier is created, and the canonical decision register is not
edited.

## Required synchronization

- Measurement code: canonical percentage only.
- Registry and glossary: describe the percentage measurement, not RMSZ.
- Tests and semantic locks: assert the canonical metric and re-pin the registry
  hash after the synchronized documentation change.
- Audit evidence: record that verdict differences caused by the removed
  experimental branch are diagnostic only. They are not validation of a new
  threshold and do not authorize v25 calibration.

## Delivery boundary

This ruling does not authorize merge, deployment, v25 calibration, held-out
evaluation, threshold tuning, cohort changes, or golden-artifact regeneration.
OD-16 remains active.
