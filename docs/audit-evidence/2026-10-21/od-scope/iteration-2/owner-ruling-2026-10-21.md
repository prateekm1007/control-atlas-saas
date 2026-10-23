# Owner ruling recorded for the OD-scope review candidate

Recorded: 2026-10-21T21:13:45Z
Branch: `review/od-scope-v25`
Base: `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
Canonical register reviewed: `docs/decisions/OWNER_DECISIONS.md` blob
`12906bf5997a3e2c6b5c41e5d5c60fee14613274`, SHA-256
`59ccd5e84176a8a610c8c234cbb95fec4734478cc1f35668599ad2e3aef3172a`.

## Controlling decision

The canonical register text controls. The historical local clarification at
`evidence/engine/owner-clarification-od-03-04-20.md` on
`engine/od-scope-v25` (`97d62a8`) is preserved as evidence and is **not** an
amendment.

- `OD-03`: pLDDT is advisory only for `LAW-130` and `LAW-200`; it does not
  independently create a favorable verdict. `LAW-155` is retired from
  enforced semantics. No new pLDDT threshold is authorized.
- `OD-04`: do not invent a disorder threshold or let low-confidence/missing
  measurements create a favorable result. Preserve truthful
  `NOT_MEASURED`/`MEASUREMENT_FAILED`/`INDETERMINATE` states.
- `OD-20`: keeps NMR separate from resolution-based strata; it is not
  authority for a `LAW-182` policy change.

## Authorized candidate scope

- Keep the `LAW-105` residue-coverage gate and total-residue scope. Do not
  return or gate on atom coverage.
- Keep `LAW-110` on all residues in original order.
- Do not apply a new pLDDT core scope to `LAW-130`, `LAW-155`, or `LAW-200`.
  `LAW-130` and `LAW-200` may expose pLDDT as advisory metadata only.
- Retain the `LAW-155` ID for traceability, mark it retired/not enforced, and
  exclude it from verdict and scoring semantics.
- Keep `LAW-182` classification unchanged; test that advisory/heuristic
  signals cannot independently veto.
- No new decision IDs and no edit to the canonical register.

## Delivery rule

The pre-ruling staged tree is not the delivery candidate. After the required
engine and calibration-tool fixes and tests pass, create one full-SHA commit
on this isolated review branch and push only that feature branch, without
force-push or changes to `main`. Do not merge, deploy, run v25 calibration,
run held-out evaluation, change thresholds/cohorts, or regenerate golden
artifacts.
