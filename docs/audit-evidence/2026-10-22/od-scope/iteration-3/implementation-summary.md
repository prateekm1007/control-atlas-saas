# OD-scope iteration 3 - implementation summary

Recorded: 2026-10-22T13:28:53Z
Branch: `review/od-scope-v25`
Base: `63da0e347002ed9b77071ff1ed0cefda7c0ce88a`
Parent / preserved delivery candidate: `9412e943014df6804e34c670d9cc9496e1f34c9c`
Code candidate SHA: `3029e21196dedd7cecd3927598639c9bd1a22e6c`.
Final evidence/archive candidate SHA: pending the evidence commit (the
engine tested here remains the code candidate above).

## Scope

This is a focused follow-up patch on the existing review branch. It
preserves `9412e94` and all prior evidence, and does not rewrite history.
It authorizes no merge, deployment, v25 calibration, held-out evaluation,
threshold or cohort change, or golden-artifact regeneration. OD-16 remains
active.

## Implemented corrections

### Measurement truthfulness

- Zero-observation measurements now report `NOT_MEASURED` instead of
  fabricating an observed zero with `PASS`.
- `NOT_APPLICABLE` is used only where the law cannot apply to the residue
  composition: LAW-145 glycine-only centres, LAW-150 no eligible chi1
  residue types, LAW-195 no cysteine or no disulfide pair, LAW-182 no
  hydrophobic residues.
- `LAW-130` and `LAW-200` report `NOT_MEASURED` for fewer than two atomic
  coordinates; the historical one-atom bounding-box `PASS` is removed.
- The adjudicator normalizes unknown or internally contradictory rows to
  `MEASUREMENT_FAILED` and preserves the normalization reason. Measured
  rows require a finite observation and positive sample; missing rows
  cannot carry an observation; `NOT_APPLICABLE` requires no observation
  and a zero sample.
- The response schema accepts only the declared status vocabulary and
  enforces the same evidence consistency rules.

### LAW-100 and LAW-120

- LAW-100 uses only the canonical percentage of bonds with `|z| > 4 sigma`
  and the `<= 5.0%` decision rule in every modality. The unregistered
  experimental RMSZ path and its `8.0` comparison are removed.
- The registry text and test reference describe the canonical percentage
  metric; the registry digest is re-pinned.
- LAW-120 now measures intra-residue `N-CA-C`, inter-residue
  `O(i)-C(i)-N(i+1)`, and inter-residue `CA-C-N`; the former
  intra-residue `O-C-N` geometry used the wrong sequential mapping.
  The `18.0 deg` threshold is unchanged.

### Confidence contract

- Client-supplied confidence is no longer hard-coded to `90.0`.
- `detect_confidence_source` and `decompose_confidence` derive a mean from
  profile observations or atom B-factors, and return an explicit
  unavailable state when no source distribution exists.
- `confidence_score` and `ConfidenceMeta.mean` are nullable; the versioned
  metadata is `confidence-meta/1.0.0`.
- PDF output labels the top bar as measured LAW-105 residue coverage (not a
  pLDDT distribution) and discloses the confidence source or its absence.
- NKG display renders unavailable confidence as `UNAVAILABLE`, without a
  percent sign.

### Calibration evaluator

- `policy_binding` requires the code SHA, canonical hash, law semantic
  hash, modality-matrix digest, and live registry digest; the running
  engine recomputes and verifies all policy identity fields.
- `run_inputs` requires digests for `run_manifest.json`,
  `coordinate_manifest.json`, and `calibration_results.json`; the optional
  `--run-dir` path verifies those bytes and the exact code/policy binding
  recorded in the run inputs.
- Label provenance and cohort/label/source mappings are required. B1
  entries require a defect-manifest reference and digest.
- `measurement_failed` is a distinct missing-evidence class and is
  counted in missing rates; `PASS` with missing deterministic evidence is
  rejected; `INDETERMINATE` requires at least one missing cause.
- Split manifests are `split-manifest/1.1.0`. Sequence hashes, unique
  coordinate digests, family disjointness, historical exclusions, and
  sequence identity are checked. Identity is independently recomputed with
  the pinned normalized-global-edit-identity method; near-identity pairs
  are rejected before freezing.
- Diagnostic reports explicitly carry `report_type="diagnostic"`,
  `calibration_status="not_approved"`, and
  `supports_threshold_validation=false`.

### Calibration runner

- `v25-runner/1.2.0` validates every redirect hop before requesting it,
  limits redirects to five, and restricts hosts and ports. Credentials in
  URLs and non-HTTPS final URLs are rejected.
- Authoritative AFDB acquisition requires an exact entity, version,
  source URL, PDB format, byte count, and SHA-256 pin. A pin is never
  promoted to a newer runtime version. Incomplete or mismatched pins fail;
  unpinned AFDB inputs are recorded as non-authoritative diagnostics.
- The exact `F#` entity suffix is preserved.
- Registry and law-semantic locks are verified before the output directory
  is created. Existing run directories are refused.
- Run summaries distinguish `not_measured`, `measurement_failed`,
  `unpinned_diagnostic_targets`, `failure_count`,
  `acquisition_failures`, and `measurement_failures`.

### Semantic lock and CI

- `governance/law_semantic.sha256` pins the semantic law payload:
  `d2b3b8501b1c73f67e7839662f70b68596d3c1e25b51d824d6b8a78325816dd5`.
- `governance/law_registry.sha256` pins the registry:
  `81b57f4c68164522b6a27f810da33c2108c1af9682d018ba15c9b8a92a8d0123`.
- The Toscanini gate independently recomputes the semantic lock on the
  exact candidate.
- The evaluator and runner independently recompute the same semantic
  payload, so a registry edit without a re-pin fails closed.

## Validation for this patch

Command (from `backend`, with a local Windows `fcntl` shim because the
platform lacks the module):

```text
python -c "import sys, types, pytest; m=types.ModuleType('fcntl'); m.LOCK_EX=1; m.LOCK_UN=2; m.flock=lambda *a, **k: None; sys.modules['fcntl']=m; sys.exit(pytest.main(['tests','-q','--disable-warnings']))"
```

Observed result at 2026-10-22T13:28Z:

```text
459 passed, 5 skipped, 1262 warnings in 32.13s
```

The five skips are pre-existing platform-specific skips. The suite command
is `pytest tests`; invoking `pytest` over the whole backend collects a nested
`tos/saas/test_key_store.py` with a test-only import path and is not the
supported suite command.

Additional focused result after the last evaluator/runner corrections:

```text
37 passed in 1.77s
```

`git diff --check` reports no whitespace errors.

## Truthful limitations

- v24 has 50 manifest targets; 40 pinned coordinate sets are verified as
  acquired and 10 historical AlphaFold `model_v4` targets are unavailable.
  No placeholder or substituted coordinates are used.
- `AF-P0DTC2-F1` remains excluded under OD-11 and `AF-P62988-F1` remains
  unavailable; both stay visibly unresolved.
- The referenced evaluation protocol file
  `backend/calibration/protocol/evaluation_protocol_v1.md` is not present
  in this branch. No protocol document is being invented in this patch.
- `docs/owner-decisions` remains absent from the live remote. The local
  owner register is read-only and was not edited.
- The C5 diagnostic is generated after this code candidate is committed,
  against base `63da0e3`, using the same verified coordinate bytes and
  digests. It is diagnostic evidence only.
- No v25 calibration, held-out evaluation, threshold/cohort change, merge,
  or deployment was performed.
