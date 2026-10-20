# SECRET_SCAN_BASELINE — EC-23 reviewed fingerprint allowlist (Round-20, 2026-10-20)

This document is the reviewed record behind
`docs/evidence/ec16-ec26/secret_scan_baseline.json`, the machine-readable
allowlist enforced by `scripts/ci/secret_baseline_check.py` in the
Governance Acceptance workflow's `secret-scan-all-refs` job.

**Defect being corrected (auditor Round-20 direction 3).** The Round-19
all-refs classification tolerated every `generic-api-key` finding under
`docs/` or `evidence/` as an "evidence-digest false-positive class". That
category exemption was not sufficiently narrow: a newly introduced credential
committed under those directories would also have passed. The Round-19
positive control additionally treated any nonzero scanner exit as proof of
detection, so a scanner crash could masquerade as a successful control (fixed
by `scripts/ci/secret_control_check.py` and AT-22).

**Replacement policy.** Every finding of the all-refs scan must match an
EXPLICIT, REVIEWED fingerprint in the allowlist. Each entry carries its source
(how and when the finding was observed), its classification, and the reason
the finding is not a live credential. Unknown fingerprints FAIL the scan
wherever the file lives — including `docs/` and `evidence/`. A fingerprint
that matches while its rule or file contradict the reviewed entry also FAILS.
There is no directory-based exemption and no `.gitleaksignore` on this branch.

**Source scan for this revision.** gitleaks `8.30.1` (release binary
SHA-256 `551f6fc83ea457d62a0d98237cbad105af8d557003051f41f3e7ca7b3f2470eb`,
the same checksum the workflow pins), invoked as
`gitleaks detect --source . --log-opts="--all" --redact --report-format json`
over the full fetched remote ref set (all 38 remote branches as of
2026-10-20T15:43Z, before this round's commits). The scan reported **48 rows /
47 distinct fingerprints**, all classified below. Matched-line structures were
inspected with values masked; **no secret value is reproduced in this file or
in the JSON baseline** — findings are identified by gitleaks fingerprint
(`<commit>:<file>:<rule>:<line>`) only.

## Classification summary

| Class | Fingerprints | Meaning |
|---|---|---|
| DOCUMENTED-HISTORICAL-A-02 | 2 | the former provider key in pre-release history — compromised until issuer-side revocation is verified; tolerated ONLY while register A-02 stays open; deployment blocked; history never rewritten as a substitute for revocation |
| FALSE-POSITIVE-PROVENANCE-DIGEST | 41 | matched text is a recorded SHA-256 content digest of a public API response (`"api_response_sha256"` / `"api_pdb_digest"` fields in evidence artifacts) — structural false positive, not a credential |
| FALSE-POSITIVE-IDENTIFIER-LIST | 1 | matched text is an internal task-identifier list in the 2026-10-19 preflight record |
| FALSE-POSITIVE-DOCUMENTATION-QUOTE | 2 | historical revisions of the Round-18 findings-disposition document quoted scanner fingerprints in their own classification tables, re-creating the structural match; those files are absent from the current tree and history is never rewritten |
| FALSE-POSITIVE-PUSH-LOG-EVIDENCE-QUOTE | 1 | a historical push-log confirmation (concurrent writer's activation branch) cites `commits API: <40-hex>` — a git commit SHA quoted after that phrase, not a credential |

## DOCUMENTED-HISTORICAL-A-02 (disclosed on every scan run; A-02 open)

| Fingerprint | Rule | File:line | Basis |
|---|---|---|---|
| `b66f0e3bc8b5e86f5810ff52b598d917420ac7a3:probe_models.py:gcp-api-key:5` | gcp-api-key | `probe_models.py:5` | A-02 former provider key; issuer-side revocation NOT verified; tree removal (H-1 branch) is not revocation; deployment stays blocked |
| `b66f0e3bc8b5e86f5810ff52b598d917420ac7a3:system_diagnosis.py:gcp-api-key:8` | gcp-api-key | `system_diagnosis.py:8` | same A-02 class; both entries are pre-release history, immutable by policy |

## FALSE-POSITIVE-PROVENANCE-DIGEST (41 fingerprints)

40 fingerprints on
`docs/audit-evidence/2026-10-20/modres-census-v2/n7_od16_40target_table.json`
(commit `dea6826b9e9649cdc8ae4e2dabda4647fe40ece0`, lines 82–1869: the
`"api_response_sha256"` field of every census row) and 1 fingerprint on
`evidence/incident-push/raw/hardstop_evidence_snapshot.json` (commit
`6b46988940723ae173d3cc9a2905822d9fc4c0f4`, line 636: the `"api_pdb_digest"`
field). Basis for all 41: the matched text is a 64-hex SHA-256 content digest
of a public API response recorded as provenance under the OD-07 digest
discipline — a digest of public data, not a credential; no issuer exists. The
full per-line fingerprint list is machine-readable in
`secret_scan_baseline.json` (entries with this classification).

## FALSE-POSITIVE-IDENTIFIER-LIST (1 fingerprint)

| Fingerprint | Rule | File:line | Basis |
|---|---|---|---|
| `dfdaacf259febea9a65dce63040c5ef2e2113901:docs/audit-evidence/2026-10-19/preflight.json:generic-api-key:5` | generic-api-key | `docs/audit-evidence/2026-10-19/preflight.json:5` | matched text is the internal task-identifier list (`"task_ids": ["P0-credential-containment", ...]`) of the 2026-10-19 preflight record — identifiers, not credentials |

## FALSE-POSITIVE-DOCUMENTATION-QUOTE (2 fingerprints)

| Fingerprint | Rule | File:line | Basis |
|---|---|---|---|
| `f29ea8c18dab3b104489a9fa0f0ec5f0aa252c83:docs/evidence/ec16-ec26/GITLEAKS_FINDINGS_DISPOSITION.md:generic-api-key:77` | generic-api-key | (historical revision, file absent from the current tree) | the Round-18 disposition document quoted scanner fingerprints in its classification table, re-creating the structural match; preserved as history, never rewritten |
| `9da33ef48a19fd54ed7072df0868aa216480982e:docs/evidence/ec16-ec26/GITLEAKS_FINDINGS_DISPOSITION.md:generic-api-key:76` | generic-api-key | (historical revision, reported twice by the scanner — one distinct fingerprint) | same class, earlier revision of the same table |

## FALSE-POSITIVE-PUSH-LOG-EVIDENCE-QUOTE (1 fingerprint)

| Fingerprint | Rule | File:line | Basis |
|---|---|---|---|
| `af713557a49574df659b11fefa0d22d3a9e165ce:docs/evidence/push-log/PUSH_LOG.md:generic-api-key:361` | generic-api-key | `docs/evidence/push-log/PUSH_LOG.md:361` (historical revision on branch `activation/ec16-ec20-consolidated`) | the concurrent writer's `P-034-CONFIRM` evidence field cites `commits API: <40-hex>` — a git commit SHA quoted after that phrase; a commit identity, not a credential |

## Maintenance discipline

- A new finding anywhere in the ref set (including `docs/` and `evidence/`)
  fails the scan until a reviewed entry with source, classification, and basis
  is added to the JSON baseline through a reviewed change — never by editing
  this file alone, never by a path exemption.
- A fingerprint whose finding disappears from history (branch deletion,
  history replacement) should be removed by a reviewed change; stale entries
  are reported by the checker as unused and reviewed at each activation.
- The two A-02 entries leave the allowlist only when register A-02 closes
  with issuer-side revocation evidence AND the pre-release history containing
  them is deliberately rewritten by an owner decision — neither condition
  holds today, so they stay, disclosed on every run.
