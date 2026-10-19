# DELIVERY README — Round-11 (U1–U9) evidence outputs

**Delivery route.** This `outputs/` folder is pushed to GitHub on the new
branch `evidence/round9-r2-2026-10-20` (base: `origin/main` @
`2195a048…`), explicitly authorized by the owner for this round. It is
the durable evidence path; the coder's workspace `download/` folder is
not owner-reachable. Nothing else was pushed: no existing branch, tag,
or ref moved; nothing merged or deployed; PR #16 remains open/unmerged.

## SHA-256 values (repeat of the report)

| Artifact | SHA-256 | Size |
|---|---|---|
| r1 evidence archive (sealed, unchanged) | `42d8e8819bc7024bfad1d9a511c7707046c3c2a0b3ff38a07fc49d4861297528` | 22,560,251 B |
| **r2 evidence archive (this delivery, byte-identical to the Round-10 build)** | `5ac52b9332c2fe2c2510b443c6aeb26ce785eccf61ad6d44cebce6151a3257cd` | 22,571,130 B |

The r2 archive was NOT rebuilt this round: the U3 TIPS-array test passed
on the sealed bytes, so U3's rebuild condition ("if any key is wrong")
was not triggered. r1 and r2 both remain byte-identical to their
Round-9/Round-10 builds.

## Contents

| Path | What it is | SHA-256 |
|---|---|---|
| `control-atlas-evidence-2026-10-20-round9-r2.tar.gz` | The r2 evidence archive (186 tar members; constitution bytes, OD-16/OD-17 cohort audit with 101 raw RCSB API responses, N7 comparison, register rev 3, N3 count correction, 6 verified git bundles of the pinned remote tips, 33 patches, lost-work register, scripts, verifier) | `5ac52b93…` (full value above) |
| `verifier.sh` | The r2 verifier, byte-identical to the copy inside the archive (SHA-256 `b05972de3f49ebfb95e5f85d653f206af591f73979b2114c608c40445e95880d`). Extract the tar.gz and run it from the extraction directory: `bash verifier.sh` (offline) or `bash verifier.sh --network` | `b05972de…` |
| `OWNER_DECISIONS.md` | **The register** — consolidated rev 4: the auditor's OD-01…OD-21 table verbatim, all statuses OPEN, recommendations in a separate advisory field, the 20-of-40 correction attached to OD-02 and OD-14 (U7), OD-21 held (U8) | `b68faee9…` |
| `INTERNAL_CROSSWALK.md` | INTERNAL annex — the only file carrying retired D-xx labels and the numbering-collision history (U5) | `f36e133c…` |
| `READING_NOTE_session3_2026-10-20.md` | U1 reading note: both constitution blob citations, every EC §1–§15 + OBS-153 + CRIT-154 ID enumerated with enforcement status, counting conventions, fetch/read/edit sequence timestamps | `401b8b43…` |
| `MODRES_CENSUS.md` | U9 census documentation: 40/40 targets digest-bound; **0 modified residues in coordinates**; 11 MODRES annotation records in 5A63 characterized (glycosylation-site annotations on standard ASN) | `d6b85e07…` |
| `evidence/verifier-tips/` | U3 closure: `cat -A` dumps of both sealed verifiers, `od -c` byte proofs of the `[main.bundle]=` line, the six-key TIPS test (`test_tips_array.sh`) and its passing runs against BOTH sealed archives, and a full verifier pass on an independent r2 extraction | see `SHA256SUMS` |
| `evidence/od-register-audit/` | U5 audit: D-xx retirement + OD single-meaning grep audit (PASSED) | see `SHA256SUMS` |
| `evidence/modres-census/` | census script, run log, machine-readable per-target results (with digests) | see `SHA256SUMS` |
| `evidence/constitution-counts/` | session-3 count re-derivation script + run log (both constitutions MATCH the pinned counts under the single pinned convention) | see `SHA256SUMS` |
| `patches-runner-provenance/` | U9: runner provenance patch (cover + 0001, base `science/p0-correctness-2026-10-19` @ `ebc1fd79`) + the **demonstration manifest** produced by running the committed patch itself (`code_sha f2a87fe…`, clean tree, sanitized remotes) | see `SHA256SUMS` |
| `patches-n8-per-residue-altlocs/` | U9/U4: per-residue altloc (N8) feature patch (cover + 0001, base `ebc1fd79`) — a FRESH implementation; the lost tip `e9985c6` is not restored or claimed | see `SHA256SUMS` |
| `SHA256SUMS` | SHA-256 of every file in this folder (except itself); verify with `sha256sum -c SHA256SUMS` | — |

## U1–U9 completion status

- **U1 (constitution reading):** DONE — both documents read in full from
  git objects at the pinned blobs after a fresh fetch (2026-10-19T21:48:59Z,
  no ref movement, no new remote branch); reading note cites the blobs and
  enumerates every EC/OBS/CRIT ID with enforcement status.
- **U2 (durable delivery):** DONE — this folder on the pushed evidence
  branch; SHA-256 values repeated above and in `SHA256SUMS`.
- **U3 (TIPS array):** CLOSED — `cat -A` + `od -c` byte evidence and a
  passing six-key test on BOTH sealed archives; the `ain.bundle]=`
  sighting is a display-layer artifact of the `[m` byte pair; no key was
  wrong; r2 not rebuilt, r1 untouched.
- **U4 (unpushed work as patches):** PARTIAL BY NECESSITY — the docs
  branch tip `94b0f7c` and N8 tip `e9985c6` are **not present in this
  environment** (lost with the prior sandbox reset; see the lost-work
  register). Bundles/patch series for those exact tips cannot be produced
  from available objects, and none were fabricated. What is delivered:
  the register content (this branch), a fresh N8 feature patch, and the
  six verified bundles of the pinned REMOTE tips inside the r2 archive.
  Both lost tips remain absent from GitHub (verifiable via `git ls-remote`).
- **U5 (one register):** DONE — consolidated rev 4; D-xx retired to the
  internal crosswalk; grep audit PASSED (each OD ID one meaning).
- **U6 (all OPEN):** DONE — every decision OPEN; DECIDED requires the
  owner's written text, name, and date; recommendations live only in the
  advisory field.
- **U7 (20-of-40 correction into decisions):** DONE — OD-02 and OD-14
  carry the corrected context (artifact blob `6cda70e2…`: 0 of 40, 38
  LAW-120@10°, 2 LAW-170 (1I0T, 1DDL); the "at most 20 of 40" figure was
  the auditor's un-recomputed bound and is labeled as such).
- **U8 (word target held):** DONE — OD-21 OPEN, no number chosen, no
  document length changed.
- **U9 (unblocked work):** DONE — per-residue altlocs on a feature patch
  (default policy preserved; owner-gated to flip); MODRES census recorded
  as documentation (with the 5A63 annotation finding recorded, not
  smoothed); runner provenance as a patch with a demonstration manifest.

## What was NOT done (owner-gated, unchanged)

No v25 run, no held-out evaluation, no threshold or cohort changes, no
merges, no deployments, no pushes beyond this evidence branch, no
golden-artifact changes. The two code patches are patch-files only —
`git am` them onto `science/p0-correctness-2026-10-19` if the owner
chooses; nothing is applied or merged here.

## Credential note

The PAT used for the one push of this branch exists only in a transient
command line. It does not appear in any file in this delivery (checkable
against `SHA256SUMS`), in the archives, in the demonstration manifest
(remotes are recorded credential-stripped), or in any patch. Per the
owner's instruction this round, no revocation reminders are issued here;
register entry A-01 continues to track the action itself.

## How to verify this delivery

```bash
git fetch origin && git checkout evidence/round9-r2-2026-10-20
cd outputs
sha256sum -c SHA256SUMS                      # every file
sha256sum control-atlas-evidence-2026-10-20-round9-r2.tar.gz
#   -> 5ac52b9332c2fe2c2510b443c6aeb26ce785eccf61ad6d44cebce6151a3257cd
mkdir /tmp/r2x && tar -xzf control-atlas-evidence-2026-10-20-round9-r2.tar.gz -C /tmp/r2x
cd /tmp/r2x && bash verifier.sh              # offline: ALL CHECKS PASSED
bash evidence/verifier-tips/test_tips_array.sh /tmp/r2x/verifier.sh /tmp/r2x
#   -> TIPS ARRAY TEST PASSED (six keys, pinned tips)
git ls-remote origin | grep -cE '94b0f7c|e9985c6'   # -> 0 (lost tips absent)
```
