# Constitution Reading Note — Session 3 (Round-11, auditor directions U1–U9)

**Date:** 2026-10-20 (session-local, UTC+05:30) / 2026-10-19 21:48–21:53 UTC
**Read by:** coder (this session, end-to-end, from git objects — not from
working-tree copies)

## 1. Documents read, pinned identities, and sequence evidence

Both constitutions were read in full this session, dumped directly from git
objects at the audited remote tips after a fresh fetch:

| Document | Read from | Commit | Blob SHA-1 (pinned) | Bytes | Chars | Words | Lines (LF/split/non-blank) |
|---|---|---|---|---|---|---|---|
| `constitution.md` (Toscanini MDI Constitution v6.1-Hardened, domain) | `origin/main` | `2195a048e87dcec4aa63e8dc7981b3bb6d237ca0` | `211ad2c5d4b40059d9f3f03de12e0701f9235bc7` | 622 | 615 | 84 | 11 / 12 / 9 |
| `ENGINEERING_CONSTITUTION.md` (SaaS Engineering Constitution v1.0.0, proposed canonical) | `origin/policy/engineering-constitution` | `4993cb3b6954ffb8ecad281817b2c59f1f8fc041` | `b6bb35e84956d84a8c9d9e65734c212baf3ed075` | 20572 | 20564 | 2784 | 178 / 179 / 129 |

Counting convention (single pinned convention, T1): **bytes** = octets;
**chars** = UTF-8 decoded length; **words** = whitespace-split
(`len(text.split())`, `wc -w` equivalent); **lines** = LF count,
`split('\n')` element count, and non-blank-line count, reported separately.
Locale pinned `LC_ALL=C.UTF-8`, `LANG=C.UTF-8`. Re-derived this session by
`scripts/session3_constitution_counts.py` from the git-object dumps
(`work/session3-read/`), byte-verified with `git hash-object` (both MATCH
the pinned blobs above).

**Sequence evidence (read precedes every edit this round):**

- Fresh `git fetch origin` — started 2026-10-19T21:48:59Z, exit 0 at
  21:49:00Z. No ref movement on any audited branch; no new remote branch
  (20 heads, same set as the Round-10 record).
- Constitution dumps + `git hash-object` verification — 21:49:08Z.
- Full read of both documents (end-to-end) — completed 21:49–21:50Z.
- Count re-derivation — 21:52:12Z.
- **First content edit of this round** (staging tree creation for the
  outputs delivery) — 21:52:22Z. Nothing under `outputs/`, the register,
  or any evidence document was written before the reads above completed.

## 2. Every ID this session reads and enforces

The domain constitution carries two ACTIVE normative IDs; the engineering
constitution carries fifteen numbered sections plus a final rule. The
citation convention used in the decision register ("EC §n") refers to these
sections. **All** of them are enumerated below — nothing is skipped — with
how this round's deliverables enforce each.

### 2.1 Domain constitution IDs (both ACTIVE)

| ID | Title | Content (one line) | Enforced this round |
|---|---|---|---|
| **OBS-153** | Avidity Quantification | Bipodal success measured by contact density ρ; requires ρ(Bipod) ≥ 2.5 × ρ(Monopod), unique heavy-atom contacts within 4.5 Å | Not exercised (no avidity measurement this round). Enforced in the narrow sense: its text and thresholds are preserved byte-identical (blob verified above) and no threshold, ID, or meaning is altered by any deliverable here (EC §13). |
| **CRIT-154** | Structural Admissibility | Lead validation requires two-stage warhead-RMSD < 1.0 Å against clinical benchmark PDB 5NIU; veto on non-physical target deformation (steric clash < 1.2 Å) | Not exercised (no lead-validation measurement this round). Same preservation enforcement as OBS-153. |

Honesty note: neither ID is *measured* by this round's work; claiming
"enforced" beyond preservation-of-definition would overstate. Their ACTIVE
status is recorded and unchanged.

### 2.2 Engineering constitution sections (EC §1–§15 + Final Rule)

| ID | Subject | Enforced this round |
|---|---|---|
| **EC §1** Purpose and authority | Engineering constitution governs engineering/security/product integrity; `constitution.md` governs the scientific claims it actually states; coder MUST read both before coding; conflicts must be recorded, not silently resolved; fail closed on ambiguity | YES — both documents read in full before any edit (§1 timestamps above); the OD numbering conflict and the 20-of-40 correction are recorded, not smoothed over. |
| **EC §2** Mandatory pre-coding protocol | Confirm repo/branch/base SHA; read both constitutions, current implementation, contracts, tests, callers before editing; record expected behavior/risks/acceptance | YES — base SHAs recorded (§1); for each code patch (N8, runner provenance) the current implementation and its tests were read from the science-branch tip before writing the change. |
| **EC §3** Security is a release requirement | Fail closed; deny by default; authorize objects; **no live credential committed to source, fixture, image, report, or logs; a credential in Git history is compromised until its issuer confirms revocation**; validate inputs; protect callbacks/data | YES — the PAT used for the evidence push exists only in transient command lines, never in any committed file, manifest, archive, or log; the demo provenance manifest sanitizes git remote URLs of userinfo; A-01/A-02/A-03 remain recorded owner actions. |
| **EC §4** Product and scientific integrity | Provenance binding for authoritative results; explicit non-success states; measured vs inferred never conflated; law-registry discipline; no silent renumbering of law IDs | YES — evidence documents keep VERIFIED / CODER_REPORTED / NOT_REPRODUCIBLE distinctions; the 0-of-40 correction is labeled historical; the MODRES census binds every counted file to its recorded SHA-256. |
| **EC §5** Architecture and maintainability | Minimal supported workflow; defined component contracts; no duplicate sources of truth; smallest complete change | YES — the N8 patch adds one policy with tests, no unrelated changes; the provenance patch extends the existing runner rather than adding a parallel mechanism. |
| **EC §6** Testing and evidence | Tests deterministic and isolated; claimed results include exact command, commit SHA, runtime versions, and counts; never invent counts | YES — every check in this delivery (TIPS test, census, verifier) records its exact command and output; counts in this note are re-derived, not copied. |
| **EC §7** CI/CD and supply chain | Feature branches + PRs; direct pushes to main disabled; checks on the exact merge candidate | RESPECTED — no push to `main`, no merge, no existing ref moved. The one push this round is a NEW evidence-delivery branch (`evidence/round9-r2-2026-10-20`) carrying `outputs/` only, explicitly owner-authorized in the Round-11 direction; it changes no product code path. |
| **EC §8** Reliability, observability, operations | SLIs/SLOs, bounded timeouts, structured logs, runbooks | Not exercised (no deployment this round). |
| **EC §9** Customer experience, accessibility, cost | Clear workflow communication; abuse controls on costly endpoints | Not exercised. |
| **EC §10** LLM and external-provider boundaries | Narrative advisory; authoritative results independent; provider failure = honest degraded state | Not exercised. |
| **EC §11** Change, release, incident discipline | Every change has problem statement/criteria/risk/tests/rollback; **never erase audit evidence or rewrite history** | YES — r1 archive preserved byte-identical; r2 archive preserved byte-identical (this round's materials ship beside it, not inside it); lost-work records kept, not rewritten. |
| **EC §12** Definition of done | Done only with recorded SHAs, passing checks on the proposed commit, docs updated; "not run" never reported as "passed" | YES — the delivery README separates run / not-run per item. |
| **EC §13** Prohibited shortcuts | No coding before the protocol; no claimed testing without evidence; no gate weakening for green; **no secrets in reports, artifacts, logs, prompts, or issue text**; no mock data as production evidence; no unapproved scientific-semantic changes; no history rewrites | YES — the PAT rule is applied to every artifact; no engine threshold, law meaning, or published claim is altered by this round's patches; per-residue altloc ships as an opt-in policy with the default preserved. |
| **EC §14** Required task artifacts | Preflight record + baseline SHA; scope/dependency map; machine-readable results bound to commit; changed-file manifest + patch hash; impact note; final checklist | YES — the DELIVERY_README + per-patch provenance blocks provide these for each patch; the r2 archive's MANIFEST covers the bundle. |
| **EC §15** External engineering references | SSDF / ASVS / SLSA as references; no compliance or certification claims | YES — no certification language anywhere in this delivery. |
| **Final Rule** | Build trust before features; prove behavior before claiming success; preserve evidence before changing authority; automate every checkable invariant | YES — this round's entire shape (verify-then-deliver, sealed archives, tests-as-evidence) follows it. |

Sections not exercised this round (§8, §9, §10) are recorded as such rather
than claimed.

## 3. What this round's work is authorized to do

Per the Round-11 directions: constitution reading with ID citation (U1);
durable `outputs/` delivery of the r2 evidence package (U2); the TIPS-array
closure with `cat -A` evidence and a six-key test (U3); honest handling of
the unrecoverable docs/N8 branch tips (U4); adoption of the consolidated
OD-01…OD-21 register with all statuses OPEN (U5, U6); the 20-of-40
correction attached to OD-02 and OD-14 as context (U7); the word target held
(OD-21 OPEN, no length change) (U8); and the three unblocked work items as
patches — per-residue altlocs (N8) on a feature patch, the MODRES census
recorded as documentation, and runner provenance as a patch with a
demonstration manifest (U9). Everything else remains owner-gated.
