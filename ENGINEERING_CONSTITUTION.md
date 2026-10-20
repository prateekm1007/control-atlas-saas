# SaaS Engineering Constitution

**Status:** Proposed canonical engineering policy  
**Version:** 1.0.0  
**Applies to:** `prateekm1007/control-atlas-saas` and its production services  
**Companion domain policy:** `constitution.md`  
**Change policy:** Amendments require an explicit, reviewed change and must never be smuggled in through an implementation PR.

## 1. Purpose and authority

This Constitution defines the minimum guarantees required for every product change. The goal is a secure, reliable, maintainable, testable, explainable SaaS—not merely code that appears to work on a happy path.

The existing `constitution.md` governs the scientific/domain claims it actually states. This document governs engineering, security, product integrity, testing, operations, and release behavior. The coder MUST read both before coding. The coder MUST also read the current implementation, contracts, tests, and deployment configuration affected by the proposed change.

When a claim in a README, handoff, glossary, registry, test, paper, API schema, or implementation conflicts with another source, the coder MUST NOT choose whichever source is easiest to edit. They MUST record the conflict, identify the effective runtime path and evidence, preserve existing behavior unless a deliberate decision authorizes a change, and fail closed where ambiguity could compromise security, data integrity, a customer decision, or a scientific verdict.

The words **MUST**, **MUST NOT**, **SHOULD**, and **MAY** express requirements, recommendations, and options in their ordinary normative sense.

## 2. Mandatory pre-coding protocol

Before editing any file, the coder MUST:

1. Confirm the repository, branch, working-tree state, and exact base commit SHA. Never assume `main` has not moved.
2. Read this Constitution and `constitution.md` in full. Record their paths, versions or hashes, and the rules relevant to the task.
3. Read the current README and handoff, architecture/contracts, target implementation files, related tests, dependencies, deployment/build files, and every caller or consumer affected by the change.
4. Map the change across inputs, trust boundaries, authorization, data persistence, asynchronous jobs, output contracts, observability, and deployment. Search for duplicate implementations and alternate entry points before modifying a behavior.
5. Record the expected behavior, relevant risks, files expected to change, regression tests, rollout/rollback strategy, and acceptance criteria in a machine-readable preflight artifact attached to the task/PR or CI run.
6. Identify contradictions, missing evidence, stale docs, unclear ownership, and untested behavior before coding. If a conflict could change security posture, customer-visible results, scientific meaning, data retention, or an external contract, stop that part of the change and report the evidence rather than guessing.

A preflight attestation is not proof that the document was understood. Automated checks MUST enforce the technical invariants below; a checkbox or statement that the Constitution was read is never a substitute for passing tests and gates.

## 3. Security is a release requirement

Security is part of correctness. A feature that works only by weakening authentication, privacy, isolation, or validation is not complete.

- **Fail closed.** Production MUST NOT default to open authentication, development secrets, permissive CORS, disabled validation, or insecure debug behavior. Missing/invalid security configuration MUST produce a clear startup or request failure.
- **Deny by default.** Every API route MUST have an explicit access policy. Unless explicitly classified and tested as public, a route MUST require an authenticated user or authenticated service identity. Service callbacks may use an appropriately signed, expiring, audience-bound service credential; they MUST NOT trust caller-supplied JSON, user IDs, audit IDs, or network location alone.
- **Authorize objects, not just requests.** Authentication alone is insufficient. Every read, write, job, certificate, history record, audit, and coordinate artifact MUST be scoped to the authorized tenant/owner or to a specifically authorized service role. Cross-tenant access MUST be tested and denied.
- **Manage secrets correctly.** No live credential, private key, password, token, or production secret may be committed to source, a fixture, an image, a generated report, or logs. Use managed secret storage or environment injection with rotation. Treat a credential found in Git history as compromised until its issuer confirms revocation/rotation; deleting the text alone does not revoke it. Secret-scan the working tree and, where appropriate, history without printing secret values.
- **Validate inputs at boundaries.** Validate type, size, format, encoding, content, and allowed values before expensive work. Protect archive processing against traversal, excessive expansion, oversized members, and resource-exhaustion attacks. Apply rate limits, quotas, concurrency bounds, and timeouts to costly work.
- **Protect callbacks and state transitions.** Use signature or token verification, expiry, replay protection, authorization, and idempotency where relevant. State-changing operations MUST be safe under retries and concurrent requests.
- **Protect sensitive data.** Treat user-submitted protein coordinates, structures, identifiers, history, audit reports, and generated artifacts as potentially proprietary. Minimize collection and exposure, enforce retention/deletion policies, encrypt transport and supported persistent storage, and keep sensitive payloads out of logs and public error responses.
- **Use security verification standards as measurable guidance.** Map applicable requirements to an explicitly chosen OWASP ASVS level and maintain evidence for each applicable control. Do not claim formal compliance or certification without an appropriate assessment.

## 4. Product and scientific integrity

The product MUST distinguish what was measured from what was inferred, missing, simulated, or failed.

- Every authoritative result MUST have provenance binding it to the input/content digest, engine/version, applicable law-registry version/hash, measurements, configuration, and run status.
- Acquisition failure, invalid input, empty input, unsupported format, incomplete measurement, failed simulation, timeout, and persistence failure MUST be represented as explicit non-success states. They MUST NOT be converted into a PASS, a VETO, a zero-valued measurement, or a plausible-looking score unless the semantics of that result explicitly justify it.
- Synthetic/mock data MUST be test-only by default. It MUST NOT receive an authoritative production verdict, enter the production negative-knowledge store as observed evidence, or be presented as real acquired material.
- Every measurement MUST distinguish `not measured` from `measured and passed`, `measured and failed`, and `measurement failed`. An absent measurement MUST NOT silently become a favorable zero or a scientific failure.
- The system MUST distinguish deterministic checks, heuristics, advisory signals, ML estimates, and user-provided claims. A heuristic MUST NOT be described as a deterministic physical law. Narrative/LLM output MUST NOT override the authoritative engine verdict or invent evidence.
- Pre- and post-refinement/minimization measurements and verdicts MUST be stored and presented separately. A failed simulation MUST NOT reuse the pre-run result as if post-run validation succeeded.
- A law registry MUST declare the ID, title, principle, class/tier, active status, threshold, operator, unit, scope, metric, implementation, tests, citations, provenance, and known conflicts. Only executable, tested entries explicitly marked as enforced may be described as enforced. Declared, advisory, retired, or unmeasured entries MUST be visibly distinguished.
- Law IDs cited in papers MUST NOT be silently renumbered or redefined. Conflicts between the runtime canon, JSON, glossaries, legacy ledgers, handoffs, tests, and papers MUST remain explicit until deliberately resolved. The canonical semantic hash MUST cover meanings and thresholds, not only the list of IDs.

## 5. Architecture and maintainability

- Maintain a documented, minimal supported customer workflow. For the current product, prioritize validated upload/acquisition → authoritative audit → provenance-backed explainable verdict → certificate/export.
- A component MUST have a defined responsibility and contract. Avoid circular imports, hidden cross-layer side effects, duplicated business rules, unowned background work, and multiple competing sources of truth.
- Before deleting, moving, merging, or duplicating code, automatically inspect import graphs, runtime entry points, Docker build contexts, workflows, tests, and deployment references. Similar names or identical content do not prove two modules are interchangeable.
- API schemas, persistent records, event names, and public response semantics are contracts. Changes MUST be versioned or backward compatible, or deliberately migrated with tested consumer updates.
- Configuration MUST be centralized, validated at startup, documented, and safely defaulted. Behavior must not vary silently with the current working directory or a local developer's machine.
- Dependencies MUST be justified, bounded/pinned to a suitable reproducibility policy, scanned for vulnerabilities, and installed from the correct component manifest. Backend and dashboard dependencies MUST not be confused.
- Database/schema changes MUST have tested migrations and a rollback or forward-recovery plan. Do not perform destructive migrations or broad data deletion as incidental cleanup.
- Prefer the smallest complete change that solves a verified problem. Scope expansion and broad rewrites require a separate justification and regression plan.

## 6. Testing and evidence

A change is not proven by a successful import, a green unit test alone, or a manually observed happy path.

- Reproduce a bug with a failing regression test before fixing it whenever feasible.
- Cover unit behavior, integration boundaries, API/response contracts, authorization and tenant isolation, invalid and hostile inputs, concurrency/replay where relevant, dependency failure, timeouts, persistence failure, and production error paths.
- Tests MUST be deterministic, isolated, repeatable in CI, and must not require live production credentials or mutate production data. External services MUST be replaced with explicit controlled test doubles in tests only.
- A test must assert the intended contract and meaningful failure behavior. Do not weaken, delete, skip, xfail, or hide a test merely to obtain green CI. Every skipped/xfail test needs a reason and owner or resolution condition.
- A claimed test result MUST include the exact command, commit SHA, dependency/runtime versions where material, and pass/fail/skip/error counts. Never invent counts or state that tests passed when they were not run.
- Tests are necessary but not sufficient. Review security boundaries, scientific validity, migration safety, and production behavior separately.

## 7. CI/CD and software supply chain

- All changes MUST be developed on a feature branch and reviewed through a pull request. Direct pushes to `main` MUST be disabled for normal development.
- Required checks MUST run on the exact proposed merge candidate and must pass before merge. A check on a different SHA is not evidence for the release commit.
- Deployment MUST depend on successful tests, validation, and build for the same source SHA or a verifiably equivalent immutable artifact. An independent push-triggered deployment MUST NOT bypass a failed quality gate.
- CI MUST include relevant unit/integration/contract/security regressions; formatting, lint/type checks appropriate to the codebase; secret detection; dependency vulnerability checks; container/build checks; law-registry and generated-artifact drift checks; and migration validation when applicable.
- Production artifacts MUST be traceable to source revision and build process. Adopt SLSA provenance and verification practices incrementally where supported. Do not claim a SLSA level until its required controls have been independently verified.
- Pin third-party CI actions and critical tooling to immutable versions/SHAs according to the repository's maintenance policy. Minimize credentials available to untrusted pull-request jobs and build steps.
- Releases MUST have a tested health/readiness check, failure signal, and documented rollback or recovery path. A deployment command succeeding is not proof the service is healthy.

## 8. Reliability, observability, and operations

- Define service-level indicators and realistic service-level objectives for the supported customer workflow. Measure availability, latency, error rate, job completion, and dependency failures before making reliability claims.
- Use bounded timeouts, controlled retries with backoff, idempotency keys, concurrency/resource limits, and cancellation for expensive or asynchronous operations.
- Health, readiness, dependency health, liveness, and business success are different signals; do not use a superficial HTTP 200 as proof that core functionality works.
- Emit structured logs, metrics, and traces with stable correlation IDs. Logs MUST not expose credentials, tokens, secrets, private structures, or unnecessary personal data.
- Persist authoritative results transactionally where required. Optional telemetry failure may be non-fatal only if it is truly optional and its failure is observable. Required persistence failure MUST not be represented as a completed durable audit.
- Every production-critical workflow MUST have actionable alerts and a runbook for failure, rollback, credential rotation, restore, and incident communication.
- Restore and disaster recovery claims MUST be supported by tested recovery exercises, not merely the existence of backups.

## 9. Customer experience, accessibility, and cost control

- The core workflow MUST communicate progress, input constraints, errors, limitations, verdict provenance, and next steps clearly. Never present a fabricated result or silently convert an error into success for a smoother UI.
- Accessibility, keyboard operation, responsive behavior, comprehensible errors, and clear data-retention/privacy statements are product quality requirements, not optional polish.
- Costly features MUST have explicit resource budgets, quotas, per-tenant limits, cancellation, timeout behavior, and measured cost drivers. Do not expose an expensive endpoint without abuse controls.
- Prioritize a secure, dependable core customer journey over expanding half-tested integrations. Non-core providers and experimental surfaces must be feature-gated, clearly labeled, or disabled until their security and operating requirements are met.

## 10. LLM and external-provider boundaries

- AI-generated narrative is advisory. The deterministic/authoritative result, its evidence, and its provenance MUST remain independently available when a model is unavailable or returns malformed output.
- Treat external content, uploaded documents, model output, tool output, and provider responses as untrusted input. Validate schemas, constrain privileges, and prevent prompt injection from crossing authorization or execution boundaries.
- Send only the minimum data needed to an external provider under the approved data-handling policy. Do not transmit proprietary user structures, secrets, or customer data to a model/provider without an explicit legal/product basis and disclosed policy.
- Provider unavailability MUST produce an honest degraded state. It MUST NOT fabricate upstream evidence or silently change the scientific verdict.

## 11. Change, release, and incident discipline

Every change MUST have a clear problem statement, acceptance criteria, a risk assessment, tests, documentation updates, and a rollback/recovery plan proportionate to its impact.

Security, authorization, tenant isolation, scientific thresholds, law meanings, customer data handling, database schemas, billing/credits, and production deployment controls are high-risk surfaces. They require focused regression tests and explicit review evidence.

Incidents MUST be handled by containment first, followed by evidence preservation, root-cause analysis, remediation, and regression prevention. Never erase audit evidence or rewrite history as a substitute for credential revocation or incident response.

## 12. Definition of done

A task is DONE only when all applicable conditions below are satisfied:

1. The change addresses a verified requirement and has not silently broadened scope.
2. The exact base and proposed commit SHAs are recorded; the working tree and changed-file manifest are known.
3. Security, data ownership, scientific integrity, compatibility, and failure behavior have been reviewed for the affected paths.
4. New or changed behavior has regression, negative-path, and contract tests appropriate to its risk.
5. All required tests and CI checks pass on the proposed commit, with exact results saved as evidence.
6. Relevant docs, registry/generated artifacts, configuration examples, and operational runbooks are updated.
7. Deployment prerequisites, monitoring, and rollback/recovery are verified for production-impacting changes.
8. The final report lists changes, evidence, unresolved risks, blocked checks, and what remains unverified. "Not run" is never reported as "passed."
9. No secrets are exposed; no unapproved scientific definition/threshold changes occurred; no unverified feature is advertised as production-ready.

Any applicable P0 security vulnerability, cross-tenant data exposure, fabricated authoritative verdict, data corruption, failed mandatory gate, or undocumented scientific-semantic change blocks release.

## 13. Prohibited shortcuts

The coder MUST NOT:

- code before completing the pre-coding protocol;
- claim testing, auditing, key rotation, branch protection, or deployment verification without evidence;
- bypass, weaken, or delete a gate solely to make CI green;
- commit secrets or include them in reports, artifacts, logs, prompts, or issue text;
- permit implicit open authentication or unowned access to customer data;
- use mock/synthetic data as production evidence;
- convert errors, absent metrics, timeouts, or failed simulations into normal verdicts;
- alter law meanings, thresholds, IDs, or published claims without an explicit decision and regression evidence;
- delete or relocate modules before automated dependency/build-reference analysis;
- merge/deploy a commit that did not pass the required checks; or
- call a task complete while material tests or checks remain unrun without prominently stating that limitation.

## 14. Required task artifacts

For each meaningful change, generate a task artifact set that includes:

- preflight record and exact baseline SHA;
- scope and dependency/call-site map for affected code;
- machine-readable test/check results bound to the proposed commit;
- changed-file manifest and patch hash;
- security/data/scientific impact note, including unresolved conflicts;
- migration, deployment, and rollback plan where applicable; and
- final completed/blocked/unverified checklist.

Artifacts MUST be produced by repeatable commands where possible. They must be redacted for secrets and sensitive customer data. Preserve evidence that matters to review and reproducibility without filling the repository with redundant generated files.

## 15. External engineering references

Use these as living engineering references and map applicable requirements into repository checks. They are not claims that this project is certified or compliant:

- NIST Secure Software Development Framework (SSDF), SP 800-218: https://csrc.nist.gov/pubs/sp/800/218/final
- OWASP Application Security Verification Standard (ASVS), stable version 5.0.0 as of this document's preparation: https://owasp.org/www-project-application-security-verification-standard/ and https://github.com/OWASP/ASVS/releases
- SLSA specification v1.2 for source/build provenance and artifact verification: https://slsa.dev/spec/v1.2/.

## Final rule

**Build trust before features. Prove behavior before claiming success. Preserve evidence before changing authority. Automate every invariant that can be checked. If a change cannot be shown to be safe, correct, and reproducible, it is not ready for production.**
