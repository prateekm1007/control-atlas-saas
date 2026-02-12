# Entry 105 — KRAS G12D × Chai-1 Terminal Veto

## Classification
- Type: TERMINAL VETO
- Scope: Generator-Specific
- Tier(s): 1

## Context
- Trigger Event: Repeated sub-2.3 Å steric overlap across v6–v8
- Prior Entries Referenced: Entry 098, Entry 101
- Objective: Determine whether failure was motif-specific or manifold-level

## Test Space
- Generator: Chai-1
- Target: KRAS G12D (Switch II)
- Variations Explored:
  - Motif length: 4–8 aa
  - Buffer length: GG → GGGG
  - Side-chain swap: Arg → Lys
  - Tail deletion: DVP removed

## Observations
- Heavy-atom distances: 1.7–2.2 Å
- Safety threshold: 2.5 Å
- Failure invariant across all variants

## Verdict
- Status: TERMINAL VETO
- Confidence: High

## Extracted Law
1. Chai-1 hallucinates a penetrable Switch II pocket in KRAS G12D
2. Increasing linker entropy does not resolve the clash
3. This is a generator-specific manifold collapse

## Governance Action
- Blacklist: KRAS_G12D × Chai-1
- Enforcement Rule: Tier 0 veto on generator-target pair

## Downstream Implications
- Saved Compute: ~400 GPU-hours
- Affected Pipelines: All KRAS programs using Chai-1
- Next Allowed Actions: Change generator or target only

## Status
- Ledger State: LOCKED
- Eligible for Revisit: No
