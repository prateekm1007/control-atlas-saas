# Control Atlas: A Physics-First Platform for Allosteric Control Discovery

## Executive Summary
Control Atlas is a deterministic, law-enforcing platform for discovering and validating chemical control sites on difficult protein targets (e.g., KRAS). Unlike traditional docking or black-box AI, it enforces **physics-based constraints** before candidate generation.

## Core Thesis
> "Geometry finds candidates; physics decides control."

## Methodology

### Phase 1: Mechanism Mapping
- Deterministic analysis of target dynamics (KRAS G12C/WT).
- Identification of "Switch" motifs as control surfaces.

### Phase 2: Control Discovery (The Falsification)
- **Geometry is Insufficient**: Pure backbone RMSD search yielded false positives (Ubiquitin).
- **Physics is Required**: Only concave, addressable pockets with specific polarity gradients enable control.
- **The Manifold**: Valid control sites occupy a constrained deformation manifold (closed → open).

### Phase 3: The Gatekeeper (Productization)
- **Chemistry Grammar**: Rules derived from pocket physics (Entry 018).
- **Resistance Awareness**: Pre-emptive rejection of scaffolds vulnerable to known escape mutations (Entry 017).
- **Deterministic Validation**: A strict PASS/REJECT compiler for candidates (Entry 020/021).

## Key Results
- **Sotorasib Recovery**: The platform deterministically validates the clinical drug without special casing.
- **Decoy Rejection**: >90% of random library compounds are rejected by physics/grammar rules.
- **Resistance Prediction**: Correctly flags R68S, Y96D, and Q99L as resistance risks.

## What Control Atlas Is Not
- **Not a docking-first screening platform**
- **Not a black-box machine learning model**
- **Not a random molecule generator**

Control Atlas is a **law-enforcing compiler** that constrains chemical search *before* optimization or simulation.

## Conclusion
Control Atlas converts protein structure into a **computable search space**, enabling high-throughput screening with physics-grade fidelity.
