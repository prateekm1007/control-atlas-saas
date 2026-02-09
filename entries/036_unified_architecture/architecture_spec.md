# Unified Theory of Computational Falsification
**Architecture Specification v1.0**

## Core Philosophy
We do not score molecules. We apply successive constraints to **falsify** invalid paths. The remaining space is "The Clearing" — the high-probability route.

## The Stack

### Layer 1: Physics (Spacetime)
- **Input:** Structure (CIF), Pockets (fpocket)
- **Action:** Geometric Hard Veto
- **Logic:** Volume, Enclosure, pLDDT Stability
- **Output:** `PocketExists` (True/False)

### Layer 2: Chemistry (Matter)
- **Input:** Molecule (SMILES), Pocket Constraints
- **Action:** Energetic Viability Check
- **Logic:** Polar/Hydrophobic Mismatch, Steric Clash, Strain
- **Output:** `ChemicallyTractable` (True/False)

### Layer 3: Biology (Context)
- **Input:** Target ID, Disease Context
- **Action:** Relevance Filter
- **Logic:** Expression (GTEx), Mutation (ClinVar), Essentiality (DepMap)
- **Output:** `BiologicallyRelevant` (True/False)

### Layer 4: Math (Confidence)
- **Input:** Constraint Margins from Layers 1-3
- **Action:** Robustness Calculation
- **Logic:** Bayesian integration of uncertainty
- **Output:** `NavigationPath` (Confidence + Advice)

## Execution Model
The system runs sequentially. A failure at Layer N stops execution immediately (Fail Fast).

1. `Physics.check()` -> Fail? -> **STOP** (Reason: "No Pocket")
2. `Chemistry.check()` -> Fail? -> **STOP** (Reason: "Untractable")
3. `Biology.check()` -> Fail? -> **STOP** (Reason: "Irrelevant")
4. `Math.integrate()` -> **Success** (Output: "Pursue this path")
