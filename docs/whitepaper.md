# Control Atlas
## A Unified Framework for Constraint-Based Navigation in Drug Discovery

### Abstract
Drug discovery pipelines fail disproportionately due to prolonged exploration of chemically and biologically infeasible hypotheses rather than an inability to generate candidate molecules. We present **Control Atlas**, a constraint-based, physics-first computational framework that formalizes *falsification* as the primary driver of search-space reduction. Instead of ranking candidates by predicted activity, Control Atlas incrementally collapses infeasible regions of protein–ligand space using deterministic physical, chemical, and biological constraints, yielding a mathematically defined **navigation signal** representing robustness rather than probability of success.

---

## 1. Conceptual Motivation
Modern computational drug discovery systems predominantly operate in a **generative or scoring paradigm**:
* docking scores
* ML-predicted affinities
* end-to-end black-box predictions

These systems optimize *likelihood of activity* but rarely quantify **distance to failure**.

In contrast, experimental drug discovery is dominated by **constraint violations**:
* absence of a stable binding pocket
* thermodynamic incompatibility
* lack of biological relevance

Control Atlas reframes drug discovery as a **constraint satisfaction and navigation problem**, not a prediction problem.

---

## 2. Core Principle: Falsification as Navigation
We adopt the following axiom:

> *The probability of drug discovery success increases primarily through early and rigorous elimination of infeasible hypotheses.*

Rather than asking *"Which molecule binds best?"*, we ask:
> **"Which regions of chemical–target space are not yet ruled out by first principles?"**

The surviving region defines the **path of least resistance**.

---

## 3. Layered Constraint Architecture
Control Atlas decomposes feasibility into four orthogonal but composable layers.

### 3.1 Physics Layer — Geometric Viability
**Purpose:** Determine whether a stable binding geometry exists.
**Inputs:** AlphaFold CIF structures, per-residue confidence (pLDDT), pocket geometry (fpocket-derived).
**Constraint Examples:** Pocket volume bounds, enclosure metrics, structural confidence thresholds.
**Property:** Deterministic, binary veto, topological collapse on failure. If no physically stable pocket exists, downstream chemistry and biology are undefined.

### 3.2 Chemistry Layer — Energetic and Tractability Constraints
**Purpose:** Determine whether matter can plausibly occupy the surviving geometric space.
**Inputs:** Molecular descriptors (LogP, MW, polarity), ensemble conformers (ETKDGv3, MMFF), pocket physicochemical features.
**Key Design Choice:** Chemistry does **not** predict binding; it falsifies *energetically hostile* hypotheses. Bias-laden metrics (e.g., QED) are treated as **warnings**, not vetoes.

### 3.3 Biology Layer — Functional Relevance
**Purpose:** Determine whether binding is causally meaningful in a disease context.
**Inputs:** Essentiality (DepMap proxies), expression context (GTEx/HPA), driver/passenger classification.
**Interpretation:** A physically and chemically valid interaction is discarded if it lacks biological consequence.

### 3.4 Math Layer — Robustness Integration
**Purpose:** Quantify **distance to constraint violation**, not predicted success.
Each layer emits standardized **Constraint objects**:
* `margin > 0` → feasible
* `margin ≈ 0` → fragile
* `margin < 0` → collapse

Robustness is computed as a function of surviving margins across layers. Importantly, failures collapse space; math never resurrects invalid paths.

---

## 4. Constraint Unification: The Critical Advance
All layers emit a common formal object:
```python
Constraint(layer, parameter, threshold, actual, margin, status, collapses_space, provenance)
