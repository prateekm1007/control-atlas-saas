# The Viability Calculus: A Formal Unification of Drug Discovery Constraints

**Objective:**
To unify Physics, Chemistry, Biology, and Mathematics into a single invariant expression of hypothesis viability.

---

## 1. The Fundamental Identity

We define the viability $\mathcal{V}$ of a drug discovery hypothesis $h$ as:

$$ \mathcal{V}(h) = \left[ \prod_{k \in \{P, C, B\}} \Theta(\text{margin}_k(h)) \right] \times \mathcal{R}(h) $$

Where:
*   $P, C, B$ represent the domains of Physics, Chemistry, and Biology.
*   $\text{margin}_k(h)$ is the normalized distance of hypothesis $h$ from the failure boundary in domain $k$.
*   $\Theta(x)$ is the Heaviside step function (Topology):
    *   $\Theta(x) = 0$ if $x < 0$ (Collapse/Rejection)
    *   $\Theta(x) = 1$ if $x \ge 0$ (Survival)
*   $\mathcal{R}(h)$ is the Robustness metric (Geometry), representing the integrated safety margin of the surviving path.

---

## 2. Domain Mapping

| Symbol | Domain | Operator | Meaning |
| :--- | :--- | :--- | :--- |
| $\Theta_P$ | **Physics** | **Space** | Does a stable pocket exist? (0 = No Space) |
| $\Theta_C$ | **Chemistry** | **Energy** | Is binding thermodynamically allowed? (0 = Forbidden) |
| $\Theta_B$ | **Biology** | **Causality** | Does binding cause relevant effect? (0 = Irrelevant) |
| $\mathcal{R}$ | **Math** | **Certainty** | How far are we from any boundary? |

---

## 3. Implications

### 3.1 The Collapse Principle
If any single domain fails ($\text{margin}_k < 0$), the entire viability $\mathcal{V}(h)$ collapses to zero.
$$ \exists k : \Theta_k = 0 \implies \mathcal{V}(h) = 0 $$
This formalizes the **Hard Veto**. There is no compensation. A biological zero cannot be saved by a perfect chemical score.

### 3.2 The Robustness Gradient
For viable hypotheses ($\mathcal{V}(h) > 0$), the value of $\mathcal{V}(h)$ is proportional to the robustness.
$$ \mathcal{V}(h) \propto \mathcal{R}(h) $$
We navigate by maximizing $\mathcal{R}$ within the feasible region defined by $\prod \Theta$.

---

## 4. Conclusion
This identity reveals that drug discovery is not a summation of features, but an **intersection of constraints**. The goal of the system is to define the boundaries of the intersection (The Clearing) and navigate to its center (Maximum Robustness).
