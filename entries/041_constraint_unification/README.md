# Entry 041 — Constraint Unification

All layers emit `Constraint` objects with numeric margins.

## Why This Matters
- Math layer measures robustness, not probabilities
- Failures collapse space deterministically
- Enables ensemble and uncertainty reasoning
- Makes navigation auditable

## Rule
No layer returns booleans or bespoke result types.
Everything is a Constraint.
