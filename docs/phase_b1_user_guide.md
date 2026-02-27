# Toscanini Phase B1 — User Guide
## Managed Refinement Orchestration (Beta)

### Overview

Toscanini B1 closes the structural governance loop:

Upload PDB → Audit → Diagnose → Prescribe → Refine (local) → Callback → Compare

You run compute locally (Rosetta/OpenMM on your hardware or Kaggle/Colab).
Toscanini orchestrates the before/after comparison automatically.

### Token Details

| Property | Value |
|----------|-------|
| Expiry | 7 days from generation |
| Uses | Single-use (consumed after upload) |
| Max callbacks per audit | 3 |
| Format | JWT (HS256) |

### Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| Token already used | Token was consumed | Generate new token |
| Token expired | >7 days since generation | Re-generate |
| Maximum callbacks reached | 3 callbacks used | Re-audit original |
| Invalid token | Copy/paste error | Check for spaces |
