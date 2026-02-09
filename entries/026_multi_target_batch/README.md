# Entry 026 — Multi-Target Batch Screening

Deterministic multi-target orchestration layer.

## Guarantees
- Universal Gatekeeper is sole decision authority
- One initialization per target
- No duplicated physics or grammar
- Long-form CSV output (compound × target)

## Usage

Screen all validated targets:
```bash
python multi_target_screen.py \
  --smi library.smi \
  --validated-only \
  --out results.csv
# Close the previously opened heredoc cleanly
