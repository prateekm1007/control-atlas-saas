# Entry 029 — End-to-End Sequence Screening

Single command: Sequence to ranked compound hits.

## Usage

python screen_sequence.py --sequence "MKTVRQ..." --smi library.smi --out results.csv

## Pipeline

1. Structure Prediction (ESMFold) - Entry 028
2. Pocket Detection (fpocket + physics) - Entry 027  
3. Compound Screening (Gatekeeper) - Entry 020

## Confidence Propagation

Final confidence = structure_confidence x pocket_confidence x compound_confidence

All rejections include reasons. No silent failures.
