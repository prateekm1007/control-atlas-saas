# Entry 027 — General Pocket Detection

Automatic druggable pocket detection for any PDB structure.

## Requirements

Install fpocket: conda install -c conda-forge fpocket

## Usage

Detect pockets: python pocket_detector.py --pdb structure.pdb --json

With confidence: python pocket_detector.py --pdb structure.pdb --confidence 0.93

## Output

Each pocket is classified as:
- VALIDATED: High confidence, passes all physics gates
- CANDIDATE: Moderate confidence, may need manual review  
- REJECTED: Fails physics thresholds (with reasons)

## Integration

Output feeds directly into Universal Gatekeeper (Entry 020) for compound screening.
