# Usage Examples

## Geometric Audit
python tools/native_audit.py structures/champ005/structure.cif --list-chains
python tools/native_audit.py structures/champ005/structure.cif --target A --binder B

## Energy Minimization
python tools/kinetic_rescue.py structures/champ005/structure.cif output.pdb

## View Results
cat structures/validation_summary.json | python -m json.tool
