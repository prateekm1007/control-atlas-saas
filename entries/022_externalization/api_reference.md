# Control Atlas API Reference

## Candidate Validation (Gatekeeper)

### Single Molecule
```bash
python3 entries/020_candidate_validation/validate_candidate.py --smiles "<SMILES>"
```

### Batch Filtering
```bash
python3 entries/021_library_filtering/filter_library.py --input <FILE.smi> --output <DIR>
```

## Scaffold Generation
```bash
python3 entries/019_generative_control/define_scaffold_space.py
```

## Resistance Simulation
```bash
python3 entries/017_failure_simulation/simulate_resistance.py
```
