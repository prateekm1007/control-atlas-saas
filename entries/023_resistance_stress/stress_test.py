#!/usr/bin/env python3
import sys
import os
import csv
import json

sys.path.insert(0, os.path.expanduser('~/control-atlas/entries/020_candidate_validation'))
from validate_candidate import validate

MUTATIONS = {
    'WT': {},
    'R68S': {'rule': 'none'},
    'Y96D': {'rule': 'clogp_max', 'threshold': 3.5, 'reason': 'Hydrophobic anchor incompatible with Asp96'},
    'Q99L': {'rule': 'tpsa_max', 'threshold': 100, 'reason': 'Polar scaffold destabilized by Leu99'},
    'G60W': {'rule': 'mw_max', 'threshold': 500, 'reason': 'Steric clash at Trp60'}
}

SCAFFOLDS = {
    'Sotorasib': 'CC(C)N1CCN(CC1)C2=NC=NC3=C2C(=CN3C4=C(C=C(C=C4F)F)C5=C(C=CC=N5)F)N6CCC(CC6)C(=O)C=C',
    'Novel_Scaffold_A': 'CC(C)N1CCN(CC1)C2=NC=NC3=C2C(=CN3C4=C(C=C(C=C4F)F)C5=C(C=CC=N5)F)N6CCC(CC6)C(=O)C=C',
    'Small_Frag_B': 'C=CC(=O)N1CCC(CC1)C2=NC=NC3=C2C(=CN3C4=CCCC4)N5CCCCC5'
}

def stress_test():
    print('=== RESISTANCE STRESS CAMPAIGN (EVIDENCE-DRIVEN) ===')
    results = []

    for scaf_name, smiles in SCAFFOLDS.items():
        print(f'\nTesting Scaffold: {scaf_name}')
        baseline = validate(smiles)

        if baseline['status'] != 'VALID':
            print(f'  [SKIP] Baseline Invalid')
            continue

        metrics = baseline.get('metrics', {})
        print(f"  Metrics: MW={metrics.get('mw')}, cLogP={metrics.get('clogp')}, TPSA={metrics.get('tpsa')}")
        
        for mut_name, rules in MUTATIONS.items():
            if mut_name == 'WT':
                continue

            outcome = 'ROBUST'
            reason = '-'
            rule = rules.get('rule')
            limit = rules.get('threshold', 0)

            if rule == 'mw_max' and metrics.get('mw', 0) > limit:
                outcome = 'VULNERABLE'
                reason = rules['reason']
            elif rule == 'clogp_max' and metrics.get('clogp', 0) > limit:
                outcome = 'VULNERABLE'
                reason = rules['reason']
            elif rule == 'tpsa_max' and metrics.get('tpsa', 0) > limit:
                outcome = 'VULNERABLE'
                reason = rules['reason']

            print(f'  vs {mut_name}: {outcome} ({reason})')
            results.append({
                'scaffold': scaf_name,
                'mutation': mut_name,
                'outcome': outcome,
                'reason': reason
            })

    out_path = os.path.expanduser('~/control-atlas/entries/023_resistance_stress/resistance_map.csv')
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['scaffold', 'mutation', 'outcome', 'reason'])
        w.writeheader()
        w.writerows(results)

    print(f'\nSaved map to {out_path}')

if __name__ == '__main__':
    stress_test()
