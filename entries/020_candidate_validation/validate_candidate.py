#!/usr/bin/env python3
import sys
import json
import os
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

GRAMMAR_PATH = os.path.expanduser('~/control-atlas/entries/018_chemistry_grammar/chemistry_grammar_v2.json')

def load_grammar():
    try:
        with open(GRAMMAR_PATH, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f'ERROR: Grammar file not found at {GRAMMAR_PATH}')
        sys.exit(1)

def check_quantitative_limits(mol, grammar):
    limits = grammar['quantitative_limits']
    reasons = []
    
    mw = Descriptors.MolWt(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rot = Descriptors.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    clogp = Descriptors.MolLogP(mol)

    if not (limits['molecular_weight']['min'] <= mw <= limits['molecular_weight']['max']):
        reasons.append(f'MW {mw:.1f} out of range')
    if not (limits['hbd']['min'] <= hbd <= limits['hbd']['max']):
        reasons.append(f'HBD {hbd} out of range')
    if not (limits['hba']['min'] <= hba <= limits['hba']['max']):
        reasons.append(f'HBA {hba} out of range')
    
    return len(reasons) == 0, reasons, {'mw': mw, 'hbd': hbd, 'hba': hba}

def check_warhead(mol):
    WARHEADS = {
        'acrylamide': 'C(=O)C=[C,c]',
        'vinyl_sulfonamide': 'C=CS(=O)(=O)N',
        'chloroacetamide': 'ClCC(=O)N'
    }
    for name, smarts in WARHEADS.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            return True, []
    return False, ['Missing electrophile warhead']

def validate(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return {'status': 'ERROR'}
    
    grammar = load_grammar()
    report = {'status': 'VALID', 'notes': []}
    
    # 1. Quant Limits
    q_pass, q_notes, metrics = check_quantitative_limits(mol, grammar)
    if not q_pass:
        report['status'] = 'REJECT'
        report['notes'].extend(q_notes)
    report['metrics'] = metrics

    # 2. Warhead
    w_pass, w_notes = check_warhead(mol)
    if not w_pass:
        report['status'] = 'REJECT'
        report['notes'].extend(w_notes)
        
    return report

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles', required=True)
    args = parser.parse_args()
    
    result = validate(args.smiles)
    print(f'VERDICT: {result["status"]}')
    if result['status'] == 'REJECT':
        for n in result['notes']: print(f'  X {n}')
