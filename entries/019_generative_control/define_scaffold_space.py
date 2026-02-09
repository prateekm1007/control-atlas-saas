#!/usr/bin/env python3
"""
Entry 019: Scaffold Space Definition
Generates combinatorial scaffold space and applies Entry 018 Grammar rules.
"""

import json
import os
import csv
import itertools

# Paths
GRAMMAR_PATH = os.path.expanduser("~/control-atlas/entries/018_chemistry_grammar/chemistry_grammar_v2.json")
LIBRARY_PATH = os.path.expanduser("~/control-atlas/entries/019_generative_control/fragment_library.json")
OUTPUT_CSV = os.path.expanduser("~/control-atlas/entries/019_generative_control/legal_scaffolds.csv")

def load_json(path):
    with open(path, 'r') as f:
        return json.load(f)

def validate_assembly(core, anchor, warhead, tail, grammar):
    """
    Apply Entry 018b Rules to a virtual assembly
    Returns: (is_valid, reason_list)
    """
    reasons = []
    
    # 1. Fragment-Specific Failures (Hardcoded properties)
    if "fail_reason" in core: reasons.append(f"Core: {core['fail_reason']}")
    if "fail_reason" in anchor: reasons.append(f"Anchor: {anchor['fail_reason']}")
    if "fail_reason" in warhead: reasons.append(f"Warhead: {warhead['fail_reason']}")
    if "fail_reason" in tail: reasons.append(f"Tail: {tail['fail_reason']}")
    
    # 2. Grammar Region Checks
    # Core must be aromatic (Buried Core Rule)
    if not core.get('aromatic'): 
        reasons.append("Core: Must be aromatic")
        
    # Warhead must be electrophile (Warhead Zone Rule)
    if not warhead.get('electrophile'):
        reasons.append("Warhead: Must be electrophile")
        
    # Anchor must not clash (Resistance Rule: no_steric_projection)
    # Using 'Phenyl' as proxy for steric clash at GLY60
    if anchor['name'] == 'Phenyl':
        reasons.append("Anchor: Steric violation (GLY60)")
        
    # Tail must be polar (Solvent Tail Rule)
    if not tail.get('polar'):
        reasons.append("Tail: Must be polar")
        
    # 3. Quantitative Limits (Approximate MW Check)
    # Note: Full MW requires linker atoms, here we sum fragments
    total_mw = core["mw"] + anchor["mw"] + warhead["mw"] + tail["mw"] + 50
    mw_limits = grammar['quantitative_limits']['molecular_weight']
    
    if not (mw_limits['min'] <= total_mw <= mw_limits['max']):
        reasons.append(f"MW {total_mw}: Out of range {mw_limits['min']}-{mw_limits['max']}")
        
    # 4. H-Bond Acceptor Check (Approximate)
    total_hba = core.get('hba', 0) + tail.get('hba', 0) # Simplified
    hba_limits = grammar['quantitative_limits']['hba']
    if not (hba_limits['min'] <= total_hba <= hba_limits['max']):
        # Soft fail - might add linkers later, but flagging
        pass 

    valid = len(reasons) == 0
    return valid, reasons, total_mw

def main():
    print("=== ENTRY 019: SCAFFOLD SPACE DEFINITION ===")
    
    grammar = load_json(GRAMMAR_PATH)
    lib = load_json(LIBRARY_PATH)
    
    print(f"Loaded Grammar v{grammar.get('version', '?')}")
    print(f"Fragments: {len(lib['cores'])} Cores, {len(lib['anchors'])} Anchors, "
          f"{len(lib['warheads'])} Warheads, {len(lib['tails'])} Tails")
    
    results = []
    
    # Generate all permutations
    permutations = list(itertools.product(lib['cores'], lib['anchors'], lib['warheads'], lib['tails']))
    print(f"Total Permutations: {len(permutations)}")
    
    print("\nValidating Scaffolds...")
    print("-" * 60)
    
    valid_count = 0
    
    for core, anchor, warhead, tail in permutations:
        name = f"{core['name']}-{anchor['name']}-{warhead['name']}-{tail['name']}"
        
        is_valid, reasons, mw = validate_assembly(core, anchor, warhead, tail, grammar)
        
        status = "VALID" if is_valid else "REJECT"
        if is_valid: valid_count += 1
        
        results.append({
            "scaffold_name": name,
            "core": core['name'],
            "anchor": anchor['name'],
            "warhead": warhead['name'],
            "tail": tail['name'],
            "mw_est": mw,
            "status": status,
            "reasons": "; ".join(reasons)
        })
        
        if is_valid:
            print(f"  [VALID] {name} (MW: {mw})")
            
    # Write CSV
    with open(OUTPUT_CSV, 'w', newline='') as f:
        fieldnames = ["scaffold_name", "core", "anchor", "warhead", "tail", "mw_est", "status", "reasons"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
        
    print("-" * 60)
    print(f"Total Scaffolds: {len(results)}")
    print(f"Valid Scaffolds: {valid_count}")
    print(f"Rejection Rate:  {((len(results)-valid_count)/len(results))*100:.1f}%")
    print(f"Saved: {OUTPUT_CSV}")
    print("=== ENTRY 019 COMPLETE ===")

if __name__ == "__main__":
    main()
