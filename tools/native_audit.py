"""
Sovereign Judge: Geometric Audit for PPI Interfaces
Version: 6.1-Hardened

Methodology:
- Heavy atoms only (excludes H)
- Protein atoms only (ATOM records, excludes HETATM)
- Explicit chain ID specification (no guessing)
- Clash threshold: 2.5 Å (MDI LAW-148)
- Contact threshold: 4.5 Å (standard for polar/hydrophobic)
"""

import numpy as np
from scipy.spatial.distance import cdist
import json
import sys

try:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
except ImportError:
    print("ERROR: Biopython required. Run: pip install biopython")
    sys.exit(1)


class SovereignJudge:
    
    # MDI-compliant thresholds
    CLASH_THRESHOLD = 2.5  # Angstroms - LAW-148
    CONTACT_THRESHOLD = 4.5  # Angstroms - standard interface cutoff
    
    def __init__(self, clash_threshold=None, contact_threshold=None):
        self.clash_threshold = clash_threshold or self.CLASH_THRESHOLD
        self.contact_threshold = contact_threshold or self.CONTACT_THRESHOLD

    def load_structure(self, path):
        """Load MMCIF and extract atom data."""
        d = MMCIF2Dict(path)
        
        return {
            'coords': np.stack([
                np.array(d['_atom_site.Cartn_x'], dtype=float),
                np.array(d['_atom_site.Cartn_y'], dtype=float),
                np.array(d['_atom_site.Cartn_z'], dtype=float)
            ], axis=1),
            'chain_ids': np.array(d['_atom_site.label_asym_id']),
            'elements': np.array(d.get('_atom_site.type_symbol', ['C'] * len(d['_atom_site.Cartn_x']))),
            'record_types': np.array(d.get('_atom_site.group_PDB', ['ATOM'] * len(d['_atom_site.Cartn_x']))),
            'residue_ids': np.array(d.get('_atom_site.label_seq_id', ['0'] * len(d['_atom_site.Cartn_x']))),
            'atom_names': np.array(d.get('_atom_site.label_atom_id', ['CA'] * len(d['_atom_site.Cartn_x'])))
        }

    def get_chain_mask(self, data, chain_id):
        """Create mask for heavy protein atoms in specified chain."""
        chain_mask = data['chain_ids'] == chain_id
        heavy_mask = data['elements'] != 'H'
        protein_mask = data['record_types'] == 'ATOM'
        
        combined = chain_mask & heavy_mask & protein_mask
        
        if not np.any(combined):
            available = sorted(set(data['chain_ids']))
            raise ValueError(f"No atoms found for chain '{chain_id}'. Available: {available}")
        
        return combined

    def calculate_contact_density(self, path, target_chain, binder_chain):
        """
        Calculate interface metrics between two chains.
        
        Returns:
            dict with rho, min_distance, clashes, atom counts, status
        """
        data = self.load_structure(path)
        
        target_mask = self.get_chain_mask(data, target_chain)
        binder_mask = self.get_chain_mask(data, binder_chain)
        
        tar_coords = data['coords'][target_mask]
        bin_coords = data['coords'][binder_mask]
        
        # Distance matrix
        dists = cdist(tar_coords, bin_coords)
        
        min_dist = float(np.min(dists))
        contacts = int(np.sum(dists < self.contact_threshold))
        clashes = int(np.sum(dists < self.clash_threshold))
        
        # Interface residues (for verification)
        contact_mask = dists < self.contact_threshold
        tar_residues = data['residue_ids'][target_mask]
        bin_residues = data['residue_ids'][binder_mask]
        
        interface_target = sorted(set(tar_residues[contact_mask.any(axis=1)]))
        interface_binder = sorted(set(bin_residues[contact_mask.any(axis=0)]))
        
        # Status determination
        if clashes > 0:
            status = "CLASH_VETO"
        elif contacts < 50:
            status = "LOW_CONTACT_WARNING"
        else:
            status = "PASS"
        
        return {
            'rho': contacts,
            'min_distance_A': round(min_dist, 2),
            'clashes': clashes,
            'target_atoms': int(np.sum(target_mask)),
            'binder_atoms': int(np.sum(binder_mask)),
            'interface_target_residues': interface_target[:10],  # First 10 for brevity
            'interface_binder_residues': interface_binder[:10],
            'status': status,
            'thresholds': {
                'clash': self.clash_threshold,
                'contact': self.contact_threshold
            }
        }

    def list_chains(self, path):
        """List available chains and their atom counts."""
        data = self.load_structure(path)
        chains = sorted(set(data['chain_ids']))
        
        result = {}
        for c in chains:
            mask = self.get_chain_mask(data, c) if c in data['chain_ids'] else np.zeros(len(data['chain_ids']), dtype=bool)
            try:
                mask = (data['chain_ids'] == c) & (data['elements'] != 'H') & (data['record_types'] == 'ATOM')
                result[c] = int(np.sum(mask))
            except:
                result[c] = 0
        
        return result


def main():
    """Command-line interface for Sovereign Judge."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Sovereign Judge: Geometric Audit")
    parser.add_argument("cif_path", help="Path to MMCIF structure file")
    parser.add_argument("--target", "-t", help="Target chain ID")
    parser.add_argument("--binder", "-b", help="Binder chain ID")
    parser.add_argument("--list-chains", "-l", action="store_true", help="List available chains")
    
    args = parser.parse_args()
    
    judge = SovereignJudge()
    
    if args.list_chains:
        chains = judge.list_chains(args.cif_path)
        print("Available chains (heavy protein atoms):")
        for chain, count in sorted(chains.items(), key=lambda x: -x[1]):
            print(f"  {chain}: {count} atoms")
        return
    
    if not args.target or not args.binder:
        print("ERROR: Must specify --target and --binder chain IDs")
        print("Use --list-chains to see available chains")
        sys.exit(1)
    
    result = judge.calculate_contact_density(args.cif_path, args.target, args.binder)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
