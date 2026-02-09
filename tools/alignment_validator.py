"""
Alignment Validator: Local Motif RMSD Calculation
Version: 1.0

Methodology:
- Align on target protein backbone (chain A/target)
- Calculate RMSD on specified motif residues only
- Backbone atoms: N, CA, C, O
"""

import numpy as np
import sys

try:
    from Bio.PDB import PDBParser, MMCIFParser, Superimposer
except ImportError:
    print("ERROR: Biopython required. Run: pip install biopython")
    sys.exit(1)


def load_structure(path):
    """Load PDB or MMCIF structure."""
    if path.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("structure", path)[0]


def get_backbone_atoms(chain, residue_ids):
    """
    Extract backbone atoms for specified residues.
    
    Parameters:
        chain: Bio.PDB Chain object
        residue_ids: list of residue sequence numbers (int)
    
    Returns:
        list of Atom objects
    """
    atoms = []
    backbone_names = ['N', 'CA', 'C', 'O']
    
    for res_id in residue_ids:
        # Try different residue ID formats
        for insert_code in [' ', '', 'A']:
            try:
                res = chain[(' ', res_id, insert_code)]
                for atom_name in backbone_names:
                    if atom_name in res:
                        atoms.append(res[atom_name])
                break
            except KeyError:
                continue
    
    return atoms


def validate_warhead_rmsd(
    predicted_path, 
    reference_path,
    pred_chain,
    ref_chain,
    motif_residues_pred,
    motif_residues_ref
):
    """
    Calculate local RMSD for a motif (e.g., YWPTG warhead).
    
    Parameters:
        predicted_path: path to predicted structure
        reference_path: path to reference structure (e.g., 5NIU)
        pred_chain: chain ID in predicted structure
        ref_chain: chain ID in reference structure
        motif_residues_pred: list of residue IDs in predicted
        motif_residues_ref: list of residue IDs in reference
    
    Returns:
        dict with RMSD, atom counts, methodology
    """
    pred_struct = load_structure(predicted_path)
    ref_struct = load_structure(reference_path)
    
    try:
        pred_atoms = get_backbone_atoms(pred_struct[pred_chain], motif_residues_pred)
        ref_atoms = get_backbone_atoms(ref_struct[ref_chain], motif_residues_ref)
    except KeyError as e:
        raise ValueError(f"Chain not found: {e}")
    
    if len(pred_atoms) == 0:
        raise ValueError(f"No backbone atoms found in predicted structure for residues {motif_residues_pred}")
    if len(ref_atoms) == 0:
        raise ValueError(f"No backbone atoms found in reference structure for residues {motif_residues_ref}")
    if len(pred_atoms) != len(ref_atoms):
        raise ValueError(f"Atom count mismatch: predicted={len(pred_atoms)}, reference={len(ref_atoms)}")
    
    sup = Superimposer()
    sup.set_atoms(ref_atoms, pred_atoms)
    
    return {
        'rmsd_A': round(sup.rms, 4),
        'atom_count': len(pred_atoms),
        'residue_count': len(motif_residues_pred),
        'methodology': 'Backbone atoms (N, CA, C, O) of specified motif residues',
        'predicted_residues': motif_residues_pred,
        'reference_residues': motif_residues_ref
    }


def main():
    import argparse
    import json
    
    parser = argparse.ArgumentParser(description="Warhead RMSD Validator")
    parser.add_argument("predicted", help="Predicted structure path")
    parser.add_argument("reference", help="Reference structure path (e.g., 5NIU.pdb)")
    parser.add_argument("--pred-chain", required=True, help="Chain ID in predicted")
    parser.add_argument("--ref-chain", required=True, help="Chain ID in reference")
    parser.add_argument("--pred-residues", required=True, help="Comma-separated residue IDs in predicted")
    parser.add_argument("--ref-residues", required=True, help="Comma-separated residue IDs in reference")
    
    args = parser.parse_args()
    
    pred_res = [int(x.strip()) for x in args.pred_residues.split(',')]
    ref_res = [int(x.strip()) for x in args.ref_residues.split(',')]
    
    result = validate_warhead_rmsd(
        args.predicted,
        args.reference,
        args.pred_chain,
        args.ref_chain,
        pred_res,
        ref_res
    )
    
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
