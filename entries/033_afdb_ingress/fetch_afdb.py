#!/usr/bin/env python3
"""
AlphaFold DB fetcher (2025 reliable method).
Returns BOTH CIF (for pLDDT) and PDB (for fpocket).
"""
import requests
import gzip
import shutil
from pathlib import Path

API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"

def cif_to_pdb(cif_path: Path) -> Path:
    """Convert CIF to PDB using gemmi."""
    pdb_path = cif_path.with_suffix(".pdb")
    if pdb_path.exists() and pdb_path.stat().st_size > 10000:
        return pdb_path
    
    try:
        import gemmi
        structure = gemmi.read_structure(str(cif_path))
        structure.write_pdb(str(pdb_path))
        print(f"    Converted to PDB: {pdb_path.name}")
        return pdb_path
    except Exception as e:
        print(f"[!] CIF→PDB conversion failed: {e}")
        return None

def fetch_structure(uniprot_id: str, output_dir: Path) -> dict:
    """
    Fetch structure and return BOTH paths.
    Returns: {"cif": path, "pdb": path} or None
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    api_url = f"{API_BASE}/{uniprot_id.upper()}"
    try:
        response = requests.get(api_url, timeout=60)
        if response.status_code != 200:
            print(f"[!] API failed for {uniprot_id}: HTTP {response.status_code}")
            return None
        
        data = response.json()
        if not data:
            print(f"[!] No models for {uniprot_id}")
            return None
        
        model = data[0]
        cif_url = model.get("cifUrl")
        
        if not cif_url:
            print(f"[!] No CIF URL for {uniprot_id}")
            return None
        
        cif_path = output_dir / Path(cif_url).name
        
        # Check cache
        pdb_path = cif_path.with_suffix(".pdb")
        if cif_path.exists() and pdb_path.exists() and pdb_path.stat().st_size > 10000:
            return {"cif": str(cif_path), "pdb": str(pdb_path)}
        
        # Download CIF
        cif_response = requests.get(cif_url, timeout=120)
        if cif_response.status_code != 200:
            print(f"[!] CIF download failed for {uniprot_id}")
            return None
            
        with open(cif_path, "wb") as f:
            f.write(cif_response.content)
        print(f"[+] Downloaded CIF: {cif_path.name}")
        
        # Convert to PDB
        pdb_path = cif_to_pdb(cif_path)
        if not pdb_path:
            return None
        
        return {"cif": str(cif_path), "pdb": str(pdb_path)}
        
    except Exception as e:
        print(f"[!] Error for {uniprot_id}: {e}")
        return None
