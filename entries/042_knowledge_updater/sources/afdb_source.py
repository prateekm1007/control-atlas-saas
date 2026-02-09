"""
Checks for AlphaFold DB updates.
"""
import requests
import re

AFDB_API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"

def get_latest_version(uniprot_id):
    """
    Check the latest model version for a UniProt ID.
    Returns: version string (e.g. 'v4', 'v6') or None
    """
    url = f"{AFDB_API_BASE}/{uniprot_id.upper()}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return None
        
        data = response.json()
        if not data:
            return None
            
        # Parse version from download URL
        # e.g. https://.../AF-P01116-F1-model_v4.cif
        cif_url = data[0].get("cifUrl", "")
        match = re.search(r"model_(v\d+)\.cif", cif_url)
        
        if match:
            return match.group(1)
        return "unknown"
        
    except Exception:
        return None
