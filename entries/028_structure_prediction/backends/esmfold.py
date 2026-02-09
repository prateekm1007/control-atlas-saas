"""
ESMFold API Backend
"""

import requests

ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

def predict_structure(sequence: str, timeout: int = 300) -> dict:
    sequence = sequence.upper().strip()
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    
    invalid_chars = [aa for aa in sequence if aa not in valid_aa]
    if invalid_chars:
        return {
            "status": "ERROR",
            "pdb_string": None,
            "confidence_global": 0.0,
            "confidence_per_residue": [],
            "error": f"Invalid amino acid characters: {set(invalid_chars)}"
        }
    
    if len(sequence) < 10:
        return {
            "status": "ERROR",
            "pdb_string": None,
            "confidence_global": 0.0,
            "confidence_per_residue": [],
            "error": "Sequence too short (minimum 10 residues)"
        }
    
    if len(sequence) > 400:
        return {
            "status": "ERROR",
            "pdb_string": None,
            "confidence_global": 0.0,
            "confidence_per_residue": [],
            "error": "Sequence too long for API (maximum 400 residues)"
        }

    try:
        response = requests.post(
            ESMFOLD_API_URL,
            data=sequence,
            headers={"Content-Type": "text/plain"},
            timeout=timeout
        )
        
        if response.status_code != 200:
            return {
                "status": "ERROR",
                "pdb_string": None,
                "confidence_global": 0.0,
                "confidence_per_residue": [],
                "error": f"API returned status {response.status_code}"
            }
        
        pdb_string = response.text
        
        if not pdb_string.strip():
            return {
                "status": "ERROR",
                "pdb_string": None,
                "confidence_global": 0.0,
                "confidence_per_residue": [],
                "error": "API returned empty response"
            }
        
        has_pdb_content = any(
            line.startswith(("HEADER", "ATOM", "REMARK", "MODEL", "TITLE"))
            for line in pdb_string.split("\n")[:50]
        )
        
        if not has_pdb_content:
            return {
                "status": "ERROR",
                "pdb_string": None,
                "confidence_global": 0.0,
                "confidence_per_residue": [],
                "error": "API returned non-PDB content"
            }
        
        # Extract pLDDT from B-factor column
        plddt_scores = []
        for line in pdb_string.split("\n"):
            if line.startswith("ATOM") and " CA " in line:
                try:
                    bfactor = float(line[60:66].strip())
                    plddt_scores.append(bfactor)
                except (ValueError, IndexError):
                    continue
        
        if plddt_scores:
            max_score = max(plddt_scores)
            if max_score > 1.0:
                # 0-100 scale
                confidence_global = sum(plddt_scores) / len(plddt_scores) / 100.0
                confidence_per_residue = [round(s / 100.0, 3) for s in plddt_scores]
            else:
                # Already 0-1 scale
                confidence_global = sum(plddt_scores) / len(plddt_scores)
                confidence_per_residue = [round(s, 3) for s in plddt_scores]
        else:
            confidence_global = 0.5
            confidence_per_residue = []
        
        return {
            "status": "SUCCESS",
            "pdb_string": pdb_string,
            "confidence_global": round(confidence_global, 3),
            "confidence_per_residue": confidence_per_residue,
            "error": None
        }
        
    except requests.exceptions.Timeout:
        return {
            "status": "ERROR",
            "pdb_string": None,
            "confidence_global": 0.0,
            "confidence_per_residue": [],
            "error": f"API request timed out after {timeout}s"
        }
    except requests.exceptions.RequestException as e:
        return {
            "status": "ERROR",
            "pdb_string": None,
            "confidence_global": 0.0,
            "confidence_per_residue": [],
            "error": f"API request failed: {str(e)}"
        }
