"""
Physics Metrics for Pocket Validation
Computes volume, exposure, hydrophobicity using YOUR locked physics.
"""

# Hydrophobic residues (standard classification)
HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "TYR"}

# Quality thresholds (from locked physics engine)
THRESHOLDS = {
    "volume_min": 150.0,      # Å³ - too small = undruggable
    "volume_max": 1500.0,     # Å³ - too large = non-specific
    "hydrophobic_min": 0.30,  # At least 30% hydrophobic for small molecule binding
    "hydrophobic_max": 0.90,  # Too greasy = solubility issues
    "residue_min": 8,         # Minimum residues for real pocket
}


def compute_physics(pocket_data: dict, structure_confidence: float = 1.0) -> dict:
    """
    Compute physics metrics for a pocket.
    
    Args:
        pocket_data: From fpocket_runner
        structure_confidence: pLDDT from structure prediction
    
    Returns:
        {
            "volume": float,
            "exposure": float,
            "hydrophobic_pct": float,
            "residue_count": int,
            "confidence": float,
            "status": "VALIDATED" | "CANDIDATE" | "REJECTED",
            "rejection_reasons": list
        }
    """
    residues = pocket_data.get("residues", [])
    volume = pocket_data.get("volume", 0.0)
    druggability = pocket_data.get("druggability", 0.0)
    
    # Count residues
    residue_count = len(residues)
    
    # Compute hydrophobicity
    hydrophobic_count = 0
    for res in residues:
        # Format: "A:ALA123" or "ALA123"
        if ":" in res:
            res_name = res.split(":")[1][:3]
        else:
            res_name = res[:3]
        
        if res_name.upper() in HYDROPHOBIC_RESIDUES:
            hydrophobic_count += 1
    
    hydrophobic_pct = hydrophobic_count / residue_count if residue_count > 0 else 0.0
    
    # Exposure proxy (inverse of druggability score - higher druggability = more buried)
    exposure = 1.0 - min(druggability, 1.0)
    
    # Compute confidence (structure confidence * pocket quality)
    pocket_quality = druggability  # fpocket druggability as quality proxy
    confidence = structure_confidence * pocket_quality
    
    # Apply quality gates
    rejection_reasons = []
    
    if residue_count < THRESHOLDS["residue_min"]:
        rejection_reasons.append(f"too_few_residues ({residue_count} < {THRESHOLDS['residue_min']})")
    
    if volume < THRESHOLDS["volume_min"]:
        rejection_reasons.append(f"volume_too_small ({volume:.1f} < {THRESHOLDS['volume_min']})")
    
    if volume > THRESHOLDS["volume_max"]:
        rejection_reasons.append(f"volume_too_large ({volume:.1f} > {THRESHOLDS['volume_max']})")
    
    if hydrophobic_pct < THRESHOLDS["hydrophobic_min"]:
        rejection_reasons.append(f"too_polar ({hydrophobic_pct:.2f} < {THRESHOLDS['hydrophobic_min']})")
    
    if hydrophobic_pct > THRESHOLDS["hydrophobic_max"]:
        rejection_reasons.append(f"too_hydrophobic ({hydrophobic_pct:.2f} > {THRESHOLDS['hydrophobic_max']})")
    
    # Determine status
    if rejection_reasons:
        status = "REJECTED"
    elif confidence >= 0.5 and druggability >= 0.5:
        status = "VALIDATED"
    else:
        status = "CANDIDATE"
    
    return {
        "volume": round(volume, 1),
        "exposure": round(exposure, 3),
        "hydrophobic_pct": round(hydrophobic_pct, 3),
        "residue_count": residue_count,
        "confidence": round(confidence, 3),
        "status": status,
        "rejection_reasons": rejection_reasons
    }
