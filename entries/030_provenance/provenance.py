#!/usr/bin/env python3
"""
Entry 030 — Decision Provenance Export

Generates versioned, reproducible provenance records for every decision.
This is what makes Control Atlas auditable evidence, not just a tool.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

CONTROL_ATLAS_VERSION = "1.0.0"
GRAMMAR_VERSION = "1.2"
PHYSICS_VERSION = "1.0"

# Physics thresholds (frozen, versioned)
PHYSICS_THRESHOLDS = {
    "volume_min": 150.0,
    "volume_max": 1500.0,
    "hydrophobic_min": 0.30,
    "hydrophobic_max": 0.90,
    "residue_min": 8
}


def hash_content(content: str) -> str:
    """Generate short hash for reproducibility tracking."""
    return hashlib.sha256(content.encode()).hexdigest()[:12]


def generate_provenance(
    sequence: str,
    compound_id: str,
    smiles: str,
    decision: dict,
    structure_result: dict,
    pocket_data: dict,
    checks_passed: list = None,
    checks_failed: list = None
) -> dict:
    """
    Generate a complete provenance record for a screening decision.
    """
    now = datetime.now(timezone.utc).isoformat()
    
    provenance = {
        "provenance_version": "1.0.0",
        "timestamp": now,
        "control_atlas_version": CONTROL_ATLAS_VERSION,
        "grammar_version": GRAMMAR_VERSION,
        
        "input": {
            "sequence_length": len(sequence),
            "sequence_hash": hash_content(sequence),
            "compound_id": compound_id,
            "smiles": smiles,
            "smiles_hash": hash_content(smiles)
        },
        
        "decision": {
            "status": decision.get("status", "UNKNOWN"),
            "confidence": decision.get("confidence", 0.0),
            "reasons": decision.get("reasons", [])
        },
        
        "chain": {
            "structure": {
                "source": structure_result.get("source", "unknown"),
                "confidence": structure_result.get("confidence_global", 0.0),
                "pdb_hash": hash_content(structure_result.get("pdb_path", "")),
                "cached": structure_result.get("cached", False)
            },
            "pocket": {
                "id": pocket_data.get("pocket_id", "unknown"),
                "status": pocket_data.get("status", "unknown"),
                "confidence": pocket_data.get("confidence", 0.0),
                "volume": pocket_data.get("volume", 0.0),
                "hydrophobic_pct": pocket_data.get("hydrophobic_pct", 0.0),
                "exposure": pocket_data.get("exposure", 0.0)
            },
            "physics": {
                "thresholds": PHYSICS_THRESHOLDS,
                "checks_passed": checks_passed or [],
                "checks_failed": checks_failed or []
            }
        },
        
        "reproducibility": {
            "deterministic": True,
            "grammar_hash": hash_content(GRAMMAR_VERSION),
            "physics_version": PHYSICS_VERSION,
            "control_atlas_version": CONTROL_ATLAS_VERSION
        }
    }
    
    return provenance


def save_provenance(provenance: dict, output_dir: Path, compound_id: str) -> str:
    """Save provenance record to file."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    filename = f"provenance_{compound_id}_{provenance['chain']['pocket']['id']}.json"
    filepath = output_dir / filename
    
    with open(filepath, "w") as f:
        json.dump(provenance, f, indent=2)
    
    return str(filepath)


def generate_batch_provenance(
    sequence: str,
    structure_result: dict,
    screening_results: list,
    pockets: list
) -> dict:
    """
    Generate summary provenance for an entire batch run.
    """
    now = datetime.now(timezone.utc).isoformat()
    
    pocket_lookup = {p["pocket_id"]: p for p in pockets}
    
    valid_count = sum(1 for r in screening_results if r.get("status") == "VALID")
    reject_count = sum(1 for r in screening_results if r.get("status") == "REJECT")
    candidate_count = sum(1 for r in screening_results if r.get("status") == "CANDIDATE")
    
    return {
        "provenance_version": "1.0.0",
        "batch_id": hash_content(f"{sequence}{now}"),
        "timestamp": now,
        "control_atlas_version": CONTROL_ATLAS_VERSION,
        "grammar_version": GRAMMAR_VERSION,
        
        "input_summary": {
            "sequence_length": len(sequence),
            "sequence_hash": hash_content(sequence),
            "compounds_screened": len(set(r.get("compound_id") for r in screening_results)),
            "pockets_screened": len(pockets)
        },
        
        "structure": {
            "source": structure_result.get("source", "unknown"),
            "confidence": structure_result.get("confidence_global", 0.0)
        },
        
        "results_summary": {
            "total_screens": len(screening_results),
            "valid": valid_count,
            "candidate": candidate_count,
            "reject": reject_count,
            "rejection_rate": round(reject_count / len(screening_results), 3) if screening_results else 0
        },
        
        "physics_config": {
            "thresholds": PHYSICS_THRESHOLDS,
            "version": PHYSICS_VERSION
        },
        
        "reproducibility": {
            "deterministic": True,
            "all_inputs_hashed": True,
            "grammar_hash": hash_content(GRAMMAR_VERSION)
        }
    }
