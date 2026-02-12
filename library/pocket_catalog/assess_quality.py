#!/usr/bin/env python3
"""
Pocket Quality Assessment
Classifies pockets based on physics thresholds
"""

import os
import json

CATALOG = os.path.expanduser("~/control-atlas/library/pocket_catalog")

# Quality Thresholds (derived from KRAS reference + literature)
THRESHOLDS = {
    "volume_A3": {
        "min": 300,      # Too small = not druggable
        "max": 5000,     # Too large = not specific
        "optimal": (800, 3000)
    },
    "exposure": {
        "min": 10,       # Too buried = not accessible
        "max": 50,       # Too exposed = not selective
        "optimal": (15, 40)
    },
    "hydrophobic_pct": {
        "min": 10,       # Need some hydrophobic anchor
        "max": 60,       # Too hydrophobic = solubility issues
        "optimal": (15, 50)
    },
    "atom_count": {
        "min": 50        # Minimum for meaningful pocket
    }
}

def assess_quality(metrics):
    """
    Returns: (status, score, issues)
    Status: VALIDATED, CANDIDATE, REJECTED
    Score: 0-100
    """
    if metrics.get("status") != "computed":
        return "REJECTED", 0, ["Physics not computed"]
    
    issues = []
    score = 100
    
    vol = metrics.get("volume_A3", 0)
    exp = metrics.get("exposure", 0)
    hydro = metrics.get("hydrophobic_pct", 0)
    atoms = metrics.get("atom_count", 0)
    
    # Volume checks
    if vol < THRESHOLDS["volume_A3"]["min"]:
        issues.append(f"Volume too small ({vol} < {THRESHOLDS['volume_A3']['min']})")
        score -= 30
    elif vol > THRESHOLDS["volume_A3"]["max"]:
        issues.append(f"Volume too large ({vol} > {THRESHOLDS['volume_A3']['max']})")
        score -= 20
    elif not (THRESHOLDS["volume_A3"]["optimal"][0] <= vol <= THRESHOLDS["volume_A3"]["optimal"][1]):
        score -= 10
    
    # Exposure checks
    if exp < THRESHOLDS["exposure"]["min"]:
        issues.append(f"Too buried ({exp} < {THRESHOLDS['exposure']['min']})")
        score -= 25
    elif exp > THRESHOLDS["exposure"]["max"]:
        issues.append(f"Too exposed ({exp} > {THRESHOLDS['exposure']['max']})")
        score -= 15
    elif not (THRESHOLDS["exposure"]["optimal"][0] <= exp <= THRESHOLDS["exposure"]["optimal"][1]):
        score -= 5
    
    # Hydrophobicity checks
    if hydro < THRESHOLDS["hydrophobic_pct"]["min"]:
        issues.append(f"Low hydrophobicity ({hydro} < {THRESHOLDS['hydrophobic_pct']['min']})")
        score -= 20
    elif hydro > THRESHOLDS["hydrophobic_pct"]["max"]:
        issues.append(f"Too hydrophobic ({hydro} > {THRESHOLDS['hydrophobic_pct']['max']})")
        score -= 15
    
    # Atom count
    if atoms < THRESHOLDS["atom_count"]["min"]:
        issues.append(f"Insufficient atoms ({atoms})")
        score -= 25
    
    # Classify
    score = max(0, score)
    if score >= 80 and len(issues) == 0:
        status = "VALIDATED"
    elif score >= 50:
        status = "CANDIDATE"
    else:
        status = "REJECTED"
    
    return status, score, issues

def assess_catalog():
    print("=== POCKET QUALITY ASSESSMENT ===")
    print("")
    print(f"{'TARGET':<20} {'STATUS':<12} {'SCORE':>6} {'ISSUES'}")
    print("-" * 70)
    
    results = {"VALIDATED": [], "CANDIDATE": [], "REJECTED": []}
    
    for entry in sorted(os.listdir(CATALOG)):
        physics_path = os.path.join(CATALOG, entry, "physics_metrics.json")
        frame_path = os.path.join(CATALOG, entry, "pocket_frame.json")
        
        if not os.path.exists(physics_path):
            continue
        
        with open(physics_path) as f:
            metrics = json.load(f)
        
        status, score, issues = assess_quality(metrics)
        results[status].append(entry)
        
        issue_str = "; ".join(issues) if issues else "-"
        print(f"{entry:<20} {status:<12} {score:>6} {issue_str}")
        
        # Update pocket frame with quality status
        if os.path.exists(frame_path):
            with open(frame_path) as f:
                frame = json.load(f)
            frame["quality"]["physics_consistency"] = score / 100
            frame["quality"]["status"] = status
            frame["quality"]["issues"] = issues
            with open(frame_path, "w") as f:
                json.dump(frame, f, indent=2)
    
    print("")
    print("=== SUMMARY ===")
    print(f"VALIDATED:  {len(results['VALIDATED'])}")
    print(f"CANDIDATE:  {len(results['CANDIDATE'])}")
    print(f"REJECTED:   {len(results['REJECTED'])}")
    
    # Save thresholds
    thresh_path = os.path.join(CATALOG, "quality_thresholds.json")
    with open(thresh_path, "w") as f:
        json.dump(THRESHOLDS, f, indent=2)
    print(f"\nThresholds saved to {thresh_path}")

if __name__ == "__main__":
    assess_catalog()
