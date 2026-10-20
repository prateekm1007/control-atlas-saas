import json

def load_ci():
    with open("entries/095_complement/complement_index_v1.json") as f:
        return json.load(f)["constraints"]

def feasible(candidate):
    ci = load_ci()
    target = candidate["target"]

    if target not in ci:
        return False, "No positive constraints defined"

    rules = ci[target]

    # Scaffold geometry check
    if candidate.get("scaffold") not in rules["allowed_scaffolds"]:
        return False, "Scaffold geometry not allowed"

    # Anchor chemistry (placeholder for later parsing)
    return True, "Candidate within feasible envelope"

if __name__ == "__main__":
    import sys
    candidate = {
        "target": sys.argv[1],
        "scaffold": sys.argv[2]
    }
    ok, reason = feasible(candidate)
    if ok:
        print("✅ FEASIBLE")
        print(" ", reason)
    else:
        print("⛔ INFEASIBLE")
        print(" ", reason)
