import json

def load_skeletons():
    with open("entries/096_design_skeletons/pd_l1_v15_skeletons.json") as f:
        return json.load(f)["skeletons"]

def scaffold_family(scaffold):
    if scaffold.startswith("helical"):
        return "helical"
    if scaffold.startswith("beta_hairpin"):
        return "beta_hairpin"
    return scaffold

def guard(candidate):
    skeletons = load_skeletons()
    allowed = [
        (s["buffer"], s["anchor"], s["scaffold"])
        for s in skeletons
    ]

    key = (
        candidate["buffer"],
        candidate["anchor"],
        scaffold_family(candidate["scaffold"])
    )

    return key in allowed

if __name__ == "__main__":
    import sys
    candidate = {
        "buffer": sys.argv[1],
        "anchor": sys.argv[2],
        "scaffold": sys.argv[3]
    }
    if guard(candidate):
        print("✅ GENERATION ALLOWED")
    else:
        print("⛔ GENERATION BLOCKED")
        exit(1)
