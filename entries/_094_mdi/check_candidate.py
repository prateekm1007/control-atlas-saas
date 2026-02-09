from .load_mdi import load_mdi
from .scaffold_registry import classify_scaffold

def normalize(s):
    return s.lower().replace("-", "").replace("_", "") if s else ""

def law_applies(law, candidate):
    applies = law.get("applies_to", {})
    tgt = normalize(candidate.get("target"))
    gen = normalize(candidate.get("generator"))

    if normalize(applies.get("target")) != tgt:
        return False

    if "generator" in applies and normalize(applies["generator"]) != gen:
        return False

    if "scaffold_geometry" in applies:
        info = classify_scaffold(candidate.get("scaffold"))
        if not info or info.get("geometry") != applies["scaffold_geometry"]:
            return False

    if "anchor_requires" in applies:
        # Anchor inspection can be wired later
        return True

    return True

def violates_mdi(candidate):
    laws = sorted(load_mdi(), key=lambda l: l.get("priority", 0), reverse=True)
    for law in laws:
        if law_applies(law, candidate):
            return law
    return None

def all_violations(candidate):
    """
    Returns a list of ALL violated laws, sorted by priority.
    """
    violations = []
    laws = sorted(load_mdi(), key=lambda l: l.get("priority", 0), reverse=True)

    for law in laws:
        if law_applies(law, candidate):
            violations.append(law)

    return violations
