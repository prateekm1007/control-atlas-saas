"""
Law glossary — dynamically sourced from the law registry (WP-8).

Authority chain: station_sop.LAW_CANON (runtime, hash-pinned) ->
governance/law_registry.yaml (documentation) -> this glossary.

The registry is YAML; PyYAML is only needed when this module is used at
runtime. Without it, explanations fall back to the pinned LAW_CANON so
the product never breaks. Tests always parse the registry itself.
"""
from ..governance.station_sop import LAW_CANON

_REGISTRY_PATHS = (
    "governance/law_registry.yaml",
    "../governance/law_registry.yaml",
)


def _load_registry():
    """Load governance/law_registry.yaml if PyYAML is available; else None."""
    from pathlib import Path
    try:
        import yaml  # provided by requirements-dev.txt / PyYAML
    except ImportError:
        return None
    here = Path(__file__).resolve()
    for rel in _REGISTRY_PATHS:
        p = (here.parent / rel).resolve()
        if not p.exists():
            # also try relative to the repo root (backend/ cwd)
            p = (here.parent.parent.parent / rel).resolve()
        if p.exists():
            with open(p, "r", encoding="utf-8") as f:
                return yaml.safe_load(f)
    return None


def list_all_law_ids():
    """Return all law IDs from the canonical source (enforced + retired)."""
    reg = _load_registry()
    if reg and "laws" in reg:
        return sorted(entry["id"] for entry in reg["laws"])
    return sorted(LAW_CANON.keys())


def get_law_explanation(law_id: str) -> dict:
    """Return a customer-safe explanation for a given law ID.

    Enforced laws: registry text (what the code actually measures),
    with pinned threshold/unit values from LAW_CANON as fallback.
    Retired laws: explicit retirement notice.
    """
    reg = _load_registry()
    entry = None
    if reg and "laws" in reg:
        entry = next((e for e in reg["laws"] if e["id"] == law_id), None)

    if entry is not None:
        if entry.get("status") == "retired":
            return {
                "law_id": law_id,
                "title": entry["title"],
                "status": "retired",
                "explanation": (entry.get("measures") or "").strip(),
            }
        canon = LAW_CANON.get(law_id, {})
        return {
            "law_id": law_id,
            "title": entry["title"],
            "status": entry.get("status", "enforced"),
            "tier": entry.get("tier"),
            "principle": (entry.get("measures") or "").strip(),
            "threshold": canon.get("threshold", entry.get("threshold")),
            "unit": canon.get("unit", entry.get("units")),
        }

    # Fallback: pinned canon only (no registry available).
    canon = LAW_CANON.get(law_id)
    if not canon:
        return {"error": f"Unknown law: {law_id}"}
    return {
        "law_id": law_id,
        "title": canon["title"],
        "principle": "See governance/law_registry.yaml for the authoritative description.",
        "threshold": canon.get("threshold", "N/A"),
        "unit": canon.get("unit", "N/A"),
    }
