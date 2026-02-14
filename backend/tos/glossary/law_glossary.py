"""
Law glossary — dynamically sourced from the canonical LAW_CANON.
"""
from ..governance.station_sop import LAW_CANON


def list_all_law_ids():
    """Return all law IDs from the canonical source."""
    return sorted(LAW_CANON.keys())


def get_law_explanation(law_id: str) -> dict:
    """Return full explanation for a given law ID."""
    canon = LAW_CANON.get(law_id)
    if not canon:
        return {"error": f"Unknown law: {law_id}"}
    return {
        "law_id": law_id,
        "title": canon["title"],
        "principle": canon["principle"],
        "threshold": canon.get("threshold", "N/A"),
        "unit": canon.get("unit", "N/A")
    }
