import logging
from typing import Dict, Any

logger = logging.getLogger("toscanini.governance")

def unify_architecture_authority(arch: Dict[str, Any]) -> Dict[str, Any]:
    is_override = arch.get("is_override", False)
    if is_override:
        arch["authoritative_category"] = arch.get("user_intent", "SCAFFOLD")
        arch["authority_reason"] = "user_asserted_override"
    else:
        arch["authoritative_category"] = arch.get("derived_category", "SCAFFOLD")
        arch["authority_reason"] = "physics_derived"
    return arch
