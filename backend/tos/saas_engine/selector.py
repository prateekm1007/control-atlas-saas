class ArchitectureSelector:
    """PILLAR 13: 7-Warhead Architecture Lab."""
    CATEGORIES = ["LINEAR", "MULTIVALENT", "ENGAGEMENT", "MOTIFS", "LINKERS", "CONFORMATIONAL", "METAL", "NONE"]

    @staticmethod
    def select(structure, user_intent="NONE"):
        atom_count = len(structure.atoms)
        
        # Internal system guess (for telemetry only)
        if atom_count < 500: derived = "LINEAR"
        elif atom_count < 2000: derived = "CONFORMATIONAL"
        else: derived = "MULTIVALENT"

        # PILLAR 07: If user intent is NONE, we refuse to characterize the warhead.
        is_none = user_intent == "NONE" or user_intent == "NONE OF THE ABOVE"
        
        return {
            "derived_category": derived,
            "user_intent": user_intent,
            "is_override": not is_none and (derived != user_intent),
            # CRITICAL FIX: If NONE, authority is SCAFFOLD, not the system's guess.
            "authoritative_category": "SCAFFOLD (UNCHARACTERIZED)" if is_none else user_intent,
            "derivation_confidence": 0.0 if is_none else 0.85,
            "rationale": "Base scaffold identified. Waiting for design intent." if is_none else "User-asserted design intent."
        }
