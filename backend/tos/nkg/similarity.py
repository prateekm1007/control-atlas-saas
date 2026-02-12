from typing import Dict, Any, List

class SimilarityEngine:
    """
    PILLAR 16: Deterministic similarity for NKG lookups.
    Matches current candidates against historical wet-lab failures.
    """
    WEIGHTS = {
        "architecture": 0.30,
        "tier1_profile": 0.25,
        "ca_span": 0.15,
        "aspect_ratio": 0.15,
        "secondary_structure": 0.15
    }

    @staticmethod
    def calculate_jaccard(list1: List[str], list2: List[str]) -> float:
        s1, s2 = set(list1), set(list2)
        if not s1 and not s2: return 1.0
        return len(s1 & s2) / len(s1 | s2)

    @staticmethod
    def normalized_delta(v1: float, v2: float) -> float:
        if v1 == 0 and v2 == 0: return 1.0
        return 1.0 - (abs(v1 - v2) / max(v1, v2, 1.0))

    @classmethod
    def compute_similarity(cls, f1: Dict[str, Any], f2: Dict[str, Any]) -> float:
        """Computes S = Σ (wi * sim_i)"""
        try:
            s_arch = 1.0 if f1['architecture'] == f2['architecture'] else 0.0
            s_t1 = cls.calculate_jaccard(f1['t1_profile'], f2['t1_profile'])
            s_span = cls.normalized_delta(f1['ca_span'], f2['ca_span'])
            s_ratio = cls.normalized_delta(f1['aspect_ratio'], f2['aspect_ratio'])
            
            # Secondary structure (Helix/Sheet/Coil) L1 distance
            ss1, ss2 = f1['sec_str'], f2['sec_str']
            l1_dist = abs(ss1['h']-ss2['h']) + abs(ss1['s']-ss2['s']) + abs(ss1['c']-ss2['c'])
            s_ss = max(0.0, 1.0 - (l1_dist / 2.0))

            score = (
                cls.WEIGHTS["architecture"] * s_arch +
                cls.WEIGHTS["tier1_profile"] * s_t1 +
                cls.WEIGHTS["ca_span"] * s_span +
                cls.WEIGHTS["aspect_ratio"] * s_ratio +
                cls.WEIGHTS["secondary_structure"] * s_ss
            )
            return round(score, 3)
        except KeyError:
            return 0.0
