from ..governance.constants import AuditConstants

class StrategicPredictor:
    @staticmethod
    def calculate(phys_score, confidence, weight_deriv):
        base = AuditConstants.EPI_BASE_PRIOR
        physics_bonus = 0.25 if phys_score == 100 else 0.0
        ml = (confidence / 100) * AuditConstants.ML_WEIGHT
        s6 = (base + physics_bonus + ml) * weight_deriv
        s8 = 0.0
        p_final = max(0.05, min(0.95, s6 * (1 - s8)))
        return {"probability": round(p_final * 100, 1), "certainty": "HIGH" if p_final > 0.6 else "EMERGING",
                "s6_raw": round(s6 * 100, 1), "s8_penalty": s8}
