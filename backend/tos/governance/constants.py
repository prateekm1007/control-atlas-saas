from typing import Final

class AuditConstants:
    OVERRIDE_CONFIDENCE: Final = 0.85
    S8_VETO_TAX: Final = 0.25
    SMALL_THRESHOLD: Final = 500
    MEDIUM_THRESHOLD: Final = 2000
    HASH_LEN: Final = 12
    EPI_BASE_PRIOR: Final = 0.15
    ML_WEIGHT: Final = 0.40

# Coverage threshold for LAW-105 (Reliability Coverage)
LAW_105_THRESHOLD = 70.0
