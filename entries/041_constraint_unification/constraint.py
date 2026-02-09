from dataclasses import dataclass, field
from typing import Dict, Optional

@dataclass
class Constraint:
    layer: str              # "Physics" | "Chemistry" | "Biology"
    parameter: str          # e.g. "Volume", "LogP", "Essentiality"
    threshold: float        # Rejection boundary
    actual: float           # Observed value
    margin: float           # Normalized distance to threshold (>0 PASS, <0 FAIL)
    status: str             # "PASS" | "FAIL" | "WARNING"
    collapses_space: bool   # True = hard veto
    provenance: Dict = field(default_factory=dict) # Evidence, source, metadata
