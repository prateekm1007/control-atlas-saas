"""
Entry 036 — Unified Atlas Architecture
Defines the interfaces for the Multi-Layer Navigation System.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Dict, Optional

@dataclass
class Constraint:
    """Represents a restriction on the search space."""
    layer: str  # "Physics", "Chemistry", "Biology"
    parameter: str # e.g. "Volume", "LogP", "Expression"
    threshold: float
    actual: float
    status: str # "PASS", "FAIL", "WARNING"

@dataclass
class NavigationResult:
    """The output of the Navigation Engine."""
    status: str # "CLEARED", "BLOCKED"
    confidence: float
    constraints: List[Constraint] = field(default_factory=list)
    path_advice: str = ""
    provenance: Dict = field(default_factory=dict)

class Layer(ABC):
    """Abstract base class for all layers."""
    
    @abstractmethod
    def evaluate(self, context: Dict) -> Constraint:
        """
        Evaluate the context against this layer's rules.
        Returns a Constraint object (Pass/Fail).
        """
        pass

class PhysicsLayer(Layer):
    """
    Layer 1: The Hard Veto.
    Checks geometric existence and stability.
    """
    pass

class ChemistryLayer(Layer):
    """
    Layer 2: Tractability.
    Checks if matter can exist comfortably in the geometric space.
    """
    pass

class BiologyLayer(Layer):
    """
    Layer 3: Relevance.
    Checks if the target matters in the disease context.
    """
    pass

class MathLayer:
    """
    Layer 4: Integration.
    Calculates robustness of the remaining path.
    """
    pass
