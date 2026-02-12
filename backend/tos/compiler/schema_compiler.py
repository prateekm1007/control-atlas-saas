import json
from typing import Dict, Any, Optional, List
from pathlib import Path

class SchemaCompiler:
    def __init__(self):
        self.contracts_dir = Path(__file__).parent.parent / "contracts"

    def compile_executive_summary(self, v_binary, phys_score, conf_score):
        """Typesetting truth: No new claims permitted."""
        if v_binary == "PASS":
            return f"Audit confirms {phys_score}% Physical Integrity. Adherence to 10-Law Canon is absolute. Model Confidence ({conf_score}%) suggests high conformational credibility."
        return f"Audit VETOED. Physical Integrity ({phys_score}%) indicates critical violations. Biological synthesis is impossible as coordinates represent infinite energy gradients."

    def normalize_s8_report(self, text):
        """Maps free-text to taxonomy keywords."""
        text_lower = text.lower()
        with open(self.contracts_dir / "failure_taxonomy.json") as f:
            tax = json.load(f)["failure_classes"]
        for cls, keywords in tax.items():
            if any(k in text_lower for k in keywords): return cls
        return "UNKNOWN_PATHOLOGY"

_compiler = SchemaCompiler()
def get_compiler(): return _compiler
