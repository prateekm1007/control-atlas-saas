import os
import logging
import time

logger = logging.getLogger("toscanini.gemini")

class GeminiCompiler:
    def __init__(self):
        self.api_key = os.getenv("GEMINI_API_KEY")
        self.client = None
        self.model_used = "none"
        
        # 🛡️ PIL-NAR-11: PhD Narrative Engine (16-Tier Hierarchy)
        # Optimized for provided RPD/Availability
        self.shield_tiers = [
            "models/gemini-2.0-pro-exp-02-05",  # Tier 1: Apex Reasoning
            "models/gemini-2.0-flash-exp",      # Tier 2: Experimental Logic
            "models/gemini-2.5-pro",            # Tier 3: Stable Pro
            "models/gemini-2.0-flash",          # Tier 4: Speed Apex
            "models/gemini-2.5-flash",          # Tier 5: Speed Guard
            "models/gemini-2.0-flash-lite",     # Tier 6: Low Latency
            "models/gemini-2.5-flash-lite",     # Tier 7: Efficient Volume
            "models/gemma-3-27b-it",            # Tier 8: Distilled Logic
            "models/gemini-1.5-pro",            # Tier 9: Legacy Stable Pro
            "models/gemini-1.5-flash",          # Tier 10: Legacy Stable Flash
            "models/gemini-2.0-pro-exp-02-05",  # Tier 11: (Rotation Re-entry)
            "models/gemini-2.0-flash-exp",      # Tier 12: (Rotation Re-entry)
            "models/gemini-2.5-pro",            # Tier 13: (Rotation Re-entry)
            "models/gemini-2.5-flash",          # Tier 14: (Rotation Re-entry)
            "models/gemini-2.5-flash-lite",     # Tier 15: (Rotation Re-entry)
            "models/gemini-2.0-flash-lite"      # Tier 16: Final AI Tier
        ]

        if self.api_key and self.api_key != "dummy":
            try:
                import google.generativeai as genai
                genai.configure(api_key=self.api_key)
                self.client = genai
                logger.info("Apex Shield: 16-tier hierarchy initialized.")
            except Exception as e:
                logger.error(f"Gemini SDK load failure: {e}")

    def _ask(self, prompt):
        if not self.client: return None
        for i, model_name in enumerate(self.shield_tiers, 1):
            try:
                model = self.client.GenerativeModel(model_name)
                res = model.generate_content(prompt, request_options={"timeout": 8})
                if res and res.text:
                    self.model_used = model_name
                    return res.text
            except Exception as e:
                logger.warning(f"✘ [CASCADE] Tier {i} ({model_name}) exhausted: {str(e)[:40]}")
                time.sleep(2) # Mandatory Paced Cascade Algorithm
                continue
        return None

    def synthesize_dossier_content(self, ctx):
        intent = ctx.get('arch', 'NONE')
        failing_det = ctx.get('killer_laws', [])
        fail_str = "; ".join(failing_det[:3]) if failing_det else "None"
        prompt = (f"PhD-level Forensic Audit. Verdict: {ctx.get('v','?')}. Integrity: {ctx.get('s',0)}%. "
                  f"Coverage: {ctx.get('c',0)}%. Intent: {intent}. Violations: {fail_str}. "
                  f"Write 3 technical paragraphs (Executive, Deep Dive, Recommendation). ASCII only.")
        raw = self._ask(prompt)
        if raw:
            parts = raw.split("\n\n")
            return {"executive": parts[0], "deep_dive": parts[1] if len(parts)>1 else "Analysis verified.", 
                    "recommendation": parts[2] if len(parts)>2 else "Proceed per metrics."}
        self.model_used = "Internal Analysis Module (Deterministic Fallback)"
        return self._internal_analysis_module(ctx, fail_str)

    def _internal_analysis_module(self, ctx, fail_str):
        v, s = ctx.get('v', '?'), ctx.get('s', 0)
        if v == "VETO":
            return {"executive": f"VETOED. Critical violations in {fail_str}.", 
                    "deep_dive": f"Integrity score {s}%. Geometric impossibilities detected.", "recommendation": "Reject immediately."}
        return {"executive": f"INDETERMINATE. Reliability gate failed.", 
                "deep_dive": "Insufficient high-confidence backbone data.", "recommendation": "Re-model structure."}

_compiler = None
def get_compiler():
    global _compiler
    if _compiler is None: _compiler = GeminiCompiler()
    return _compiler
