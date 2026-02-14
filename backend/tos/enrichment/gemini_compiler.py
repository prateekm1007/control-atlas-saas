import os, logging, time
logger = logging.getLogger("toscanini.gemini")

class GeminiCompiler:
    def __init__(self):
        self.api_key = os.getenv("GEMINI_API_KEY")
        self.client = None
        self.model_used = "none"
        # SPEED-OPTIMIZED: 3 tiers only, fastest-first for responsiveness
        # Deep models available but only tried if fast ones fail
        self.fast_tier = [
            "models/gemini-2.0-flash",
            "models/gemini-2.5-flash",
            "models/gemini-2.5-pro",
        ]
        if self.api_key and self.api_key != "dummy":
            try:
                import google.generativeai as genai
                genai.configure(api_key=self.api_key)
                self.client = genai
                logger.info("Gemini Compiler: 3-tier fast hierarchy active.")
            except Exception as e:
                logger.error(f"Gemini API offline: {e}")

    def _ask(self, prompt):
        if not self.client:
            return None
        for model_name in self.fast_tier:
            try:
                model = self.client.GenerativeModel(model_name)
                res = model.generate_content(prompt, request_options={"timeout": 8})
                if res and res.text:
                    logger.info(f"✔ Response from {model_name}")
                    self.model_used = model_name
                    return res.text
            except Exception as e:
                logger.warning(f"✘ {model_name}: {e}")
                time.sleep(0.5)
                continue
        return None

    def synthesize_dossier_content(self, ctx):
        intent = ctx.get('arch', 'NONE')
        failing = ctx.get('failing_laws', [])
        fail_str = "; ".join([f"{l['law_id']}: {l['title']}" for l in failing[:3]]) if failing else "None"

        prompt = (
            f"You are a structural biology PhD. Write a clinical assessment.\n"
            f"Structure verdict: {ctx.get('v','?')}. Physical integrity: {ctx.get('s',0)}%.\n"
            f"Confidence: {ctx.get('c',0)}. EPI priority: {ctx.get('p',0)}%.\n"
            f"Architecture: {intent}. Failing laws: {fail_str}.\n"
            f"Write exactly 3 paragraphs: 1) Executive summary 2) Structural detail 3) Recommendation.\n"
            f"Clinical tone. No marketing. ASCII only. Max 400 words total."
        )
        raw = self._ask(prompt)
        if raw:
            parts = raw.split("\n\n")
            return {
                "executive": parts[0] if len(parts) > 0 else "Analysis complete.",
                "deep_dive": parts[1] if len(parts) > 1 else "Assessment verified.",
                "recommendation": parts[2] if len(parts) > 2 else "Proceed per metrics."
            }
        # Deterministic fallback — no AI needed
        self.model_used = "Internal analysis module (deterministic fallback)"
        v = ctx.get('v', '?')
        s = ctx.get('s', 0)
        c = ctx.get('c', 0)
        if v == "VETO":
            return {
                "executive": f"Structure VETOED. Deterministic physical integrity: {s}%. One or more physical invariant laws violated.",
                "deep_dive": f"Failing laws: {fail_str}. Confidence: {c}. Structure cannot advance to wet-lab consideration.",
                "recommendation": "Do not allocate resources. Address physical violations and resubmit."
            }
        return {
            "executive": f"Structure PASSED all deterministic invariants. Physical integrity: {s}%. Confidence: {c}.",
            "deep_dive": f"All deterministic laws satisfied. Architecture: {intent}. No structural violations detected.",
            "recommendation": "Structure eligible for prioritized wet-lab evaluation."
        }

_compiler = None
def get_compiler():
    global _compiler
    if _compiler is None:
        _compiler = GeminiCompiler()
    return _compiler
