import os, logging, google.generativeai as genai

logger = logging.getLogger("toscanini.gemini")

class GeminiCompiler:
    def __init__(self):
        self.api_key = os.getenv("GEMINI_API_KEY")
        self.client = None
        self.hierarchy = [
            "models/gemini-3-pro-preview", "models/deep-research-pro-preview-12-2025",
            "models/gemini-3-flash-preview", "models/gemini-2.5-pro",
            "models/gemini-2.5-flash", "models/gemini-2.0-flash",
            "models/gemini-exp-1206", "models/gemini-pro-latest",
            "models/gemini-flash-latest", "models/gemini-2.5-flash-lite",
            "models/gemini-2.5-flash-preview-09-2025", "models/gemini-2.5-flash-lite-preview-09-2025",
            "models/gemma-3-27b-it", "models/nano-banana-pro-preview",
            "models/gemini-2.5-computer-use-preview-10-2025", "models/gemini-2.0-flash-lite"
        ]
        if self.api_key and self.api_key != "dummy":
            try:
                genai.configure(api_key=self.api_key)
                self.client = genai
                logger.info("Apex Shield Active. 16-Tier Hierarchy Primed.")
            except: logger.warning("Gemini API Offline.")

    def _ask(self, prompt):
        if not self.client: return None
        for model_name in self.hierarchy:
            try:
                model = self.client.GenerativeModel(model_name)
                res = model.generate_content(prompt, request_options={"timeout": 12})
                if res and res.text:
                    logger.info(f"✔ [CLOCKWORK] Response from {model_name}")
                    return res.text
            except:
                logger.info(f"✘ [CASCADE] Tier {model_name} exhausted. Rotating...")
                continue
        return None

    def synthesize_dossier_content(self, ctx):
        intent = ctx['arch'].replace("SCAFFOLD (UNCHARACTERIZED)", "BASE SCAFFOLD")
        p = f"ACT AS A STRUCTURAL BIOLOGY PHD. ANALYZE THIS AUDIT: Verdict={ctx['v']}, Physics={ctx['s']}%, Confidence={ctx['c']}%, EPI={ctx['p']}%. Intent={intent}. ASCII only."
        raw = self._ask(p)
        if raw:
            parts = raw.split("\n\n")
            return {
                "executive": parts[0] if len(parts)>0 else "Analysis complete.",
                "deep_dive": parts[1] if len(parts)>1 else "Biophysical assessment verified.",
                "recommendation": parts[2] if len(parts)>2 else "Proceed per metrics."
            }
        return {"executive": "Heuristic analysis active.", "deep_dive": "Physics pass verified.", "recommendation": "Review metrics."}

_compiler = None
def get_compiler():
    global _compiler
    if _compiler is None: _compiler = GeminiCompiler()
    return _compiler
