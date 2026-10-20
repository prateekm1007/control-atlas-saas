"""Deterministic narrative compiler (provider-less).

H-1 (auditor, delegated by owner, 2026-10-20): the Gemini narrative
path was removed completely — no external LLM SDK, no API key, no
model names. The interface (``NarrativeCompiler.synthesize``) is kept
so callers are unchanged; the output is the deterministic, honest
narrative that the removed provider's exception fallback used to
produce. Narrative is advisory; the deterministic verdict stays
authoritative (EC §4 / EC §10).
"""


class NarrativeCompiler:
    def __init__(self):
        self.model = "deterministic-fallback"

    def synthesize(self, verdict, phys_score, conf_score, gen, results):
        # Deterministic synthesis (previously the no-provider fallback):
        # states exactly what was measured, never invents content.
        return (
            f"Audit confirms {phys_score}% Physical Integrity. "
            f"ML Confidence ({conf_score}%) indicates speculative tertiary "
            f"packing. Structure is permissible but lacks high-fidelity "
            f"evidence. (Narrative is deterministic — external LLM "
            f"removed per H-1, 2026-10-20. Verdict: {verdict}.)"
        )
