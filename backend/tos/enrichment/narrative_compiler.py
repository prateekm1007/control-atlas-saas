"""
TOSCANINI Narrative Compiler — deterministic, provider-less (v22.5.3)
PIL-NAR-07: Non-blocking narrative synthesis.

H-1 (auditor, delegated by owner, 2026-10-20): the Gemini narrative
cascade was REMOVED completely — no external LLM SDK, no API key, no
model names, no environment entries. This module keeps the exact
interface the API contract requires (`get_compiler`, `.model_used`,
`.synthesize_dossier_content`) and produces only the deterministic,
honest fallback narrative:

- EC §4 / EC §10 compliance: narrative output is advisory, clearly
  separated from the authoritative deterministic verdict, and provider
  unavailability degrades honestly (there is no provider to become
  unavailable — the narrative layer states exactly what it is).
- The engine verdict, measurements, and provenance are unchanged.
"""
import logging

logger = logging.getLogger("toscanini.narrative")

NARRATIVE_LAYER = "deterministic-fallback"


class NarrativeCompiler:
    """Deterministic narrative compiler (no external provider).

    Interface-compatible with the removed Gemini cascade compiler:
    `synthesize_dossier_content(context)` returns the three-section
    narrative dict, and `model_used` reports the narrative layer that
    produced it.
    """

    def __init__(self):
        self._model_used = (
            "Deterministic Narrative Compiler "
            "(no external LLM — Gemini removed per H-1, 2026-10-20)")

    @property
    def model_used(self) -> str:
        return self._model_used

    def synthesize_dossier_content(self, context: dict) -> dict:
        """Deterministic, honest narrative from the audit context alone."""
        return self._build_narrative(context)

    def _build_narrative(self, ctx: dict) -> dict:
        verdict = ctx.get("v", "UNKNOWN")
        score = ctx.get("s", 0)
        coverage = ctx.get("c", 0)
        killers = ctx.get("killer_laws", [])
        arch = ctx.get("arch", "NONE")
        return {
            "executive": f"Deterministic audit complete. Verdict: {verdict}. "
                         f"Score: {score}% with {coverage}% coverage.",
            "deep_dive": f"Architecture class: {arch}. "
                         + (f"Failing laws: {', '.join(killers)}."
                            if killers else "No deterministic failures detected."),
            "recommendation": "Narrative layer is deterministic (no external "
                              "LLM is configured in this build). Human review "
                              "recommended. The deterministic verdict "
                              "is authoritative regardless of narrative availability."
        }


_instance = None

def get_compiler() -> NarrativeCompiler:
    global _instance
    if _instance is None:
        _instance = NarrativeCompiler()
    return _instance
