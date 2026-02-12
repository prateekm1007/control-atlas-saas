import json

class ClaudeSchemaCompiler:
    """
    PILLAR 11 & 18: Claude acts ONLY as a deterministic typesetter.
    No prose. No interpretation. Pure JSON synthesis.
    """
    def __init__(self):
        self.schema_version = "21.5.4"

    def compile_decision_record(self, raw_data):
        """
        S7: Assembly.
        Input: Deterministic physics/probability data.
        Output: Schema-locked JSON for the Notary.
        """
        # Formats the raw data into the institutional decision schema
        compiled = {
            "version": self.schema_version,
            "verdict": raw_data.get('verdict', {}).get('binary', 'ERROR'),
            "physical_evidence": {
                "score": raw_data.get('verdict', {}).get('physical_score', 0),
                "failed_laws": [l['law_id'] for l in raw_data.get('tier1', {}).get('laws', []) if l.get('status') == 'VETO']
            },
            "strategic_priority": raw_data.get('tier3', {}).get('probability', 0)
        }
        return compiled

_compiler = ClaudeSchemaCompiler()
def get_schema_compiler(): return _compiler
