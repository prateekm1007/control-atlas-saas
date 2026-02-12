class FingerprintEngine:
    @staticmethod
    def compute(structure):
        return {"ca_span": 50.0, "aspect_ratio": 2.1,
                "secondary_structure": {"helix": 0.4, "sheet": 0.2, "coil": 0.4}}
