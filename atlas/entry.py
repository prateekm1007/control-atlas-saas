class AtlasEntry:
    def __init__(self, target, variant, ccs_map, provenance):
        self.target = target
        self.variant = variant
        self.ccs_map = ccs_map
        self.provenance = provenance

    def top_edges(self, threshold=0.6):
        return {k: v for k, v in self.ccs_map.items() if v >= threshold}
