class BipodGenerator:
    """
    Implements the TwinRod-v2 architecture.
    Stitches an Anchor to a Warhead via a Helical Vertex.
    """
    def __init__(self, anchor="KAWAKKAAAKAEAAKAEAAK"):
        self.anchor = anchor # The CHAMP-005 Chassis
        self.vertex = "PAEAAKP" # Proline-capped helical vertex (Rigid)

    def generate(self, warhead_sequence):
        """
        Creates a Bipodal FASTA.
        """
        bipod = f"{self.anchor}{self.vertex}{warhead_sequence}"
        return bipod

if __name__ == "__main__":
    gen = BipodGenerator()
    # Example: Grafting a PD-L1 distal warhead
    print(f"Generated Bipod: {gen.generate('YWPTG')}")
