import math
class Tier1Laws:
    def _law_155(self, atoms):
        min_d = 999.0
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                a, b = atoms[i], atoms[j]
                if a["chain"] == b["chain"] and abs(a["res_seq"] - b["res_seq"]) <= 1: continue
                d = math.dist(a["pos"], b["pos"])
                min_d = min(min_d, d)
        return ("VETO" if min_d < 2.5 else "PASS"), f"{min_d:.2f}A"
