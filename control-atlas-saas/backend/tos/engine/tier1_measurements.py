import numpy as np
from typing import Tuple, Dict
from .chemistry_registry import is_contextually_bonded

class Tier1Measurements:
    @staticmethod
    def check_law_155_L(structure, threshold: float = 2.0) -> Tuple[str, str, Dict]:
        """
        LAW-155-L: Voxel-Accelerated Steric Clash (Exact O(N)).
        Prunes pairs using a spatial grid and inspects only neighboring voxels.
        """
        all_atoms = structure.atoms + structure.ligands
        if not all_atoms:
            return "PASS", "No atoms.", {}

        grid: Dict[Tuple[int, int, int], list[int]] = {}
        voxel_size = max(threshold * 1.5, 2.5)
        coords = np.array([a.pos for a in all_atoms])
        threshold_sq = threshold * threshold

        for idx in range(len(all_atoms)):
            v_coord = tuple((coords[idx] // voxel_size).astype(int))
            grid.setdefault(v_coord, []).append(idx)

        neighbor_offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]

        for v_coord, atom_indices in grid.items():
            for dx, dy, dz in neighbor_offsets:
                neighbor_v = (v_coord[0] + dx, v_coord[1] + dy, v_coord[2] + dz)
                neighbor_indices = grid.get(neighbor_v)
                if not neighbor_indices:
                    continue

                for i in atom_indices:
                    for j in neighbor_indices:
                        if i >= j:
                            continue

                        a1, a2 = all_atoms[i], all_atoms[j]
                        diff = coords[i] - coords[j]
                        d_sq = float(np.dot(diff, diff))

                        if d_sq >= threshold_sq:
                            continue

                        d = np.sqrt(d_sq)
                        if a1.chain == a2.chain and abs(a1.res_seq - a2.res_seq) <= 1 and d < 1.65:
                            continue
                        if is_contextually_bonded(a1, a2, d):
                            continue
                        if a1.chain == a2.chain and a1.res_seq == a2.res_seq:
                            continue

                        anchor = {"chain": a1.chain, "res_seq": a1.res_seq, "pos": a1.pos}
                        return "VETO", f"Clash: {a1.chain}:{a1.res_name}{a1.res_seq} @ {d:.2f}A", anchor

        return "PASS", f"Verified {len(all_atoms)} atoms via Spatial Grid.", {}

    @staticmethod
    def check_law_160(structure, threshold: float = 4.5) -> Tuple[str, str, Dict]:
        cas = [a for a in structure.atoms if a.atom_name == "CA"]
        for i in range(len(cas) - 1):
            c1, c2 = cas[i], cas[i+1]
            if c1.chain == c2.chain:
                d = np.linalg.norm(np.array(c1.pos) - np.array(c2.pos))
                if d > threshold:
                    return "FAIL", f"Torn chain {c1.res_seq}-{c2.res_seq} @ {d:.2f}A", {"chain": c1.chain, "res_seq": c1.res_seq, "pos": c1.pos}
        return "PASS", "Backbone continuity verified.", {}
