"""
WP-SCI-1 regression tests for the corrected LAW-120 measurement (directive §2, §7).

The historical LAW-120 measured O(i)-C(i)-N(i) — the same residue's amide
nitrogen, which is NOT bonded to C(i) — against the Engh-Huber peptide
ideal (123.0°). Definition v2 measures:

  - N(i)-CA(i)-C(i)   intra-residue (unchanged)
  - O(i)-C(i)-N(i+1)  peptide carbonyl angle, NEXT residue's nitrogen,
                      valid sequential pairs only
  - CA(i)-C(i)-N(i+1) distinct peptide angle, sequential pairs only

The exact-geometry tests below construct backbone coordinates with
numpy INDEPENDENTLY of the engine (no import from coord_math) so a
shared bug cannot hide: all measured triples sit exactly on their
Engh-Huber ideals, so a correct engine must report RMSD ~ 0.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Atom, Structure, ConfidenceProfile
from tos.governance.station_sop import IDEAL_TABLE


def _atom(name, element, pos, res_name="ALA", res_seq=1, chain="A",
          icode="", b_iso=90.0):
    return Atom(atom_name=name, element=element,
                pos=np.array(pos, dtype=float), res_name=res_name,
                res_seq=res_seq, chain_id=chain, insertion_code=icode,
                b_iso=b_iso)


def _dir(angle_deg):
    a = math.radians(angle_deg)
    return np.array([math.cos(a), math.sin(a), 0.0])


def _ideal_two_residue_backbone(seq1=1, seq2=2):
    """Two residues whose every LAW-120 triple sits exactly on its ideal.

    Construction (independent of engine code):
      C1 at origin; O1 along +x at 1.231 Å;
      N2 at 1.329 Å from C1 at 123.0° from O1   -> O-C-N ideal
      CA1 at 1.525 Å from C1 at 123+116.6 deg   -> CA-C-N ideal
      N1 at 1.470 Å from CA1 s.t. N-CA-C = 111.0°
      CA2/C2 built off N2 with N-CA-C = 111.0°
    """
    c1 = np.zeros(3)
    o1 = IDEAL_TABLE["C-O"] * _dir(0.0)
    n2 = IDEAL_TABLE["C-N"] * _dir(123.0)
    ca1 = IDEAL_TABLE["CA-C"] * _dir(123.0 + 116.6)
    # direction CA1->C1 is dir(239.6 - 180) = dir(59.6); N1 at 111 deg from it
    n1 = ca1 + IDEAL_TABLE["N-CA"] * _dir(59.6 + 111.0)
    ca2 = n2 + IDEAL_TABLE["N-CA"] * _dir(30.0)
    # CA2->N2 is dir(210); C2 at 111 deg from that direction
    c2 = ca2 + IDEAL_TABLE["CA-C"] * _dir(210.0 - 111.0)
    o2 = c2 + IDEAL_TABLE["C-O"] * _dir(0.0)

    atoms = [
        _atom("N", "N", n1, res_seq=seq1),
        _atom("CA", "C", ca1, res_seq=seq1),
        _atom("C", "C", c1, res_seq=seq1),
        _atom("O", "O", o1, res_seq=seq1),
        _atom("N", "N", n2, res_seq=seq2),
        _atom("CA", "C", ca2, res_seq=seq2),
        _atom("C", "C", c2, res_seq=seq2),
        _atom("O", "O", o2, res_seq=seq2),
    ]
    return atoms


def _independent_angle(a, b, c):
    """Angle at b in degrees — local implementation, no engine import."""
    ba, bc = a - b, c - b
    cosang = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    return math.degrees(math.acos(np.clip(cosang, -1.0, 1.0)))


def _struct(atoms, method="predicted"):
    return Structure(atoms=atoms, audit_id="WPSCI1",
                     confidence=ConfidenceProfile([], False, method))


class TestLaw120CorrectAtomIndices:

    def test_ideal_geometry_yields_zero_rmsd(self):
        """Every measured triple at its ideal -> RMSD ~ 0, measured PASS.

        Under the corrected definition all four triples (2x N-CA-C,
        O-C-N, CA-C-N) match the ideals exactly. Under the historical
        O(i)-C(i)-N(i) definition the same coordinates would produce a
        large RMSD because the intra O-C-N pseudo-angle (~35 deg here)
        is compared against 123.0°.
        """
        atoms = _ideal_two_residue_backbone()
        # independent check of the constructed peptide angles
        assert abs(_independent_angle(atoms[3].pos, atoms[2].pos,
                                      atoms[4].pos) - 123.0) < 1e-6, \
            "fixture construction: O-C-N must be exactly ideal"
        assert abs(_independent_angle(atoms[1].pos, atoms[2].pos,
                                      atoms[4].pos) - 116.6) < 1e-6, \
            "fixture construction: CA-C-N must be exactly ideal"

        results, _, _ = Tier1Measurements.run_full_audit(_struct(atoms))
        res = results["LAW-120"]
        assert res["status"] == "PASS"
        assert res["sample"] == 4  # 2 intra + 1 O-C-N + 1 CA-C-N
        assert res["observed"] == pytest.approx(0.0, abs=0.05)

    def test_ocn_sample_requires_next_residue_nitrogen(self):
        """Removing residue 2's nitrogen removes exactly the O-C-N and
        CA-C-N samples (not replaced by same-residue N)."""
        atoms = _ideal_two_residue_backbone()
        no_n2 = [a for a in atoms if not (a.res_seq == 2 and a.atom_name == "N")]
        results, _, _ = Tier1Measurements.run_full_audit(_struct(no_n2))
        res = results["LAW-120"]
        assert res["status"] == "PASS"
        # only residue 1's intra N-CA-C remains: residue 2 loses its intra
        # angle (no N) AND both peptide angles (no next-residue nitrogen)
        assert res["sample"] == 1

    def test_non_adjacent_residues_are_not_joined(self):
        """A numbering gap (1 -> 3) yields NO peptide-angle samples."""
        atoms = _ideal_two_residue_backbone(seq2=3)
        results, _, _ = Tier1Measurements.run_full_audit(_struct(atoms))
        res = results["LAW-120"]
        assert res["status"] == "PASS"  # the two intra angles still measured
        assert res["sample"] == 2

    def test_insertion_code_successor_counts_as_sequential(self):
        """Residues 5 and 5A are immediate successors; O-C-N is measured."""
        atoms = _ideal_two_residue_backbone(seq1=1, seq2=2)
        remapped = []
        for a in atoms:
            if a.res_seq == 1:
                remapped.append(replace_res(a, seq=5, icode=""))
            else:
                remapped.append(replace_res(a, seq=5, icode="A"))
        results, _, _ = Tier1Measurements.run_full_audit(_struct(remapped))
        res = results["LAW-120"]
        assert res["sample"] == 4

    def test_cross_chain_never_joined(self):
        """The two residues on different chains must produce no peptide
        angle samples."""
        atoms = _ideal_two_residue_backbone()
        atoms = [replace_res(a, chain="B") if a.res_seq == 2 else a
                 for a in atoms]
        results, _, _ = Tier1Measurements.run_full_audit(_struct(atoms))
        res = results["LAW-120"]
        assert res["sample"] == 2


class TestLaw120EmptySample:

    def test_empty_sample_is_not_measured_never_pass(self):
        """A structure with no measurable backbone angle (single residue
        with only C and O) must yield NOT_MEASURED with observed=None —
        the historical code returned observed 0.0 and PASS."""
        atoms = [
            _atom("C", "C", [0.0, 0.0, 0.0], res_seq=1),
            _atom("O", "O", [1.231, 0.0, 0.0], res_seq=1),
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_struct(atoms))
        res = results["LAW-120"]
        assert res["status"] == "NOT_MEASURED"
        assert res["observed"] is None
        assert res["sample"] == 0
        assert "reason" in res

    def test_ca_only_residues_not_measured(self):
        atoms = [
            _atom("CA", "C", [float(i * 3.8), 0.0, 0.0], res_seq=i)
            for i in range(1, 5)
        ]
        results, _, _ = Tier1Measurements.run_full_audit(_struct(atoms))
        res = results["LAW-120"]
        assert res["status"] == "NOT_MEASURED"
        assert res["observed"] is None

    def test_measured_veto_still_possible(self):
        """Genuine distortion is still measured and vetoes at the law
        level (containment is applied downstream by the adjudicator)."""
        atoms = _ideal_two_residue_backbone()
        # stretch the first N-CA-C angle to ~170 deg by moving N1
        ca1 = atoms[1].pos
        c1 = atoms[2].pos
        dir_out = (ca1 - c1) / np.linalg.norm(ca1 - c1)
        # rotate N1 around CA1 to near-collinear with CA1->C1
        atoms[0].pos = ca1 - 1.47 * dir_out
        results, _, _ = Tier1Measurements.run_full_audit(_struct(atoms))
        res = results["LAW-120"]
        assert res["status"] == "VETO"
        assert res["observed"] > 18.0
        assert res["sample"] > 0


# ── helpers ────────────────────────────────────────────────────────────────
def dict_res(atom, icode):
    a = _clone(atom)
    a.insertion_code = icode
    return a


def replace_res(atom, seq=None, icode=None, chain=None):
    a = _clone(atom)
    if seq is not None:
        a.res_seq = seq
    if icode is not None:
        a.insertion_code = icode
    if chain is not None:
        a.chain_id = chain
    return a


def _clone(atom):
    import copy
    return copy.deepcopy(atom)
