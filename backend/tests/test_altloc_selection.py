"""
P0 item 6 (2026-10-19 audit): alternate-conformer selection tests.

Rule (one conformer per atom, applied everywhere): HIGHEST occupancy
wins; ties keep the FIRST record in file order.

Before this fix the three atom-lookup paths disagreed:
  - IngestionParser kept every record (LAW-130/155/200 double-counted
    altloc atoms);
  - coord_math.get_atom_pos took the FIRST match;
  - coord_math.get_residue_atoms and _extract kept the LAST.

Real bytes: 3ERD carries 66 altloc atom records (ASP A 321 backbone in
conformers A and B at 0.50/0.50 occupancy — a genuine tie-break case;
plus ARG A 363, CYS B 381, MET B 437 side chains). Pinned fixture
sha256 7a192385c93b90670a52e8fc478b441152a49811b93e6ca202c9a487cb4e7418.
"""
from pathlib import Path

import numpy as np

from tos.ingestion.processor import IngestionProcessor, Atom
from tos.engine.coord_math import get_atom_pos, get_residue_atoms
from tos.engine.tier1_measurements import Tier1Measurements

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_3erd.pdb"


def _parse():
    return IngestionProcessor.run(FIXTURE.read_bytes(), "3erd.pdb",
                                  "Audit", "experimental")


# ── Real 3ERD ───────────────────────────────────────────────────────────

def test_3erd_one_atom_per_name_after_selection():
    structure = _parse()
    keys = [(a.chain_id, a.res_seq, a.insertion_code, a.atom_name)
            for a in structure.atoms]
    assert len(keys) == len(set(keys)), "duplicate atom records survived"


def test_3erd_asp321_tie_keeps_first_conformer():
    structure = _parse()
    # ASP A 321 has N/CA/C/O in conformers A and B, both occupancy 0.50.
    # Tie -> first record in file order (conformer A) must win.
    n_pos = get_atom_pos(structure.atoms, "A", 321, "N", "")
    assert n_pos is not None
    assert np.allclose(n_pos, [30.044, -0.841, -5.710])  # the 'A' record


def test_3erd_full_audit_runs_on_selected_conformers():
    structure = _parse()
    results, coverage, _ = Tier1Measurements.run_full_audit(structure)
    assert results["LAW-120"]["status"] in ("PASS", "VETO")
    assert results["LAW-120"]["sample"] > 0


# ── Synthetic selection-rule unit tests ─────────────────────────────────

def _atom(name, x, occ, altloc=" ", b=20.0):
    return Atom(atom_name=name, element=name[0], pos=np.array([x, 0.0, 0.0]),
                res_name="ALA", res_seq=1, chain_id="A", insertion_code="",
                b_iso=b, altloc=altloc, occupancy=occ)


def test_highest_occupancy_wins_even_if_later_in_file():
    atoms = [_atom("CB", 1.0, 0.3, altloc="A"),
             _atom("CB", 9.0, 0.7, altloc="B")]
    pos = get_atom_pos(atoms, "A", 1, "CB", "")
    assert pos[0] == 9.0  # occupancy 0.7 beats first-in-file 0.3


def test_tie_keeps_first_record():
    atoms = [_atom("CB", 1.0, 0.5, altloc="A"),
             _atom("CB", 9.0, 0.5, altloc="B")]
    pos = get_atom_pos(atoms, "A", 1, "CB", "")
    assert pos[0] == 1.0


def test_get_residue_atoms_agrees_with_get_atom_pos():
    atoms = [_atom("N", 1.0, 0.4, altloc="A"),
             _atom("N", 8.0, 0.6, altloc="B"),
             _atom("CA", 2.0, 0.5, altloc="A"),
             _atom("CA", 7.0, 0.5, altloc="B")]
    residue = get_residue_atoms(atoms, "A", 1, "")
    single = get_atom_pos(atoms, "A", 1, "N", "")
    assert residue["N"][0] == 8.0 and single[0] == 8.0
    assert residue["CA"][0] == 2.0  # tie -> first


def test_extract_selects_best_conformer_per_name():
    class S:  # minimal structure stand-in
        atoms = [_atom("N", 1.0, 0.4, altloc="A", b=10.0),
                 _atom("N", 8.0, 0.6, altloc="B", b=50.0),
                 _atom("CA", 2.0, 0.5, altloc="A", b=20.0),
                 _atom("CA", 7.0, 0.5, altloc="B", b=60.0)]
        confidence = None
    residues = Tier1Measurements._extract(S)
    r = residues[0]
    assert set(r["_atoms"]) == {"N", "CA"}
    assert r["_atoms"]["N"][0] == 8.0      # higher occupancy
    assert r["_atoms"]["CA"][0] == 2.0     # tie -> first
    # confidence averages the RETAINED conformers' B-factors only
    assert r["_conf"] == (50.0 + 20.0) / 2


def _pdb_line(serial, name, altloc, x, occ, b, res="ALA", element=None):
    element = element or name[0]
    if name == "CA":
        element = "C"
    return (f"ATOM  {serial:5d} {name:>4s}{altloc:1s}{res:>3s} A   1 "
            f"  {x:8.3f}   0.000   0.000{occ:6.2f}{b:6.2f}"
            f"          {element:>2s}  ")


SYNTHETIC_ALTLOC = "\n".join([
    "HEADER    SYNTHETIC ALTLOC FIXTURE",
    "EXPDTA    X-RAY DIFFRACTION",
    _pdb_line(1, "N", "A", 0.000, 0.30, 10.00),
    _pdb_line(2, "N", "B", 9.000, 0.70, 40.00),
    _pdb_line(3, "CA", "B", 1.458, 0.70, 40.00),
    _pdb_line(4, "CA", "A", 1.500, 0.30, 10.00),
    "END",
]) + "\n"


def test_parser_dedups_altlocs_by_occupancy():
    structure = IngestionProcessor.run(
        SYNTHETIC_ALTLOC.encode("ascii"), "altloc.pdb", "Audit", "experimental")
    assert len(structure.atoms) == 2  # one N, one CA
    n = [a for a in structure.atoms if a.atom_name == "N"][0]
    ca = [a for a in structure.atoms if a.atom_name == "CA"][0]
    assert n.pos[0] == 9.0 and n.occupancy == 0.70   # B wins on occupancy
    assert ca.pos[0] == 1.458 and ca.occupancy == 0.70
    # confidence from retained conformers only: (40 + 40) / 2
    assert abs(structure.confidence.mean_plddt - 40.0) < 1e-9
