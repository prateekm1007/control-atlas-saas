"""
P0 item 5 (2026-10-19 audit): HETATM retention tests.

The parser previously dropped every HETATM record. Selenomethionine (MSE)
is deposited as HETATM and is a normal methionine with Se for S; dropping
it silently deletes residues from the chain and fabricates backbone gaps.

Acceptance from the audit direction: "1A8O chain A yields 70 residues.
The ATOM-only parser yields 66." The real 1A8O bytes are committed at
tests/fixtures/real_pdb_1a8o.pdb (sha256
b0b4bf0d52c313047d5d7a7b188bf07a66a57c480db8ec4680f29461f18cbe4e).
"""
from pathlib import Path

from tos.ingestion.processor import IngestionProcessor, MODIFIED_RESIDUE_PARENTS

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1a8o.pdb"


def _parse(path=FIXTURE):
    content = path.read_bytes()
    return IngestionProcessor.run(content, "1a8o.pdb", "Audit", "experimental")


def _chain_residues(structure, chain="A"):
    seen = []
    for a in structure.atoms:
        if a.chain_id != chain:
            continue
        if (a.res_seq, a.insertion_code) not in seen:
            seen.append((a.res_seq, a.insertion_code))
    return seen


# ── Acceptance: real 1A8O ───────────────────────────────────────────────

def test_1a8o_chain_a_yields_70_residues():
    structure = _parse()
    residues = _chain_residues(structure, "A")
    assert len(residues) == 70
    # 66 from ATOM records + 4 MSE HETATM residues (151, 185, 214, 215)
    assert (151, "") in residues and (215, "") in residues


def test_1a8o_mse_residues_present_and_mapped_to_met():
    structure = _parse()
    met_positions = {a.res_seq for a in structure.atoms
                     if a.chain_id == "A" and a.res_name == "MET"
                     and a.atom_name == "CA"}
    # The four MSE residues (151, 185, 214, 215) now count as MET.
    assert {151, 185, 214, 215} <= met_positions


def test_1a8o_waters_excluded():
    structure = _parse()
    assert all(a.res_name != "HOH" for a in structure.atoms)
    # 88 HOH oxygen atoms in the file; none may survive parsing.
    assert sum(1 for a in structure.atoms if a.res_name == "HOH") == 0


def test_1a8o_mse_selenium_atoms_retained():
    # The SE atoms are real heavy atoms of the retained residues.
    se = [a for a in structure_atoms(_parse()) if a.element == "SE"]
    assert len(se) == 4


def structure_atoms(structure):
    return structure.atoms


# ── Unit behavior on a synthetic file ───────────────────────────────────

SYNTHETIC = """\
HEADER    SYNTHETIC HETATM FIXTURE
EXPDTA    X-RAY DIFFRACTION
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 20.00           C
ATOM      4  N   ALA A   2       3.200   1.600   0.100  1.00 20.00           N
ATOM      5  CA  ALA A   2       4.000   2.800   0.200  1.00 20.00           C
HETATM    6  N   MSE A   3       5.100   2.700   0.300  1.00 22.00           N
HETATM    7  CA  MSE A   3       6.400   3.100   0.400  1.00 22.00           C
HETATM    8  SE  MSE A   3       7.500   4.200   0.500  1.00 24.00          SE
HETATM    9  O   HOH A 100      20.000  20.000  20.000  1.00 30.00           O
HETATM   10  ZN  ZN  A 101      21.000  21.000  21.000  1.00 30.00          ZN
END
"""


def test_synthetic_mse_retained_mapped_water_and_zinc_dropped():
    structure = IngestionProcessor.run(
        SYNTHETIC.encode("ascii"), "synthetic.pdb", "Audit", "experimental")
    residues = _chain_residues(structure, "A")
    assert residues == [(1, ""), (2, ""), (3, "")]
    names = {a.res_name for a in structure.atoms}
    assert names == {"ALA", "MET"}          # MSE -> MET; HOH and ZN gone
    assert any(a.element == "SE" for a in structure.atoms)


def test_modified_residue_map_contains_no_d_amino_acids():
    # D-amino acids must stay unmapped so LAW-145 can flag them.
    d_residues = {k for k in MODIFIED_RESIDUE_PARENTS if k.startswith("D")}
    assert d_residues == set()


def test_map_covers_the_canonical_modified_residues():
    for expected in ("MSE", "SEC", "PYL", "HYP"):
        assert expected in MODIFIED_RESIDUE_PARENTS
