"""Ingestion contract tests (P6): HETATM polymer classification, altloc
accounting, explicit failure modes, and PDB/mmCIF residue parity."""
import gemmi
import pytest

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import IngestionError, IngestionProcessor
from tos.engine.identity import compute_coordinate_hash


def test_hetatm_polymer_retained_water_and_ligand_excluded():
    pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
HETATM    5  N   MSE A   2       3.000   0.000   0.000  1.00 90.00           N
HETATM    6  CA  MSE A   2       4.458   0.000   0.000  1.00 90.00           C
HETATM    7  C   MSE A   2       5.004   1.424   0.000  1.00 90.00           C
HETATM    8  O   MSE A   2       5.445   2.573   0.000  1.00 90.00           O
HETATM    9  SE  MSE A   2       6.000   0.000   0.000  1.00 90.00           SE
HETATM   10  O   HOH A   3       9.000   0.000   0.000  1.00 90.00           O
HETATM   11  C1  LIG A   4      10.000   0.000   0.000  1.00 90.00           C
END
"""
    s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    assert s.metadata.excluded_waters == 1
    assert s.metadata.excluded_ligands == 1
    # MSE (polymer-like HETATM) is retained as a residue and measured.
    residues = {r["_name"] for r in Tier1Measurements._extract(s)}
    assert "MSE" in residues
    assert "ALA" in residues
    assert "HOH" not in residues
    assert "LIG" not in residues


def test_altloc_dominant_conformer_and_discard_accounting():
    pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA AALA A   1       1.458   0.000   0.000  0.40 90.00           C
ATOM      3  CA BALA A   1       1.500   0.000   0.000  0.60 90.00           C
ATOM      4  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      5  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
    s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    ca = next(a for a in s.atoms if a.atom_name == "CA")
    assert ca.altloc == "B"            # dominant conformer (0.60 > 0.40)
    # Two CA altlocs -> one selected; no artificial mixed residue.
    assert s.metadata.altloc_selected_atoms >= 1
    assert s.metadata.selected_atoms <= s.metadata.source_atom_records


def test_incomplete_residue_is_retained_not_fabricated_as_complete():
    pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
ATOM      5  N   GLY A   2       3.000   0.000   0.000  1.00 90.00           N
ATOM      6  C   GLY A   2       5.000   0.000   0.000  1.00 90.00           C
ATOM      7  O   GLY A   2       5.500   0.000   0.000  1.00 90.00           O
END
"""
    s = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    assert s.metadata.complete_backbone_residues == 1   # only ALA
    assert s.metadata.total_residues == 2               # GLY retained (incomplete)
    results, _, _ = Tier1Measurements.run_full_audit(s)
    assert results["LAW-105"]["complete_backbone_residues"] == 1


@pytest.mark.parametrize("bad,ext", [
    (b"", "input.pdb"),
    (b"not a structure at all", "input.txt"),
    (b'{"not": "a structure"}', "input.json"),
])
def test_malformed_or_empty_input_raises(bad, ext):
    with pytest.raises(IngestionError):
        IngestionProcessor.run(bad, ext, "Audit", "Audit")


def test_pdb_mmcif_parity_residue_identity():
    pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
ATOM      5  N   GLY A   2       3.000   0.000   0.000  1.00 90.00           N
ATOM      6  CA  GLY A   2       4.000   0.000   0.000  1.00 90.00           C
ATOM      7  C   GLY A   2       5.000   0.000   0.000  1.00 90.00           C
ATOM      8  O   GLY A   2       5.500   0.000   0.000  1.00 90.00           O
END
"""
    s_pdb = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
    cif = gemmi.read_pdb_string(pdb).make_mmcif_document().as_string().encode()
    s_cif = IngestionProcessor.run(cif, "input.cif", "Audit", "Audit")
    assert compute_coordinate_hash(s_pdb) == compute_coordinate_hash(s_cif)
    r1 = {r["_seq"]: r["_name"] for r in Tier1Measurements._extract(s_pdb)}
    r2 = {r["_seq"]: r["_name"] for r in Tier1Measurements._extract(s_cif)}
    assert r1 == r2
