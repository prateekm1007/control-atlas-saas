import gemmi
import numpy as np
import pytest

from tos.engine.identity import compute_coordinate_hash
from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import IngestionError, IngestionProcessor


PDB_WITH_HETATM = """\
HEADER    TEST
EXPDTA    X-RAY DIFFRACTION
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA AALA A   1       1.458   0.000   0.000  0.40 90.00           C
ATOM      3  CA BALA A   1       1.500   0.000   0.000  0.60 90.00           C
ATOM      4  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      5  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
ATOM      6  N   GLY A   2       3.300   1.424   0.000  1.00 90.00           N
ATOM      7  CA  GLY A   2       4.758   1.424   0.000  1.00 90.00           C
ATOM      8  C   GLY A   2       5.304   2.848   0.000  1.00 90.00           C
ATOM      9  O   GLY A   2       5.745   3.997   0.000  1.00 90.00           O
HETATM   10  N   MSE A   3       6.600   2.848   0.000  1.00 90.00           N
HETATM   11  CA  MSE A   3       8.058   2.848   0.000  1.00 90.00           C
HETATM   12  C   MSE A   3       8.604   4.272   0.000  1.00 90.00           C
HETATM   13  O   MSE A   3       9.045   5.421   0.000  1.00 90.00           O
HETATM   14  SE  MSE A   3       8.900   1.500   0.000  1.00 90.00          SE
HETATM   15  O   HOH A 101      20.000  20.000  20.000  1.00 50.00           O
HETATM   16  C1  LIG A 201      30.000  30.000  30.000  1.00 50.00           C
END
"""


def _as_mmcif(pdb_text: str) -> bytes:
    structure = gemmi.read_pdb_string(pdb_text)
    return structure.make_mmcif_document().as_string().encode()


def _by_atom(structure):
    return {(a.chain_id, a.res_seq, a.insertion_code, a.res_name, a.atom_name): a for a in structure.atoms}


class TestCanonicalIngestion:
    def test_pdb_and_mmcif_are_measurement_equivalent(self):
        pdb = IngestionProcessor.run(PDB_WITH_HETATM.encode(), "input.pdb", "Audit", "Audit")
        cif = IngestionProcessor.run(_as_mmcif(PDB_WITH_HETATM), "input.cif", "Audit", "Audit")
        assert pdb.metadata.format == "pdb"
        assert cif.metadata.format == "mmcif"
        assert compute_coordinate_hash(pdb) == compute_coordinate_hash(cif)
        assert [(a.chain_id, a.res_seq, a.atom_name) for a in pdb.atoms] == [
            (a.chain_id, a.res_seq, a.atom_name) for a in cif.atoms
        ]
        pdb_results, pdb_coverage, _ = Tier1Measurements.run_full_audit(pdb)
        cif_results, cif_coverage, _ = Tier1Measurements.run_full_audit(cif)
        assert pdb_coverage == cif_coverage
        for law_id in ("LAW-100", "LAW-105", "LAW-110", "LAW-120", "LAW-125", "LAW-135", "LAW-160"):
            assert pdb_results[law_id]["status"] == cif_results[law_id]["status"]
            assert pdb_results[law_id]["observed"] == cif_results[law_id]["observed"]
            assert pdb_results[law_id]["sample"] == cif_results[law_id]["sample"]

    def test_hetatm_polymer_is_included_and_water_ligand_are_excluded(self):
        structure = IngestionProcessor.run(PDB_WITH_HETATM.encode(), "input.pdb", "Audit", "Audit")
        residues = {(a.chain_id, a.res_seq, a.res_name) for a in structure.atoms}
        assert ("A", 3, "MSE") in residues
        assert ("A", 101, "HOH") not in residues
        assert ("A", 201, "LIG") not in residues
        assert structure.metadata.excluded_waters == 1
        assert structure.metadata.excluded_ligands == 1
        assert any(a.res_seq == 3 and a.is_hetatm for a in structure.atoms)

    def test_altloc_selection_uses_occupancy_and_reports_conflicts(self):
        structure = IngestionProcessor.run(PDB_WITH_HETATM.encode(), "input.pdb", "Audit", "Audit")
        ca = next(a for a in structure.atoms if a.res_seq == 1 and a.atom_name == "CA")
        assert ca.altloc == "B"
        assert ca.occupancy == pytest.approx(0.6)
        assert structure.metadata.altloc_conflicts == 1
        assert structure.metadata.altloc_selected_atoms == 1

    def test_altloc_tie_break_prefers_blank_then_lexicographic(self):
        pdb = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA AALA A   1       1.458   0.000   0.000  0.50 90.00           C
ATOM      3  CA BALA A   1       1.500   0.000   0.000  0.50 90.00           C
ATOM      4  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      5  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""
        structure = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        ca = next(a for a in structure.atoms if a.atom_name == "CA")
        assert ca.altloc == "A"

    def test_first_model_is_selected_and_other_models_are_counted(self):
        pdb = """\
MODEL        1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
ENDMDL
MODEL        2
ATOM      5  N   ALA A   1       9.000   0.000   0.000  1.00 90.00           N
ATOM      6  CA  ALA A   1      10.458   0.000   0.000  1.00 90.00           C
ATOM      7  C   ALA A   1      11.004   1.424   0.000  1.00 90.00           C
ATOM      8  O   ALA A   1      11.445   2.573   0.000  1.00 90.00           O
ENDMDL
END
"""
        structure = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        assert structure.metadata.excluded_models == 1
        assert all(a.pos[0] < 9.0 for a in structure.atoms)

    def test_experimental_coverage_is_source_backed_not_hardcoded(self):
        pdb = """\
HEADER    TEST
EXPDTA    X-RAY DIFFRACTION
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 20.00           O
ATOM      5  CA  GLY A   2       5.000   0.000   0.000  1.00 20.00           C
END
"""
        structure = IngestionProcessor.run(pdb.encode(), "input.pdb", "Audit", "Audit")
        results, coverage, is_experimental = Tier1Measurements.run_full_audit(structure)
        assert is_experimental is True
        assert coverage == 50.0
        assert results["LAW-105"]["coverage_basis"] == "complete_backbone_residues"
        assert results["LAW-105"]["coverage_numerator"] == 1
        assert results["LAW-105"]["coverage_denominator"] == 2
        assert results["LAW-105"]["source_atom_records"] == 5

    def test_unrecognized_source_fails_closed(self):
        with pytest.raises(IngestionError):
            IngestionProcessor.run(b"not a structure", "input.txt", "Audit", "Audit")
