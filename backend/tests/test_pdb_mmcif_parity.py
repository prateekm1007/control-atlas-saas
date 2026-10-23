"""P6 regression tests: PDB/mmCIF parity through the Gemmi ingestion path.

The integration candidate uses a Gemmi-based parser for both PDB and mmCIF.
A parser that accepts both formats is not, by itself, evidence of numerical
parity. These tests ingest the same structure in both formats and assert the
engine's residue identity, coverage, and per-law numerical outputs agree.
"""
import os
from pathlib import Path
import numpy as np
import pytest

os.environ.setdefault("TOSCANINI_API_KEYS", "test-key-xyz")
from tos.ingestion.processor import IngestionProcessor
from tos.engine.tier1_measurements import Tier1Measurements

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_1crn.pdb"


def _pdb_bytes():
    return FIXTURE.read_bytes()


def _mmcif_bytes(pdb_bytes):
    import gemmi
    st = gemmi.read_pdb_string(pdb_bytes.decode(errors="replace"))
    return st.make_mmcif_document().as_string().encode("utf-8")


def _audit(content, filename, mode="experimental"):
    structure = IngestionProcessor.run(content, filename, "Audit", mode)
    results, coverage, is_exp = Tier1Measurements.run_full_audit(structure)
    return structure, results, coverage, is_exp


class TestPDBMMcifParity:

    def test_ingestion_accepts_both_formats(self):
        pdb = _pdb_bytes()
        cif = _mmcif_bytes(pdb)
        assert IngestionProcessor.run(pdb, "1crn.pdb", "Audit", "experimental") is not None
        assert IngestionProcessor.run(cif, "1crn.cif", "Audit", "experimental") is not None

    def test_residue_identity_parity(self):
        pdb = _pdb_bytes()
        cif = _mmcif_bytes(pdb)
        s_pdb, _, _, _ = _audit(pdb, "1crn.pdb")
        s_cif, _, _, _ = _audit(cif, "1crn.cif")
        # residue identities (chain, seq, icode, name) must match
        def residue_ids(s):
            return sorted({(a.chain_id, a.res_seq, a.insertion_code,
                           a.res_name.strip().upper()) for a in s.atoms})
        assert residue_ids(s_pdb) == residue_ids(s_cif)

    def test_per_law_output_parity(self):
        pdb = _pdb_bytes()
        cif = _mmcif_bytes(pdb)
        _, r_pdb, cov_pdb, _ = _audit(pdb, "1crn.pdb")
        _, r_cif, cov_cif, _ = _audit(cif, "1crn.cif")
        assert abs(cov_pdb - cov_cif) < 1e-6
        for lid in set(r_pdb) & set(r_cif):
            a, b = r_pdb[lid], r_cif[lid]
            assert a.get("status") == b.get("status"), (lid, a.get("status"), b.get("status"))
            if a.get("observed") is not None and b.get("observed") is not None:
                assert abs(float(a["observed"]) - float(b["observed"])) < 1e-6, lid
