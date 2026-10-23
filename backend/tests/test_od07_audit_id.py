"""OD-07 audit-identity binding verification (P7).

The audit ID must be a deterministic full SHA-256 over the canonical serialization
of (coordinate_hash, canon_hash, matrix_hash): identical inputs yield identical
IDs, and changing any one required component changes the ID.
"""
import pytest

from tos.engine.identity import compute_audit_id, compute_coordinate_hash
from tos.governance.station_sop import LAW_CANON_HASH
from tos.governance.modality_matrix import compute_matrix_hash
from tos.ingestion.processor import IngestionProcessor

PDB_BASE = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.004   1.424   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       2.445   2.573   0.000  1.00 90.00           O
END
"""


def _audit():
    s = IngestionProcessor.run(PDB_BASE.encode(), "a.pdb", "Audit", "Audit")
    return s, compute_coordinate_hash(s)


def test_audit_id_is_deterministic_full_sha256():
    _, h = _audit()
    a1 = compute_audit_id(h, LAW_CANON_HASH, compute_matrix_hash())
    a2 = compute_audit_id(h, LAW_CANON_HASH, compute_matrix_hash())
    assert a1 == a2
    assert len(a1) == 64
    assert all(c in "0123456789abcdef" for c in a1)


def test_each_component_change_changes_audit_id():
    _, h = _audit()
    base = compute_audit_id(h, LAW_CANON_HASH, compute_matrix_hash())
    assert base != compute_audit_id(h, "0" * 16, compute_matrix_hash())
    assert base != compute_audit_id(h, LAW_CANON_HASH, "0" * 64)
    assert base != compute_audit_id("0" * 64, LAW_CANON_HASH, compute_matrix_hash())


def test_coordinate_change_changes_audit_id():
    s1, h1 = _audit()
    pdb2 = PDB_BASE.replace("1.458", "1.500")
    s2 = IngestionProcessor.run(pdb2.encode(), "b.pdb", "Audit", "Audit")
    h2 = compute_coordinate_hash(s2)
    assert h1 != h2
    assert compute_audit_id(h1, LAW_CANON_HASH, compute_matrix_hash()) != \
        compute_audit_id(h2, LAW_CANON_HASH, compute_matrix_hash())
