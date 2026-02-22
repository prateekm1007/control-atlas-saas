"""
Tests for batch processing engine.
Tests the processor directly (no HTTP) and validates schema compliance.
"""
import io
import json
import zipfile
import pytest
from tos.engine.batch_processor import (
    process_batch,
    _extract_pdb_files,
    _derive_candidate_id,
    BatchResult,
)
from tos.schemas.batch_v1 import BatchResponse, BatchStructureResult, BatchSummary


# ── Helpers ───────────────────────────────────────────────────────────────────

MINIMAL_PDB = (
    "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\n"
    "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C\n"
    "ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 10.00           C\n"
    "ATOM      4  O   ALA A   1       4.000   5.000   6.000  1.00 10.00           O\n"
    "ATOM      5  CB  ALA A   1       2.500   2.500   5.000  1.00 10.00           C\n"
    "ATOM      6  N   ALA A   2       4.000   5.000   6.000  1.00 10.00           N\n"
    "ATOM      7  CA  ALA A   2       5.000   6.000   7.000  1.00 10.00           C\n"
    "ATOM      8  C   ALA A   2       6.000   7.000   8.000  1.00 10.00           C\n"
    "ATOM      9  O   ALA A   2       7.000   8.000   9.000  1.00 10.00           O\n"
    "ATOM     10  CB  ALA A   2       5.500   5.500   8.000  1.00 10.00           C\n"
    "END\n"
)


def _make_zip(files: dict[str, bytes]) -> bytes:
    """Create a ZIP archive in memory from a dict of filename→content."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, content in files.items():
            zf.writestr(name, content)
    return buf.getvalue()


def _mock_physics(content_bytes, candidate_id, mode, t3_category):
    """Mock physics engine that returns a minimal valid payload."""
    return {
        "verdict": {
            "binary": "PASS",
            "deterministic_score": 80,
            "advisory_score": 100,
            "physical_score": 80,
            "confidence_score": 90.0,
            "det_passed": 10,
            "det_total": 12,
            "heur_passed": 3,
            "heur_total": 3,
            "coverage_pct": 100.0,
            "suppression_reason": None,
        },
        "governance": {
            "audit_id": "TEST1234",
            "station_version": "22.5.3",
            "timestamp_utc": "2026-01-01T00:00:00Z",
            "governance_fingerprint": {
                "canon_hash": "6a9cd4b4349b81de",
                "matrix_hash": "ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4",
                "matrix_schema_version": "1.0.0",
                "policy_ref": "PIL-CAL-03",
            },
        },
        "provenance": {"source": candidate_id, "hash": "abc123", "byte_count": len(content_bytes)},
        "tier1": {"laws": []},
        "tier3": {"probability": 80},
        "characterization": {
            "total_atoms": 10, "total_residues": 2,
            "source_type": "predicted", "resolution": None, "method": "Standard",
        },
        "confidence_meta": {
            "source_type": "pLDDT", "provenance_method": "static",
            "mean": 90.0, "data_available": True,
        },
        "strategic_math": {"s6": 0.8, "w_arch": 1.0, "m_s8": 0.0, "architecture": "NONE"},
        "witness_reports": None,
        "ai_model_used": "mock",
        "pdf_b64": "",
        "pdb_b64": "",
    }


def _mock_physics_veto(content_bytes, candidate_id, mode, t3_category):
    """Mock physics engine that returns VETO."""
    result = _mock_physics(content_bytes, candidate_id, mode, t3_category)
    result["verdict"]["binary"] = "VETO"
    result["verdict"]["deterministic_score"] = 40
    result["tier1"]["laws"] = [
        {
            "law_id": "LAW-120", "title": "Bond Angle", "status": "VETO",
            "method": "deterministic", "observed": 25.0, "threshold": 18.0,
            "operator": "<=", "units": "°", "deviation": "7.0",
            "sample_size": 100, "scope": "core residues", "principle": "N/A",
        }
    ]
    return result


# ── ZIP Extraction Tests ─────────────────────────────────────────────────────

class TestZipExtraction:

    def test_extracts_pdb_files(self):
        zip_bytes = _make_zip({"test.pdb": MINIMAL_PDB.encode()})
        files = _extract_pdb_files(zip_bytes)
        assert len(files) == 1
        assert files[0][0] == "test.pdb"

    def test_skips_non_pdb_files(self):
        zip_bytes = _make_zip({
            "structure.pdb": MINIMAL_PDB.encode(),
            "readme.txt": b"ignore me",
            "data.csv": b"1,2,3",
        })
        files = _extract_pdb_files(zip_bytes)
        assert len(files) == 1

    def test_handles_nested_directories(self):
        zip_bytes = _make_zip({
            "structures/4HHB.pdb": MINIMAL_PDB.encode(),
            "structures/1CRN.pdb": MINIMAL_PDB.encode(),
        })
        files = _extract_pdb_files(zip_bytes)
        assert len(files) == 2

    def test_skips_macosx_artifacts(self):
        zip_bytes = _make_zip({
            "test.pdb": MINIMAL_PDB.encode(),
            "__MACOSX/._test.pdb": b"artifact",
        })
        files = _extract_pdb_files(zip_bytes)
        assert len(files) == 1

    def test_rejects_invalid_zip(self):
        with pytest.raises(ValueError, match="Invalid ZIP"):
            _extract_pdb_files(b"not a zip file")

    def test_rejects_empty_zip(self):
        zip_bytes = _make_zip({"readme.txt": b"no pdbs here"})
        with pytest.raises(ValueError, match="no valid PDB"):
            _extract_pdb_files(zip_bytes)

    def test_accepts_ent_files(self):
        zip_bytes = _make_zip({"pdb4hhb.ent": MINIMAL_PDB.encode()})
        files = _extract_pdb_files(zip_bytes)
        assert len(files) == 1


class TestCandidateIdDerivation:

    def test_simple_filename(self):
        assert _derive_candidate_id("4HHB.pdb") == "4HHB"

    def test_nested_path(self):
        assert _derive_candidate_id("structures/proteins/1CRN.pdb") == "1CRN"

    def test_ent_extension(self):
        assert _derive_candidate_id("pdb4hhb.ent") == "pdb4hhb"


# ── Batch Processing Tests ───────────────────────────────────────────────────

class TestBatchProcessing:

    def test_three_structures_all_pass(self):
        zip_bytes = _make_zip({
            "a.pdb": MINIMAL_PDB.encode(),
            "b.pdb": MINIMAL_PDB.encode(),
            "c.pdb": MINIMAL_PDB.encode(),
        })
        result = process_batch(zip_bytes, _mock_physics)
        assert result.summary.total == 3
        assert result.summary.passed == 3
        assert result.summary.vetoed == 0
        assert result.summary.errors == 0
        assert all(r.success for r in result.results)

    def test_mixed_pass_and_veto(self):
        zip_bytes = _make_zip({
            "good.pdb": MINIMAL_PDB.encode(),
            "bad.pdb": MINIMAL_PDB.encode(),
        })

        call_count = [0]
        def alternating_physics(content, cid, mode, cat):
            call_count[0] += 1
            if call_count[0] % 2 == 0:
                return _mock_physics_veto(content, cid, mode, cat)
            return _mock_physics(content, cid, mode, cat)

        result = process_batch(zip_bytes, alternating_physics)
        assert result.summary.total == 2
        assert result.summary.passed == 1
        assert result.summary.vetoed == 1

    def test_corrupt_file_handled_gracefully(self):
        zip_bytes = _make_zip({
            "good.pdb": MINIMAL_PDB.encode(),
            "corrupt.pdb": MINIMAL_PDB.encode(),
        })

        call_count = [0]
        def failing_physics(content, cid, mode, cat):
            call_count[0] += 1
            if call_count[0] == 2:
                raise RuntimeError("Parse error")
            return _mock_physics(content, cid, mode, cat)

        result = process_batch(zip_bytes, failing_physics)
        assert result.summary.total == 2
        assert result.summary.passed == 1
        assert result.summary.errors == 1
        assert not result.results[1].success
        assert "Parse error" in result.results[1].error

    def test_governance_fingerprint_per_structure(self):
        zip_bytes = _make_zip({
            "a.pdb": MINIMAL_PDB.encode(),
            "b.pdb": MINIMAL_PDB.encode(),
        })
        result = process_batch(zip_bytes, _mock_physics)
        for r in result.results:
            fp = r.payload["governance"]["governance_fingerprint"]
            assert fp["canon_hash"] == "6a9cd4b4349b81de"
            assert len(fp["matrix_hash"]) == 64

    def test_max_structures_enforced(self):
        files = {f"s{i}.pdb": MINIMAL_PDB.encode() for i in range(10)}
        zip_bytes = _make_zip(files)
        with pytest.raises(ValueError, match="exceeding limit"):
            process_batch(zip_bytes, _mock_physics, max_structures=5)

    def test_common_failing_laws_tracked(self):
        zip_bytes = _make_zip({
            "a.pdb": MINIMAL_PDB.encode(),
            "b.pdb": MINIMAL_PDB.encode(),
        })
        result = process_batch(zip_bytes, _mock_physics_veto)
        assert result.summary.vetoed == 2
        assert len(result.summary.common_failing_laws) > 0
        assert result.summary.common_failing_laws[0]["law_id"] == "LAW-120"
        assert result.summary.common_failing_laws[0]["count"] == 2

    def test_mean_score_calculated(self):
        zip_bytes = _make_zip({
            "a.pdb": MINIMAL_PDB.encode(),
            "b.pdb": MINIMAL_PDB.encode(),
        })
        result = process_batch(zip_bytes, _mock_physics)
        assert result.summary.mean_deterministic_score == 80.0


# ── Schema Compliance Tests ──────────────────────────────────────────────────

class TestBatchSchemaCompliance:

    def test_result_converts_to_pydantic(self):
        zip_bytes = _make_zip({"test.pdb": MINIMAL_PDB.encode()})
        result = process_batch(zip_bytes, _mock_physics)

        # Build the Pydantic response the same way main.py does
        response = BatchResponse(
            results=[
                BatchStructureResult(
                    filename=r.filename,
                    candidate_id=r.candidate_id,
                    success=r.success,
                    response=r.payload if r.success else None,
                    error=r.error,
                )
                for r in result.results
            ],
            summary=BatchSummary(
                total=result.summary.total,
                passed=result.summary.passed,
                vetoed=result.summary.vetoed,
                indeterminate=result.summary.indeterminate,
                errors=result.summary.errors,
                mean_deterministic_score=result.summary.mean_deterministic_score,
                common_failing_laws=result.summary.common_failing_laws,
            ),
        )

        assert response.summary.total == 1
        assert response.results[0].success is True
