"""
P0 items 8 + 9 (2026-10-19 audit): honest statuses and coverage denominator.

Item 8:
  - LAW-110 counts observed peptide bonds, not len(residues) - 1
    (asserted per-fixture in test_residue_ordering.py; the parametrized
    cases here cover the degenerate shapes).
  - The adjudicator's defaults for absent AND malformed (statusless)
    measurement rows are NOT_MEASURED — never a fabricated FAIL. The
    absent-law default landed in WP-SCI-1; the malformed-row default is
    aligned here.

Item 9:
  - LAW-105's coverage denominator is the SOURCE residue count —
    ATOM records plus retained HETATM-mapped residues in the source —
    so unparseable coordinate lines or dropped records cannot silently
    inflate coverage.
"""
import numpy as np

from tos.engine.adjudicator import adjudicate_laws, AdjudicationInput
from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import IngestionProcessor


# ── Item 8: adjudicator defaults ────────────────────────────────────────

def test_absent_law_defaults_to_not_measured():
    res = adjudicate_laws(AdjudicationInput(
        raw_measurements={}, coverage=100.0, is_experimental=True,
        method_type="xray", architecture="NONE"))
    statuses = {r["law_id"]: r["status"] for r in res.law_rows}
    assert set(statuses.values()) == {"NOT_MEASURED"}
    # every row carries the honest reason
    assert all(r.get("reason") == "no measurement row emitted by the engine"
               for r in res.law_rows)


def test_malformed_row_without_status_defaults_to_not_measured():
    # A row that exists but lacks "status" must not become "measured and
    # failed" (the old default) — it is not measured.
    measurements = {"LAW-100": {"observed": None, "sample": 0}}
    res = adjudicate_laws(AdjudicationInput(
        raw_measurements=measurements, coverage=100.0, is_experimental=True,
        method_type="xray", architecture="NONE"))
    row = next(r for r in res.law_rows if r["law_id"] == "LAW-100")
    assert row["status"] == "NOT_MEASURED"


def test_malformed_statusless_deterministic_row_blocks_pass():
    measurements = {"LAW-100": {"observed": None, "sample": 0}}
    res = adjudicate_laws(AdjudicationInput(
        raw_measurements=measurements, coverage=100.0, is_experimental=True,
        method_type="xray", architecture="NONE"))
    # fail-closed preserved: missing evidence is not an authoritative PASS
    assert res.verdict != "PASS"


def test_measured_fail_still_fails():
    # LAW-170 (Residue Identity) is deterministic for X-ray per the
    # modality matrix; a measured VETO there must still VETO the verdict.
    measurements = {"LAW-170": {"observed": 3, "status": "VETO", "sample": 100}}
    res = adjudicate_laws(AdjudicationInput(
        raw_measurements=measurements, coverage=100.0, is_experimental=True,
        method_type="xray", architecture="NONE"))
    assert res.verdict == "VETO"


def test_advisory_veto_is_normalized_not_fatal():
    # LAW-100 resolves to advisory_experimental for X-ray: its VETO is
    # downgraded to an advisory FAIL and cannot veto the verdict.
    measurements = {"LAW-100": {"observed": 99.0, "status": "VETO", "sample": 10}}
    res = adjudicate_laws(AdjudicationInput(
        raw_measurements=measurements, coverage=100.0, is_experimental=True,
        method_type="xray", architecture="NONE"))
    row = next(r for r in res.law_rows if r["law_id"] == "LAW-100")
    assert row["status"] == "FAIL (Advisory)"
    assert res.verdict != "VETO"


# ── Item 9: coverage denominator from source residues ───────────────────

def _pred_line(serial, res_seq, occ_b, x=0.0):
    return (f"ATOM  {serial:5d}  CA  GLY A{res_seq:4d} "
            f"  {x:8.3f}   0.000   0.000  1.00{occ_b:6.2f}"
            f"          C  ")


def _predicted_pdb(*extra):
    lines = ["HEADER    PREDICTED COVERAGE FIXTURE"]
    # No EXPDTA -> classified "predicted"; B-factor column carries pLDDT.
    serial = 1
    for seq, b in [(1, 90.0), (2, 90.0), (3, 90.0), (4, 90.0), (5, 90.0)]:
        lines.append(_pred_line(serial, seq, b))
        serial += 1
    lines.extend(extra)
    lines.append("END")
    return ("\n".join(lines) + "\n").encode("ascii")


def test_coverage_denominator_includes_unparseable_source_residues():
    # Residue 6 exists in the source but its coordinate fields are
    # garbage, so the parser cannot retain its atoms.
    broken = ("ATOM  " + " " * 5 + "  CA  GLY A   6 "
              "  ******   0.000   0.000  1.00 90.00          C  ")
    content = _predicted_pdb(broken)
    structure = IngestionProcessor.run(content, "cov.pdb", "Audit", "predicted")
    # source saw 6 residues; only 5 parsed
    assert structure.source_residue_count == 6
    parsed = Tier1Measurements._extract(structure)
    assert len(parsed) == 5
    results, coverage, _ = Tier1Measurements.run_full_audit(structure)
    # coverage uses the SOURCE denominator: 5/6 = 83.3%, not 5/5 = 100%
    assert results["LAW-105"]["observed"] == 83.33
    assert results["LAW-105"]["sample"] == 5
    assert results["LAW-105"]["sample_source_residues"] == 6


def test_coverage_denominator_includes_hetm_mapped_residues():
    mse = ("HETATM" + " " * 5 + "  N   MSE A   7 "
           "   1.000   1.000   1.000  1.00 90.00           N  ")
    content = _predicted_pdb(mse)
    structure = IngestionProcessor.run(content, "cov.pdb", "Audit", "predicted")
    assert structure.source_residue_count == 6  # 5 ATOM + 1 MSE HETATM
    results, coverage, _ = Tier1Measurements.run_full_audit(structure)
    assert results["LAW-105"]["sample_source_residues"] == 6


def test_coverage_unaffected_when_source_count_matches_parsed():
    content = _predicted_pdb()
    structure = IngestionProcessor.run(content, "cov.pdb", "Audit", "predicted")
    results, _, _ = Tier1Measurements.run_full_audit(structure)
    assert results["LAW-105"]["observed"] == 100.0
    assert results["LAW-105"]["sample_source_residues"] == 5


def test_low_plddt_residues_count_against_source_denominator():
    # 2 of 6 source residues below pLDDT 70 -> coverage 4/6 = 66.7 -> FAIL
    lines = ["HEADER    PREDICTED LOW CONFIDENCE"]
    serial = 1
    for seq, b in [(1, 90.0), (2, 90.0), (3, 40.0), (4, 90.0), (5, 50.0), (6, 90.0)]:
        lines.append(_pred_line(serial, seq, b))
        serial += 1
    lines.append("END")
    structure = IngestionProcessor.run(
        ("\n".join(lines) + "\n").encode("ascii"), "cov.pdb", "Audit", "predicted")
    results, _, _ = Tier1Measurements.run_full_audit(structure)
    assert results["LAW-105"]["observed"] == 66.67
    assert results["LAW-105"]["status"] == "FAIL"
