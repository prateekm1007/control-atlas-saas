"""Authorized OD-03/OD-04 scope tests.

These tests pin the canonical semantics:
  - pLDDT is advisory metadata only for LAW-130 and LAW-200;
  - LAW-105 remains residue-coverage and total-residue scoped;
  - LAW-110 remains all-residue in original order;
  - LAW-155 is retired from verdict and scoring semantics;
  - missing mandatory measurements cannot become PASS.
"""
import numpy as np
import pytest

from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws
from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Atom, ConfidenceProfile, Structure
from tos.governance.station_sop import (
    HEURISTIC_COUNT,
    LAW_CANON,
    RETIRED_LAWS,
)


def _atom(name, pos, seq, b_iso, chain="A", res_name="ALA", icode=""):
    return Atom(
        atom_name=name,
        element=name[0],
        pos=np.array(pos, dtype=float),
        res_name=res_name,
        res_seq=seq,
        chain_id=chain,
        insertion_code=icode,
        b_iso=b_iso,
    )


def _structure(atoms, experimental=False):
    from tos.ingestion.processor import IngestionMetadata
    chains = {}
    for a in atoms:
        chains.setdefault(a.chain_id, set()).add((a.res_seq, a.insertion_code))
    md = IngestionMetadata(
        source_polymer_residues_by_chain={k: len(v) for k, v in chains.items()},
        source_sequence_available=True,
    )
    return Structure(
        atoms=atoms,
        audit_id="TEST",
        confidence=ConfidenceProfile([], experimental, "xray" if experimental else "predicted"),
        metadata=md,
    )


def _residue(chain, seq, c_pos, n_pos, b_iso=90.0, icode=""):
    return [
        _atom("N", n_pos, seq, b_iso, chain=chain, icode=icode),
        _atom("CA", [n_pos[0] + 1.0, n_pos[1], n_pos[2]], seq, b_iso, chain=chain, icode=icode),
        _atom("C", c_pos, seq, b_iso, chain=chain, icode=icode),
    ]


def _all_passing():
    return {
        lid: {
            "observed": 0.0,
            "deviation": "0.0",
            "sample": 100,
            "status": "PASS",
        }
        for lid in LAW_CANON
    }


def _clash_structure(extra_clash=False):
    atoms = []
    for seq in range(1, 26):
        x = seq * 20.0
        atoms.append(_atom("N", [x, 0, 0], seq, 90.0))
        atoms.append(_atom("CA", [x + 5.0, 0, 0], seq, 90.0))

    atoms[0] = _atom("N", [0.0, 0, 0], 1, 90.0)
    atoms[4] = _atom("N", [1.0, 0, 0], 3, 90.0)
    if extra_clash:
        atoms[6] = _atom("N", [80.0, 0, 0], 4, 90.0)
        atoms[10] = _atom("N", [81.0, 0, 0], 6, 90.0)
    return _structure(atoms)


class TestLaw130:
    def test_exactly_at_threshold_passes_with_inclusive_operator(self):
        struct = _clash_structure()
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        row = results["LAW-130"]
        assert row["observed"] == 20.0
        assert row["status"] == "PASS"
        assert row["sample"] == 50

    def test_above_threshold_vetoes(self):
        struct = _clash_structure(extra_clash=True)
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        row = results["LAW-130"]
        assert row["observed"] == 40.0
        assert row["status"] == "VETO"

    def test_below_threshold_passes(self):
        struct = _structure([
            _atom("N", [0, 0, 0], 1, 90.0),
            _atom("CA", [20, 0, 0], 1, 90.0),
            _atom("N", [40, 0, 0], 3, 90.0),
            _atom("CA", [60, 0, 0], 3, 90.0),
        ])
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        row = results["LAW-130"]
        assert row["observed"] == 0.0
        assert row["status"] == "PASS"

    def test_neighbor_search_failure_is_explicit(self, monkeypatch):
        import scipy.spatial

        def boom(*args, **kwargs):
            raise RuntimeError("synthetic tree failure")

        monkeypatch.setattr(scipy.spatial, "cKDTree", boom)
        results, _, _ = Tier1Measurements.run_full_audit(_clash_structure())
        row = results["LAW-130"]
        assert row["status"] == "MEASUREMENT_FAILED"
        assert row["observed"] is None
        assert "neighbor search failed" in row["reason"]

    def test_empty_core_does_not_suppress_all_atom_measurement(self):
        atoms = _residue("A", 1, [0, 0, 0], [0, 5, 0], b_iso=10.0)
        results, _, _ = Tier1Measurements.run_full_audit(_structure(atoms))
        row = results["LAW-130"]
        assert row["sample"] == 3
        assert row["status"] in ("PASS", "VETO")
        assert row["plddt_advisory"]["fraction_at_or_above_70"] == 0.0
        assert row["plddt_advisory"]["policy"].startswith("advisory only")

    def test_plddt_is_advisory_for_law_130(self):
        high = _residue("A", 1, [0, 0, 0], [0, 5, 0], b_iso=90.0)
        low = _residue("A", 1, [0, 0, 0], [0, 5, 0], b_iso=10.0)
        high_results, _, _ = Tier1Measurements.run_full_audit(_structure(high))
        low_results, _, _ = Tier1Measurements.run_full_audit(_structure(low))
        assert high_results["LAW-130"]["status"] == low_results["LAW-130"]["status"]
        assert high_results["LAW-130"]["observed"] == low_results["LAW-130"]["observed"]
        assert high_results["LAW-130"]["plddt_advisory"] != low_results["LAW-130"]["plddt_advisory"]


class TestLaw105Coverage:
    def test_returns_residue_coverage_not_atom_coverage(self):
        atoms = []
        for seq in range(1, 13):
            atoms.extend(_residue("A", seq, [seq * 20, 0, 0], [seq * 20, 5, 0], b_iso=90.0))
        for seq in range(13, 17):
            for i in range(12):
                atoms.append(_atom("CB", [seq * 20, i * 2, 0], seq, 10.0))
        results, coverage, _ = Tier1Measurements.run_full_audit(_structure(atoms))
        assert coverage == 75.0
        assert results["LAW-105"]["observed"] == 75.0
        assert "atom_coverage" not in results["LAW-105"]

    def test_residue_gate_is_authoritative_over_atom_fraction(self):
        # 75% of residues but only 50% of atoms are above the pLDDT boundary.
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=_all_passing(),
            coverage=75.0,
            is_experimental=False,
            method_type="predicted",
            architecture="NONE",
        ))
        assert adjudication.verdict == "PASS"

    def test_below_residue_gate_is_indeterminate_with_reason(self):
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=_all_passing(),
            coverage=69.9,
            is_experimental=False,
            method_type="predicted",
            architecture="NONE",
        ))
        assert adjudication.verdict == "INDETERMINATE"
        assert "LAW-105 residue coverage" in adjudication.verdict_reason


class TestLaw110Adjacency:
    def test_chain_boundary_is_not_paired(self):
        atoms = _residue("A", 1, [0, 0, 0], [0, 5, 0])
        atoms += _residue("B", 2, [50, 0, 0], [1, 0, 0])
        results, _, _ = Tier1Measurements.run_full_audit(_structure(atoms))
        assert results["LAW-110"]["status"] == "NOT_MEASURED"
        assert results["LAW-110"]["observed"] is None
        assert results["LAW-110"]["sample"] == 0

    def test_residue_gap_is_not_paired(self):
        atoms = _residue("A", 1, [0, 0, 0], [0, 5, 0])
        atoms += _residue("A", 3, [50, 0, 0], [1, 0, 0])
        results, _, _ = Tier1Measurements.run_full_audit(_structure(atoms))
        assert results["LAW-110"]["status"] == "NOT_MEASURED"
        assert results["LAW-110"]["observed"] is None
        assert results["LAW-110"]["sample"] == 0

    def test_insertion_codes_remain_sequential(self):
        atoms = _residue("A", 1, [0, 0, 0], [0, 5, 0], icode="")
        atoms += _residue("A", 1, [10, 0, 0], [5, 0, 0], icode="A")
        atoms += _residue("A", 1, [20, 0, 0], [15, 0, 0], icode="B")
        results, _, _ = Tier1Measurements.run_full_audit(_structure(atoms))
        assert results["LAW-110"]["observed"] == 2

    def test_non_core_residues_are_still_audited(self):
        atoms = _residue("A", 1, [0, 0, 0], [0, 5, 0], b_iso=10.0)
        atoms += _residue("A", 2, [50, 0, 0], [5, 0, 0], b_iso=10.0)
        results, _, _ = Tier1Measurements.run_full_audit(_structure(atoms))
        assert results["LAW-110"]["observed"] == 1
        assert results["LAW-110"]["sample"] == 1


class TestRetiredAndAdvisorySemantics:
    def test_law_155_is_advisory_and_cannot_veto(self):
        measurements = _all_passing()
        measurements["LAW-155"] = {
            "observed": 0.1, "deviation": "999", "sample": 10, "status": "VETO",
        }
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=measurements,
            coverage=100.0,
            is_experimental=False,
            method_type="predicted",
            architecture="NONE",
        ))
        row = next(r for r in adjudication.law_rows if r["law_id"] == "LAW-155")
        assert row["method"] == "heuristic"
        assert row["status"] == "FAIL (Advisory)"
        assert row["raw_status"] == "VETO"
        assert row["enforced"] is True
        assert adjudication.verdict == "PASS"
        # advisory_score counts PASS heuristic laws; LAW-155 adverse here is a
        # FAIL (Advisory), so score reflects the two passing heuristics.
        assert adjudication.advisory_score == 66

    def test_retired_law_is_traceable_when_missing(self):
        measurements = _all_passing()
        measurements.pop("LAW-155")
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=measurements,
            coverage=100.0,
            is_experimental=False,
            method_type="predicted",
            architecture="NONE",
        ))
        row = next(r for r in adjudication.law_rows if r["law_id"] == "LAW-155")
        assert row["status"] == "NOT_MEASURED"
        assert row["raw_status"] == "NOT_MEASURED"
        assert adjudication.verdict == "PASS"

    @pytest.mark.parametrize("law_id", ["LAW-182", "LAW-200"])
    def test_heuristic_veto_is_normalized_not_final(self, law_id):
        measurements = _all_passing()
        measurements[law_id] = {
            "observed": 0.0, "deviation": "999", "sample": 10, "status": "VETO",
        }
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=measurements,
            coverage=100.0,
            is_experimental=False,
            method_type="predicted",
            architecture="NONE",
        ))
        row = next(r for r in adjudication.law_rows if r["law_id"] == law_id)
        assert row["raw_status"] == "VETO"
        assert row["status"] == "FAIL (Advisory)"
        assert adjudication.verdict == "PASS"

    def test_law_130_deterministic_veto_still_vetoes(self):
        measurements = _all_passing()
        measurements["LAW-130"] = {
            "observed": 99.0, "deviation": "999", "sample": 10, "status": "VETO",
        }
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=measurements,
            coverage=100.0,
            is_experimental=False,
            method_type="predicted",
            architecture="NONE",
        ))
        row = next(r for r in adjudication.law_rows if r["law_id"] == "LAW-130")
        assert row["method"] == "deterministic"
        assert row["status"] == "VETO"
        assert adjudication.verdict == "VETO"

    def test_counts_keep_all_ids_but_retire_only_frozen(self):
        assert len(LAW_CANON) == 15
        assert RETIRED_LAWS == set()
        assert HEURISTIC_COUNT == 3  # LAW-155, LAW-182, LAW-200


class TestMissingMeasurements:
    def test_missing_deterministic_measurement_is_not_measured(self):
        measurements = _all_passing()
        measurements.pop("LAW-130")
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=measurements,
            coverage=100.0,
            is_experimental=False,
            method_type="predicted",
            architecture="NONE",
        ))
        row = next(r for r in adjudication.law_rows if r["law_id"] == "LAW-130")
        assert row["status"] == "NOT_MEASURED"
        assert row["observed"] is None
        assert adjudication.verdict == "INDETERMINATE"
        assert "LAW-130" in adjudication.verdict_reason

    def test_failed_deterministic_measurement_is_not_pass(self):
        measurements = _all_passing()
        measurements["LAW-130"] = {
            "observed": None, "deviation": None, "sample": 0,
            "status": "MEASUREMENT_FAILED", "reason": "synthetic failure",
        }
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=measurements,
            coverage=100.0,
            is_experimental=False,
            method_type="predicted",
            architecture="NONE",
        ))
        assert adjudication.verdict == "INDETERMINATE"
        assert "LAW-130" in adjudication.verdict_reason


class TestLaw100CanonicalMetric:
    def _bond_atoms(self, cb_distance=1.63):
        atoms = [
            _atom("N", [0.0, 0.0, 0.0], 1, 90.0),
            _atom("CA", [1.458, 0.0, 0.0], 1, 90.0),
            _atom("C", [2.004, 1.424, 0.0], 1, 90.0),
            _atom("O", [2.445, 2.573, 0.0], 1, 90.0),
            _atom("CB", [1.458, 0.0, cb_distance], 1, 90.0),
        ]
        return atoms

    def test_experimental_uses_percentage_not_unregistered_rmsz(self):
        struct = _structure(self._bond_atoms(), experimental=True)
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        row = results["LAW-100"]
        assert row["sample"] == 4
        assert row["observed"] == 25.0
        assert row["status"] == "VETO"

    def test_law_100_threshold_is_canonical_percentage(self):
        assert LAW_CANON["LAW-100"]["threshold"] == 5.0
        assert LAW_CANON["LAW-100"]["operator"] == "<="
        assert LAW_CANON["LAW-100"]["unit"] == "%"
        assert LAW_CANON["LAW-100"]["type"] == "Percentage"

    def test_no_bonds_is_not_measured(self):
        struct = _structure([_atom("CA", [0.0, 0.0, 0.0], 1, 90.0)])
        results, _, _ = Tier1Measurements.run_full_audit(struct)
        row = results["LAW-100"]
        assert row["status"] == "NOT_MEASURED"
        assert row["observed"] is None
        assert row["sample"] == 0


class TestConfidenceContract:
    def test_uses_profile_observations_and_source_metadata(self):
        struct = Structure(
            atoms=[_atom("CA", [0, 0, 0], 1, 10.0)],
            audit_id="CONF",
            confidence=ConfidenceProfile([71.5, 84.5], False, "predicted"),
        )
        source, mean, provenance = Tier1Measurements.detect_confidence_source(struct)
        assert source == "pLDDT"
        assert mean == 78.0
        assert provenance == "predicted"
        meta = Tier1Measurements.decompose_confidence(struct)
        assert meta["schema_version"] == "confidence-meta/1.0.0"
        assert meta["data_available"] is True
        assert meta["observation_count"] == 2

    def test_unavailable_confidence_is_explicit_not_hardcoded(self):
        struct = Structure(
            atoms=[], audit_id="EMPTY",
            confidence=ConfidenceProfile([], False, "predicted"),
        )
        source, mean, provenance = Tier1Measurements.detect_confidence_source(struct)
        assert source == "unavailable"
        assert mean is None
        assert provenance == "predicted"
        meta = Tier1Measurements.decompose_confidence(struct)
        assert meta["mean"] is None
        assert meta["data_available"] is False
        assert meta["observation_count"] == 0

    def test_experimental_b_factors_are_not_plddt_or_generic_confidence(self):
        struct = Structure(
            atoms=[_atom("CA", [0, 0, 0], 1, 20.0)],
            audit_id="BF",
            confidence=ConfidenceProfile([20.0], True, "xray"),
        )
        source, mean, provenance = Tier1Measurements.detect_confidence_source(struct)
        assert source == "B-factor"
        assert mean is None
        assert provenance == "xray"
        meta = Tier1Measurements.decompose_confidence(struct)
        assert meta["mean"] is None
        assert meta["data_available"] is False
        assert meta["b_factor_mean"] == 20.0
        assert meta["b_factor_unit"] == "B-factor"
        assert meta["b_factor_observation_count"] == 1

    def test_predicted_plddt_mean_is_reported_without_b_factor_confusion(self):
        struct = Structure(
            atoms=[_atom("CA", [0, 0, 0], 1, 20.0)],
            audit_id="PLDDT",
            confidence=ConfidenceProfile([80.0], False, "predicted"),
        )
        source, mean, provenance = Tier1Measurements.detect_confidence_source(struct)
        assert source == "pLDDT"
        assert mean == 80.0
        assert provenance == "predicted"
        meta = Tier1Measurements.decompose_confidence(struct)
        assert meta["mean"] == 80.0
        assert meta["data_available"] is True
        assert meta["observation_count"] == 1
        assert meta.get("b_factor_mean") is None

class TestNonFiniteCoordinates:
    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
    def test_malformed_coordinates_cannot_pass(self, bad):
        atoms = _residue("A", 1, [0, 0, 0], [0, 5, 0])
        atoms += _residue("A", 2, [bad, 0, 0], [5, 0, 0])
        results, coverage, is_experimental = Tier1Measurements.run_full_audit(
            _structure(atoms, experimental=True))
        assert coverage is None
        assert is_experimental is True
        assert results["LAW-105"]["status"] == "MEASUREMENT_FAILED"
        assert "non-finite" in results["LAW-105"]["reason"]
        adjudication = adjudicate_laws(AdjudicationInput(
            raw_measurements=results,
            coverage=coverage,
            is_experimental=True,
            method_type="xray",
            architecture="NONE",
        ))
        assert adjudication.verdict == "INDETERMINATE"
