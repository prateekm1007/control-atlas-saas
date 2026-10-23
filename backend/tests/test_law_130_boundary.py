"""
LAW-130 boundary and tooling-failure pins (OD-04 candidate).

The canonical comparison operator for LAW-130 is "<=" (station_sop
LAW_CANON). These tests pin:

  1. below threshold  (10.0)  -> PASS
  2. exact boundary   (20.0)  -> PASS  ( "<=" : the boundary complies)
  3. above threshold  (30.0)  -> VETO
  4. a neighbor-search tooling failure yields MEASUREMENT_FAILED (not a
     silent favorable zero), and blocks an authoritative PASS.
  5. the observation-precision contract: the REPORTED value is the raw
     score rounded to 2 decimals and the "<=" threshold is applied to that
     same reported value, so displayed value and status can never
     contradict. On measurement failure, sample is 0 and the input atom
     count is exposed separately (input_atom_count), never as a measured
     sample.
"""
import sys
import types

import numpy as np

from tos.engine.adjudicator import AdjudicationInput, adjudicate_laws
from tos.engine.tier1_measurements import Tier1Measurements
from tos.governance.station_sop import LAW_CANON
from tos.ingestion.processor import IngestionProcessor

# The measured set is always exactly 100 atoms (chain A + partner atoms),
# so each cross-chain clash contributes exactly 10.0 to the clash score.
TOTAL_ATOMS = 100


def _atom_line(serial, seq, chain, x, y, z):
    return (f"ATOM  {serial:5d} CA   GLY {chain}{seq:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.0:6.2f}{90.0:6.2f}          C  ")


def _clash_fixture(partner_chains):
    """Exactly 100 single-atom GLY CA atoms on a 4.0 A same-chain grid.

    `partner_chains` is a list of partner chain IDs, one per intended
    clash: each partner atom sits 1.2 A off a distinct chain-A residue
    (anchors 10, 60, 80 in order). The 4.0 A spacing keeps every
    same-chain sequential gap inside LAW-160's 4.5 A threshold, and the
    partner is far (>= 4.1 A) from every other chain-A atom, so the only
    sub-1.5 A pairs are the intended cross-chain ones. The engine's
    exclusion clause (same-chain |res_seq difference| <= 1) can never
    mask them, because partners are on different chains. Partner atoms
    are excluded from LAW-160's sequential scan, so chain integrity
    stays PASS in every case except the isolation-free fixture.
    """
    atoms = []
    n_partners = len(partner_chains)
    n_chain_a = TOTAL_ATOMS - n_partners
    serial = 1
    for seq in range(1, n_chain_a + 1):
        atoms.append(_atom_line(serial, seq, "A", seq * 4.0, 0.0, 0.0))
        serial += 1
    anchors = (10, 60, 80)
    for p_index, chain_id in enumerate(partner_chains):
        anchor = anchors[p_index]
        atoms.append(_atom_line(serial, 1, chain_id,
                                anchor * 4.0, 1.2, 0.0))
        serial += 1
    lines = ["HEADER    LAW-130 BOUNDARY FIXTURE"] + atoms + ["END"]
    return ("\n".join(lines) + "\n").encode("ascii")


def _audit(partner_chains):
    content = _clash_fixture(partner_chains)
    structure = IngestionProcessor.run(content, "law130.pdb", "Audit", "predicted")
    results, _, _ = Tier1Measurements.run_full_audit(structure)
    return results["LAW-130"]


def test_below_threshold_passes():
    res = _audit(["B"])  # clash score = 1/100 * 1000 = 10.0
    assert res["observed"] == 10.0
    assert res["status"] == "PASS"
    assert res["sample"] == TOTAL_ATOMS
    assert res["input_atom_count"] == TOTAL_ATOMS


def test_exact_boundary_complies_with_canonical_operator():
    res = _audit(["B", "C"])  # clash score = 2/100 * 1000 = 20.0
    assert res["observed"] == 20.0
    assert res["status"] == "PASS"  # 20.0 <= 20.0 under the canonical "<="
    assert float(LAW_CANON["LAW-130"]["threshold"]) == 20.0


def test_above_threshold_vetoes():
    res = _audit(["B", "C", "D"])  # 3/100 -> 30.0
    assert res["observed"] == 30.0
    assert res["status"] == "VETO"


class _ExplodingKDTree:
    def __init__(self, *args, **kwargs):
        raise RuntimeError("simulated neighbor-search failure")


def test_neighbor_search_failure_is_measurement_failed_not_pass(monkeypatch):
    """A tooling failure must surface as MEASUREMENT_FAILED. It must
    never degenerate into the old silent `except: pass` that scored a
    favorable zero, and must not allow an authoritative PASS."""
    fake_scipy = types.ModuleType("scipy")
    fake_spatial = types.ModuleType("scipy.spatial")
    fake_spatial.cKDTree = _ExplodingKDTree
    monkeypatch.setitem(sys.modules, "scipy", fake_scipy)
    monkeypatch.setitem(sys.modules, "scipy.spatial", fake_spatial)

    content = _clash_fixture(["B"])
    structure = IngestionProcessor.run(content, "law130.pdb", "Audit", "predicted")
    results, _, _ = Tier1Measurements.run_full_audit(structure)
    res = results["LAW-130"]

    assert res["status"] == "MEASUREMENT_FAILED"
    assert res["observed"] is None
    assert res["sample"] == 0
    # input atom count is exposed separately, never as a measured sample
    assert res["input_atom_count"] == TOTAL_ATOMS

    # fail-closed adjudication: missing evidence blocks a PASS
    verdict = adjudicate_laws(AdjudicationInput(
        raw_measurements=results,
        coverage=100.0,
        is_experimental=False,
        method_type="predicted",
        architecture="NONE",
    ))
    assert verdict.verdict == "INDETERMINATE"
    assert "LAW-130" in verdict.verdict_reason


def test_observation_precision_contract_status_equals_reported():
    """The threshold is applied to the REPORTED (rounded-to-2dp) value, so
    the displayed observed and the status always agree. We exercise the
    canonical boundary directly: 20.0 (exact) is compliant under "<="."""
    res = _audit(["B", "C"])
    assert res["observed"] == 20.0
    assert res["status"] == "PASS"
    # reported value is round(raw, 2); status is PASS iff reported <= threshold
    assert res["observed"] <= float(LAW_CANON["LAW-130"]["threshold"])
    assert (res["observed"] <= float(LAW_CANON["LAW-130"]["threshold"])) == (
        res["status"] == "PASS"
    )
