
"""
Regression tests for LAW-130 (Clashscore, neighbor exclusion).

Round-26 auditor direction: the measurement had two defects - a swallowed
calculation failure that left a zero-valued PASS, and a strict `<` comparison
against a canon that declares the inclusive operator `<=` at 20.0.

The KDTree is injected (scipy is optional in this environment and the
measurement must be deterministic): a fake `scipy.spatial.cKDTree` returns a
controlled pair set, or raises, so the boundary, the veto and the failure path
are all exercised exactly.
"""
import os
import sys
import types

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.ingestion.processor import Atom, Structure, ConfidenceProfile
from tos.engine.tier1_measurements import Tier1Measurements


def _install_fake_scipy(monkeypatch, pairs=(), raise_error=False):
    """A fake scipy.spatial whose cKDTree returns `pairs` (or raises)."""
    module = types.ModuleType("scipy")
    spatial = types.ModuleType("scipy.spatial")

    class cKDTree:  # noqa: N801 - mirrors the scipy name
        def __init__(self, coords):
            if raise_error:
                raise RuntimeError("KDTree construction failed")
            self._coords = coords

        def query_pairs(self, r):  # noqa: ARG002 - the radius is fixed by the law
            return list(pairs)

    spatial.cKDTree = cKDTree
    module.spatial = spatial
    monkeypatch.setitem(sys.modules, "scipy", module)
    monkeypatch.setitem(sys.modules, "scipy.spatial", spatial)


def _atom(x, y, z, res_seq, chain="A", name="ALA", atom="CA", element="C"):
    return Atom(atom_name=atom, element=element, pos=np.array([x, y, z], dtype=float),
                res_name=name, res_seq=res_seq, chain_id=chain, insertion_code="", b_iso=90.0)


def _structure(atoms):
    return Structure(atoms=atoms, audit_id="TEST",
                     confidence=ConfidenceProfile([], False, "predicted"))


def _run(atoms):
    results, _, _ = Tier1Measurements.run_full_audit(_structure(atoms))
    return results["LAW-130"]


def _chain(n_atoms, spacing=10.0, chain="A"):
    """n_atoms spread far apart in one chain (no accidental clashes)."""
    return [_atom(i * spacing, 0.0, 0.0, res_seq=i + 1, chain=chain) for i in range(n_atoms)]


class TestLaw130:
    def test_exactly_the_threshold_passes(self, monkeypatch):
        """The canon declares operator `<=` at 20.0: exactly 20.0 must PASS."""
        atoms = _chain(50)
        # one clash in 50 atoms => (1/50)*1000 = 20.0 exactly
        _install_fake_scipy(monkeypatch, pairs=[(0, 1)])
        res = _run(atoms)
        assert res["observed"] == pytest.approx(20.0)
        assert res["status"] == "PASS", "the declared operator is inclusive (<=)"

    def test_above_the_threshold_vetoes(self, monkeypatch):
        """A clashscore above 20.0 must VETO."""
        atoms = _chain(50)
        _install_fake_scipy(monkeypatch, pairs=[(0, 1), (2, 3)])
        res = _run(atoms)
        assert res["observed"] == pytest.approx(40.0)
        assert res["status"] == "VETO"

    def test_a_calculation_failure_is_not_measured_not_a_zero_pass(self, monkeypatch):
        """A failed calculation must never become a zero-valued PASS."""
        atoms = _chain(10)
        _install_fake_scipy(monkeypatch, raise_error=True)
        res = _run(atoms)
        assert res["status"] == "NOT_MEASURED"
        assert res.get("measured") is False
        assert res["observed"] is None, "a failed calculation must not report a number"
        assert res.get("deviation") is None
        assert res["sample"] == 0

    def test_clash_counts_and_sample_semantics(self, monkeypatch):
        """`sample` is the number of atoms the clashscore is normalised over."""
        atoms = _chain(50)
        _install_fake_scipy(monkeypatch, pairs=[(0, 1)])
        res = _run(atoms)
        assert res["sample"] == 50
        assert res["observed"] == pytest.approx(20.0)

    def test_residue_pair_reporting_is_deterministic(self, monkeypatch):
        """The failing residue pairs are reported in a deterministic order."""
        atoms = (_chain(4, chain="A") + _chain(3, chain="B", spacing=10.0))
        _install_fake_scipy(monkeypatch, pairs=[(4, 5), (0, 1)])
        res = _run(atoms)
        pairs = res["failing_residue_pairs"]
        assert pairs == [("A:1", "A:2"), ("B:1", "B:2")]
        # the same input yields the same output, and the pairs are deduplicated
        _install_fake_scipy(monkeypatch, pairs=[(0, 1), (1, 0)])
        assert _run(atoms)["failing_residue_pairs"] == pairs

    def test_a_single_atom_has_no_pair_to_measure(self, monkeypatch):
        """One atom cannot yield a clashscore: the metric is absent, not zero."""
        _install_fake_scipy(monkeypatch)
        res = _run(_chain(1))
        assert res["status"] == "NOT_MEASURED"
        assert res["observed"] is None
        assert res["sample"] == 0
