"""
Unit tests for LAW-160 (Chain Integrity).
"""
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from tos.engine.tier1_measurements import Tier1Measurements
from tos.ingestion.processor import Structure, ConfidenceProfile


def test_law_160_pass(simple_backbone):
    """simple_backbone has ~3.8-3.95Å CA-CA spacing — should PASS."""
    struct = Structure(atoms=simple_backbone, audit_id="TEST",
                       confidence=ConfidenceProfile([], False, "predicted"))
    results, coverage, _ = Tier1Measurements.run_full_audit(struct)

    res = results["LAW-160"]
    assert res["status"] == "PASS"
    # Fixture geometry gives CA-CA between 3.7 and 4.0
    assert 3.7 <= res["observed"] <= 4.0
    # Sample should be CA pairs, not total residues
    assert res["sample"] == 2


def test_law_160_fail(broken_backbone):
    """broken_backbone has a ~10Å gap — should VETO."""
    struct = Structure(atoms=broken_backbone, audit_id="TEST",
                       confidence=ConfidenceProfile([], False, "predicted"))
    results, coverage, _ = Tier1Measurements.run_full_audit(struct)

    res = results["LAW-160"]
    assert res["status"] == "VETO"
    assert res["observed"] > 9.0
