"""Pytest conftest for the acceptance suite (Round-19 fix, 2026-10-20).

Soundness fix: the suite's check() helper ACCUMULATES failures into module-level
PASS/FAIL lists (that is what the direct runner in run_all.py reports). Plain
pytest never looked at those lists, so a test function whose subprocess checks
all failed still PASSED as far as pytest was concerned -- an unsound green.
This autouse fixture runs its assertion in the teardown of EVERY test, so any
check() failure recorded during that test fails the owning pytest node with the
recorded messages. Accumulate-and-report semantics are preserved (all checks in
a group still execute); the verdict is no longer decoupled from pytest's exit
status. This applies to every pytest invocation of these node IDs, including
verify_claimed_shas.py --check-tests (EC-20) and the live CI acceptance-suite
job.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # noqa: E402


@pytest.fixture(autouse=True)
def _fail_test_when_check_recorded_failures(request):
    import test_acceptance as _ta
    before = len(_ta.FAIL)
    yield
    new = _ta.FAIL[before:]
    assert not new, (
        f"{request.node.name}: {len(new)} check() assertion(s) FAILED: " + " | ".join(new))
