import pytest
from tos.witness.contract import WitnessContractEnforcer, WitnessViolation
from tos.engine.tier3_predict import StrategicPredictor
def test_witness_cage():
    e = WitnessContractEnforcer()
    with pytest.raises(WitnessViolation): e.validate_output("No disclaimer")
def test_math_bounds():
    r = StrategicPredictor.calculate(100, 95, 0.8)
    assert 5 <= r["probability"] <= 95
