import pytest
from tos.witness.contract import WitnessContractEnforcer, WitnessViolation
def test_witness_cage():
    enforcer = WitnessContractEnforcer()
    with pytest.raises(WitnessViolation):
        enforcer.validate_output("This drug will work definitely.")
