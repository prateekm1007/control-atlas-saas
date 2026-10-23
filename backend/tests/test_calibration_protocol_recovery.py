"""Byte-level pin for the recovered calibration evaluation protocol."""
import hashlib
from pathlib import Path

PROTOCOL = (
    Path(__file__).resolve().parents[1]
    / "calibration" / "protocol" / "evaluation_protocol_v1.md"
)
PROVENANCE = PROTOCOL.with_name("PROVENANCE.md")


def test_recovered_protocol_is_byte_exact():
    data = PROTOCOL.read_bytes()
    assert len(data) == 8364
    assert hashlib.sha256(data).hexdigest() == (
        "c880086f04f59b577a34ead13eafa2de6cc5ae337059bc16b1ed87434a724806"
    )


def test_protocol_is_document_of_record_not_live_gate():
    provenance = " ".join(PROVENANCE.read_text(encoding="utf-8").split())
    assert "document of record only" in provenance
    assert "not a live gate" in provenance
    assert "OD-16" in provenance
    assert "calibration_status" in provenance
    assert "not_approved" in provenance
