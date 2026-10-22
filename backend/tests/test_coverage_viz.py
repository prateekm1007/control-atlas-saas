import pytest
from tos.forensic_artifacts.confidence_viz import coverage_segments, generate_confidence_bar


def test_coverage_segments_are_measured_percentages():
    assert coverage_segments(62.5) == (62.5, 37.5)
    assert coverage_segments(-5) == (0.0, 100.0)
    assert coverage_segments(150) == (100.0, 0.0)


def test_confidence_bar_renders_without_plddt_bins():
    png = generate_confidence_bar(62.5, 10, 12)
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert len(png) > 1000
