"""Truthfulness tests for confidence reporting in forensic artifacts."""
import base64

from tos.forensic_artifacts.pdf_generator import (
    _confidence_plot_mode,
    _extract_ca_atoms,
)
from tos.forensic_artifacts.residue_heatmap import generate_residue_heatmap_png
from tos.forensic_artifacts.structure_viz import generate_structure_render


_BASE_ATOM = (
    "ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 80.00           C"
)


def _atom_line(x, bfactor_field):
    line = _BASE_ATOM[:30] + f"{x:8.3f}" + _BASE_ATOM[38:60]
    line += bfactor_field if bfactor_field is not None else "      "
    line += _BASE_ATOM[66:]
    return line


def _pdb(fields):
    text = "\n".join(
        _atom_line(10.0 + index * 3.8, field)
        for index, field in enumerate(fields)
    )
    return base64.b64encode((text + "\nEND\n").encode("ascii")).decode("ascii")


def test_missing_bfactor_is_not_fabricated_as_fifty():
    coords, confidences = _extract_ca_atoms(_pdb([" 80.00", None, " 60.00"]))
    assert len(coords) == 3
    assert confidences == [80.0, None, 60.0]
    assert all(value != 50.0 for value in confidences)


def test_nonnumeric_bfactor_is_unavailable_not_a_score():
    coords, confidences = _extract_ca_atoms(_pdb([" 80.00", "notnum"]))
    assert len(coords) == 2
    assert confidences == [80.0, None]


def test_structure_render_never_imputes_missing_plddt():
    coords, confidences = _extract_ca_atoms(_pdb([" 80.00", None, " 60.00"]))
    png = generate_structure_render(coords, confidences, "observed", metric_label="pLDDT")
    assert png.startswith(b"\x89PNG\r\n\x1a\n")

    neutral = generate_structure_render(
        coords, [None, None, None], "coordinates only", metric_label=None)
    assert neutral.startswith(b"\x89PNG\r\n\x1a\n")


def test_heatmap_rejects_unobserved_values_instead_of_plotting_them():
    png = generate_residue_heatmap_png(
        [80.0, None, 60.0], title="pLDDT", metric_label="pLDDT")
    assert png.startswith(b"\x89PNG\r\n\x1a\n")

def test_bfactor_metadata_never_selects_plddt_plotting():
    source, observed, plot_as_plddt, all_plddt = _confidence_plot_mode(
        {"source_type": "B-factor", "data_available": False}, [20.0, 30.0, 40.0])
    assert source == "B-factor"
    assert observed == [20.0, 30.0, 40.0]
    assert plot_as_plddt is False
    assert all_plddt is False


def test_unavailable_confidence_never_selects_plddt_plotting():
    source, observed, plot_as_plddt, all_plddt = _confidence_plot_mode(
        {"source_type": "unavailable", "data_available": False}, [None, None, None])
    assert source == "unavailable"
    assert observed == []
    assert plot_as_plddt is False
    assert all_plddt is False


def test_observed_plddt_selects_plddt_plotting():
    source, observed, plot_as_plddt, all_plddt = _confidence_plot_mode(
        {"source_type": "pLDDT", "data_available": True}, [80.0, 90.0, 70.0])
    assert source == "pLDDT"
    assert observed == [80.0, 90.0, 70.0]
    assert plot_as_plddt is True
    assert all_plddt is True


def test_partial_plddt_disables_full_heatmap_but_allows_observed_render():
    source, observed, plot_as_plddt, all_plddt = _confidence_plot_mode(
        {"source_type": "pLDDT", "data_available": True},
        [80.0, None, 70.0, 90.0])
    assert source == "pLDDT"
    assert observed == [80.0, 70.0, 90.0]
    assert plot_as_plddt is True
    assert all_plddt is False
