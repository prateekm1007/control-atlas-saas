"""
P0 items 11-14 (2026-10-19 audit): AlphaFold acquisition repair.

Covers the metadata-first acquisition module (tos/generation/afdb.py),
the dispatcher rewiring, gemmi-based mmCIF ingestion, and the
fail-closed cohort cases:

  - selection by modelEntityId, never data[0]
  - highest latestVersion among matches
  - cifUrl download; format from CONTENT; JSON refused as a structure
  - URL / version / SHA-256 / field-names recorded on the artifact
  - P0DTC2: migrated identifier recorded, fail closed
  - P62988: UniProt 'Inactive' status helper

No test touches the network: the transport is injected.
"""
import json

import pytest

from tos.generation.afdb import (
    AFDBUnavailable,
    AcquiredPrediction,
    detect_structure_format,
    fetch_prediction,
    normalize_accession,
    uniprot_entry_type,
)
from tos.generation.dispatcher import GenerationDispatcher, SourceUnavailable
from tos.ingestion.processor import IngestionProcessor

CIF_BYTES = b"""data_AF-Q5VSL9-F1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 1 ? 0.000 0.000 0.000 1.0 90.0 1 A 1
ATOM 2 C CA . GLY A 1 2 ? 1.500 0.000 0.000 1.0 90.0 1 A 1
HETATM 3 SE SE . MSE B 2 3 ? 5.000 5.000 5.000 1.0 88.0 3 B 1
HETATM 4 O O . HOH W 3 . ? 9.000 9.000 9.000 1.0 30.0 5 W 1
#
"""

API_ENTRY_V6 = {
    "modelEntityId": "AF-Q5VSL9-F1",
    "latestVersion": 6,
    "allVersions": [1, 2, 3, 4, 5, 6],
    "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v6.cif",
    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v6.pdb",
    "uniprotAccession": "Q5VSL9",
    "globalMetricValue": 80.88,
}
API_ENTRY_V4 = {**API_ENTRY_V6, "latestVersion": 4,
                "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v4.cif"}
API_ENTRY_OTHER = {**API_ENTRY_V6, "modelEntityId": "AF-Q5VSL9-F2",
                   "cifUrl": "https://example/other.cif"}


class FakeResponse:
    def __init__(self, status_code=200, content=b""):
        self.status_code = status_code
        self.content = content


class FakeTransport:
    """url -> FakeResponse; records requested urls."""

    def __init__(self, mapping):
        self.mapping = mapping
        self.requested = []

    def __call__(self, url, timeout=None):
        self.requested.append(url)
        resp = self.mapping.get(url)
        if resp is None:
            return FakeResponse(404, b"{}")
        return resp


def _api_url(acc):
    return f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"


# ── Format detection from content ───────────────────────────────────────

def test_format_detection_cif_pdb_and_json_refused():
    assert detect_structure_format(CIF_BYTES) == "cif"
    assert detect_structure_format(b"HEADER    X\nATOM      1  N   ALA A   1\n") == "pdb"
    with pytest.raises(ValueError):
        detect_structure_format(json.dumps([API_ENTRY_V6]).encode())
    with pytest.raises(ValueError):
        detect_structure_format(b"")
    with pytest.raises(ValueError):
        detect_structure_format(b"\x00\x01random garbage not a structure")


# ── Metadata-first selection ────────────────────────────────────────────

def test_selects_by_entity_id_never_data0():
    # data[0] is a DIFFERENT entity (F2); the match is later in the list.
    transport = FakeTransport({
        _api_url("Q5VSL9"): FakeResponse(200, json.dumps([API_ENTRY_OTHER, API_ENTRY_V6]).encode()),
        API_ENTRY_V6["cifUrl"]: FakeResponse(200, CIF_BYTES),
    })
    pred = fetch_prediction("AF-Q5VSL9-F1", http_get=transport)
    assert pred.model_entity_id == "AF-Q5VSL9-F1"
    assert pred.fmt == "cif"


def test_selects_highest_latest_version():
    # data[0] is the SAME entity at v4; the v6 entry is second.
    transport = FakeTransport({
        _api_url("Q5VSL9"): FakeResponse(200, json.dumps([API_ENTRY_V4, API_ENTRY_V6]).encode()),
        API_ENTRY_V6["cifUrl"]: FakeResponse(200, CIF_BYTES),
    })
    pred = fetch_prediction("Q5VSL9", http_get=transport)
    assert pred.latest_version == 6
    assert pred.source_url == API_ENTRY_V6["cifUrl"]
    assert "model_v4" not in pred.source_url


def test_provenance_recorded_url_version_sha256_fields():
    transport = FakeTransport({
        _api_url("Q5VSL9"): FakeResponse(200, json.dumps([API_ENTRY_V6]).encode()),
        API_ENTRY_V6["cifUrl"]: FakeResponse(200, CIF_BYTES),
    })
    pred = fetch_prediction("Q5VSL9-F1", http_get=transport)
    import hashlib
    assert pred.sha256 == hashlib.sha256(CIF_BYTES).hexdigest()
    assert pred.uniprot_accession == "Q5VSL9"
    assert pred.global_metric_value == 80.88
    assert "modelEntityId" in pred.fields_present
    assert "cifUrl" in pred.fields_present


# ── Fail-closed paths ───────────────────────────────────────────────────

def test_p0dtc2_migrated_identifier_recorded_then_fail_closed():
    # Live-verified 2026-10-19: /api/prediction/P0DTC2 answers with
    # modelEntityId AF-0000000365840314 (provider VR3D, ColabFold), NOT
    # AF-P0DTC2-F1. The requested entity must fail closed with the
    # migrated identifier recorded.
    migrated = {**API_ENTRY_V6,
                "modelEntityId": "AF-0000000365840314",
                "uniprotAccession": "P0DTC2",
                "latestVersion": 1}
    transport = FakeTransport({
        _api_url("P0DTC2"): FakeResponse(200, json.dumps([migrated]).encode()),
    })
    with pytest.raises(AFDBUnavailable) as excinfo:
        fetch_prediction("AF-P0DTC2-F1", http_get=transport)
    assert excinfo.value.record["migrated_identifier"] == "AF-0000000365840314"
    assert "AF-0000000365840314" in str(excinfo.value)


def test_p62988_no_prediction_fails_closed():
    # Live-verified 2026-10-19: the API returns 404/empty for P62988;
    # UniProt says entryType 'Inactive' (merged entry).
    transport = FakeTransport({})
    with pytest.raises(AFDBUnavailable):
        fetch_prediction("AF-P62988-F1", http_get=transport)


def test_uniprot_entry_type_helper():
    transport = FakeTransport({
        "https://rest.uniprot.org/uniprotkb/P62988.json":
            FakeResponse(200, json.dumps({"entryType": "Inactive"}).encode()),
    })
    assert uniprot_entry_type("P62988", http_get=transport) == "Inactive"
    assert uniprot_entry_type("P00000", http_get=FakeTransport({})) is None


def test_json_payload_refused_as_structure():
    # The API entry lacks cifUrl-ish behavior in the OLD code returned
    # the JSON body itself as a 'structure'. Now: fail closed.
    entry_no_cif = {k: v for k, v in API_ENTRY_V6.items() if k != "cifUrl"}
    transport = FakeTransport({
        _api_url("Q5VSL9"): FakeResponse(200, json.dumps([entry_no_cif]).encode()),
    })
    with pytest.raises(AFDBUnavailable) as excinfo:
        fetch_prediction("Q5VSL9", http_get=transport)
    assert "cifUrl" in str(excinfo.value)

    # And a cifUrl whose download returns JSON is refused too (body is
    # padded past the 200-byte minimum so the FORMAT check is what fires).
    big_json = json.dumps({"detail": "not a structure, just a large JSON body "
                                     "that would previously have been passed "
                                     "to the engine as a .pdb file",
                           "padding": "x" * 200}).encode()
    transport2 = FakeTransport({
        _api_url("Q5VSL9"): FakeResponse(200, json.dumps([API_ENTRY_V6]).encode()),
        API_ENTRY_V6["cifUrl"]: FakeResponse(200, big_json),
    })
    with pytest.raises(AFDBUnavailable) as excinfo2:
        fetch_prediction("Q5VSL9", http_get=transport2)
    assert "JSON" in str(excinfo2.value)


def test_api_failures_fail_closed():
    assert_ok = FakeTransport({_api_url("Q5VSL9"): FakeResponse(500, b"")})
    with pytest.raises(AFDBUnavailable):
        fetch_prediction("Q5VSL9", http_get=assert_ok)
    empty = FakeTransport({_api_url("Q5VSL9"): FakeResponse(200, b"[]")})
    with pytest.raises(AFDBUnavailable):
        fetch_prediction("Q5VSL9", http_get=empty)
    garbage = FakeTransport({_api_url("Q5VSL9"): FakeResponse(200, b"not json")})
    with pytest.raises(AFDBUnavailable):
        fetch_prediction("Q5VSL9", http_get=garbage)


def test_normalize_accession():
    assert normalize_accession("AF-P01308-F1") == "P01308"
    assert normalize_accession("AF-Q5VSL9-F2") == "Q5VSL9"
    assert normalize_accession("p01308") == "P01308"


# ── Dispatcher wiring ───────────────────────────────────────────────────

def test_dispatcher_afdb_path_uses_metadata_first(monkeypatch):
    transport = FakeTransport({
        _api_url("Q5VSL9"): FakeResponse(200, json.dumps([API_ENTRY_V6]).encode()),
        API_ENTRY_V6["cifUrl"]: FakeResponse(200, CIF_BYTES),
    })
    monkeypatch.setattr("tos.generation.afdb._http_get", transport)
    content, cid, fmt = GenerationDispatcher.acquire("AF-Q5VSL9-F1", None)
    assert content == CIF_BYTES
    assert fmt == "cif"
    assert not any("model_v4" in u for u in transport.requested)


def test_dispatcher_fails_closed_through_source_unavailable(monkeypatch):
    monkeypatch.setattr("tos.generation.afdb._http_get", FakeTransport({}))
    with pytest.raises(SourceUnavailable):
        GenerationDispatcher.acquire("AF-P62988-F1", None)


# ── gemmi mmCIF ingestion (item 13) ─────────────────────────────────────

def test_cif_ingestion_semantics():
    structure = IngestionProcessor.run(CIF_BYTES, "x.cif", "Audit", "predicted")
    # ATOM GLY kept; MSE HETATM kept and mapped to MET; HOH dropped
    names = {(a.res_name, a.chain_id, a.res_seq) for a in structure.atoms}
    assert ("GLY", "A", 1) in names
    assert ("MET", "B", 3) in names
    assert all(a.res_name != "HOH" for a in structure.atoms)
    assert all(a.element != "SE" or a.res_name == "MET" for a in structure.atoms)
    assert any(a.element == "SE" for a in structure.atoms)
    # predicted classification (no _exptl in AFDB-style file)
    assert structure.confidence.method == "predicted"
    assert structure.confidence.is_experimental is False
    # source residues: GLY A 1, MSE B 3 (HOH W does not count)
    assert structure.source_residue_count == 2
    # label_seq_id carried for ordering
    assert all(a.label_seq_id is not None for a in structure.atoms
               if a.res_name == "GLY")


def test_cif_experimental_classification():
    experimental_cif = (b"data_X\n_exptl.method 'X-RAY DIFFRACTION'\n"
                        b"_refine.ls_d_res_high 1.5\n" + CIF_BYTES.split(b"\n", 1)[1])
    structure = IngestionProcessor.run(experimental_cif, "x.cif", "Audit", "experimental")
    assert structure.confidence.method == "xray"
    assert structure.confidence.is_experimental is True
    assert structure.confidence.resolution == 1.5


def test_cif_label_seq_orders_residues():
    # Same residue table, but rows arrive with auth_seq out of label
    # order; ordering must follow label_seq_id.
    rows = CIF_BYTES.decode().splitlines()
    atom_rows = [r for r in rows if r.startswith(("ATOM", "HETATM"))]
    reordered = [atom_rows[1], atom_rows[0]] + atom_rows[2:]
    content = "\n".join([r for r in rows if not r.startswith(("ATOM", "HETATM"))]
                        + reordered).encode()
    structure = IngestionProcessor.run(content, "x.cif", "Audit", "predicted")
    from tos.engine.coord_math import get_sequential_residues
    order = get_sequential_residues(structure.atoms, "A")
    assert [s for s, _, _ in order] == [1]  # single GLY residue in chain A
    # atoms within the residue are keyed by name regardless of row order
    ca = [a for a in structure.atoms if a.atom_name == "CA"]
    assert len(ca) == 1


def test_run_autodetects_format():
    pdb_bytes = (b"HEADER    T\nEXPDTA    X-RAY DIFFRACTION\n"
                 b"ATOM      1  N   ALA A   1    "
                 b"   0.000   0.000   0.000  1.00 20.00           N  \nEND\n")
    s_pdb = IngestionProcessor.run(pdb_bytes, "any.pdb", "Audit", "experimental")
    assert any(a.atom_name == "N" for a in s_pdb.atoms)
    s_cif = IngestionProcessor.run(CIF_BYTES, "any.cif", "Audit", "predicted")
    assert any(a.atom_name == "CA" for a in s_cif.atoms)
