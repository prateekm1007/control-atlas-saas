"""
N8 v2 (Round-11 audit correction, 2026-10-20): per-residue alternate-
conformer selection is the DEFAULT.

The Round-11 audit rejected N8 v1 ("instruction was per-residue default;
v1 shipped it opt-in with per_atom as default"). v2 flips the engine
default: ALTLOC_POLICY is per_residue — one altloc ID per residue,
chosen by TOTAL occupancy (first-seen ID on ties); blank-altloc atoms
are shared and always kept. The original P0 item 6 per-atom rule is
preserved as an EXPLICIT opt-in (altloc_policy="per_residue" inverse:
altloc_policy="per_atom") for byte-level comparison with pre-N8
deliveries.

The default flip is an audit-directed semantic change, recorded in the
delivery register — not a silent change. Both policies guarantee
exactly one atom per (chain, residue, atom_name).
"""
from pathlib import Path

import numpy as np

from tos.ingestion.processor import (
    ALTLOC_POLICY,
    ALTLOC_POLICY_PER_ATOM,
    ALTLOC_POLICY_PER_RESIDUE,
    IngestionProcessor,
    _select_conformers,
)
from tos.engine.coord_math import get_atom_pos
from tos.engine.tier1_measurements import Tier1Measurements

FIXTURE = Path(__file__).parent / "fixtures" / "real_pdb_3erd.pdb"


def _atom(serial, name, altloc, res, chain, seq, x, y, z, occ,
          b=20.0, el="C"):
    return (f"ATOM  {serial:5d} {name:>4s}{altloc}{res:>3s} {chain}{seq:4d} "
            f"   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {el:>2s}")


def _pdb(*lines):
    return ("\n".join(lines) + "\n").encode()


# Residue A 10: CA prefers conformer A (0.60 vs 0.40), CB prefers
# conformer B (0.70 vs 0.30). Per-atom selection mixes them (chimera);
# per-residue totals are A=0.90 vs B=1.10, so B wins for the WHOLE residue.
CHIMERA_PDB = _pdb(
    _atom(1, "N", " ", "XRS", "A", 10, 1.0, 2.0, 3.0, 1.00, el="N"),
    _atom(2, "CA", "A", "XRS", "A", 10, 4.0, 5.0, 6.0, 0.60),
    _atom(3, "CA", "B", "XRS", "A", 10, 7.0, 8.0, 9.0, 0.40),
    _atom(4, "CB", "A", "XRS", "A", 10, 10.0, 11.0, 12.0, 0.30),
    _atom(5, "CB", "B", "XRS", "A", 10, 13.0, 14.0, 15.0, 0.70),
    "END",
)

# Equal totals (1.00 vs 1.00): first-seen conformer (A) must win.
TIE_PDB = _pdb(
    _atom(1, "CA", "A", "XRS", "A", 11, 1.0, 1.0, 1.0, 0.50),
    _atom(2, "CA", "B", "XRS", "A", 11, 2.0, 2.0, 2.0, 0.50),
    _atom(3, "CB", "A", "XRS", "A", 11, 3.0, 3.0, 3.0, 0.50),
    _atom(4, "CB", "B", "XRS", "A", 11, 4.0, 4.0, 4.0, 0.50),
    "END",
)

CHIMERA_CIF = b"""data_chimera
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
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
ATOM N N . XRS A 1 . 1.000 2.000 3.000 1.00 20.0 10 A 1
ATOM C CA A XRS A 1 . 4.000 5.000 6.000 0.60 20.0 10 A 1
ATOM C CA B XRS A 1 . 7.000 8.000 9.000 0.40 20.0 10 A 1
ATOM C CB A XRS A 1 . 10.000 11.000 12.000 0.30 20.0 10 A 1
ATOM C CB B XRS A 1 . 13.000 14.000 15.000 0.70 20.0 10 A 1
"""


def _pos(structure, name, seq=10):
    p = get_atom_pos(structure.atoms, "A", seq, name, "")
    assert p is not None, f"atom {name} missing"
    return p


# ── default policy: per_residue (Round-11 audit-corrected) ─────────────

def test_module_default_is_per_residue():
    assert ALTLOC_POLICY == ALTLOC_POLICY_PER_RESIDUE


def test_default_run_matches_explicit_per_residue():
    default = IngestionProcessor.run(CHIMERA_PDB, "x.pdb", "t", "experimental")
    explicit = IngestionProcessor.run(CHIMERA_PDB, "x.pdb", "t",
                                      "experimental",
                                      altloc_policy=ALTLOC_POLICY_PER_RESIDUE)
    assert [a.pos.tolist() for a in default.atoms] == \
           [a.pos.tolist() for a in explicit.atoms]


def test_default_is_chimera_safe_on_chimera_input():
    """Under the v2 DEFAULT (per_residue), the chimera input cannot
    produce a mixed residue: conformer B has the higher TOTAL occupancy
    (1.10 vs 0.90), so BOTH CA and CB come from B."""
    s = IngestionProcessor.run(CHIMERA_PDB, "x.pdb", "t", "experimental")
    assert np.allclose(_pos(s, "CA"), [7.0, 8.0, 9.0])     # conformer B
    assert np.allclose(_pos(s, "CB"), [13.0, 14.0, 15.0])  # conformer B


def test_per_atom_optin_produces_mixed_conformers():
    """Documents the EXPLICIT opt-in per_atom behavior on the chimera
    input (the original P0 item 6 rule): CA comes from conformer A,
    CB from conformer B. This is the pre-N8 semantics, retained for
    byte-level comparison and recorded here so the difference between
    the two policies is explicit rather than implied."""
    s = IngestionProcessor.run(CHIMERA_PDB, "x.pdb", "t", "experimental",
                               altloc_policy=ALTLOC_POLICY_PER_ATOM)
    assert np.allclose(_pos(s, "CA"), [4.0, 5.0, 6.0])   # conformer A
    assert np.allclose(_pos(s, "CB"), [13.0, 14.0, 15.0])  # conformer B


# ── per_residue policy ──────────────────────────────────────────────────

def test_per_residue_keeps_whole_residue_consistent():
    s = IngestionProcessor.run(CHIMERA_PDB, "x.pdb", "t", "experimental",
                               altloc_policy=ALTLOC_POLICY_PER_RESIDUE)
    # conformer B has the higher TOTAL occupancy (1.10 vs 0.90): BOTH
    # CA and CB must come from B.
    assert np.allclose(_pos(s, "CA"), [7.0, 8.0, 9.0])
    assert np.allclose(_pos(s, "CB"), [13.0, 14.0, 15.0])


def test_per_residue_keeps_blank_altloc_shared_atoms():
    s = IngestionProcessor.run(CHIMERA_PDB, "x.pdb", "t", "experimental",
                               altloc_policy=ALTLOC_POLICY_PER_RESIDUE)
    # N has a blank altloc: shared, retained under any conformer choice.
    assert np.allclose(_pos(s, "N"), [1.0, 2.0, 3.0])


def test_per_residue_tie_keeps_first_seen_conformer():
    s = IngestionProcessor.run(TIE_PDB, "t.pdb", "t", "experimental",
                               altloc_policy=ALTLOC_POLICY_PER_RESIDUE)
    assert np.allclose(_pos(s, "CA", seq=11), [1.0, 1.0, 1.0])  # A first
    assert np.allclose(_pos(s, "CB", seq=11), [3.0, 3.0, 3.0])  # A first


def test_per_residue_guarantees_one_atom_per_name():
    s = IngestionProcessor.run(CHIMERA_PDB, "x.pdb", "t", "experimental",
                               altloc_policy=ALTLOC_POLICY_PER_RESIDUE)
    keys = [(a.chain_id, a.res_seq, a.insertion_code, a.atom_name)
            for a in s.atoms]
    assert len(keys) == len(set(keys)), "duplicate atom records survived"


def test_unknown_policy_rejected():
    import pytest
    with pytest.raises(ValueError):
        _select_conformers([], "per_whatever")


# ── mmCIF path uses the same selector ───────────────────────────────────

def test_cif_per_residue_consistent():
    s = IngestionProcessor.run(CHIMERA_CIF, "x.cif", "t", "experimental",
                               altloc_policy=ALTLOC_POLICY_PER_RESIDUE)
    assert np.allclose(_pos(s, "CA"), [7.0, 8.0, 9.0])
    assert np.allclose(_pos(s, "CB"), [13.0, 14.0, 15.0])
    assert np.allclose(_pos(s, "N"), [1.0, 2.0, 3.0])


def test_cif_default_matches_per_residue():
    """mmCIF default (v2) is per_residue: B wins the residue (1.10 vs
    0.90), so CA and CB both come from B."""
    s = IngestionProcessor.run(CHIMERA_CIF, "x.cif", "t", "experimental")
    assert np.allclose(_pos(s, "CA"), [7.0, 8.0, 9.0])
    assert np.allclose(_pos(s, "CB"), [13.0, 14.0, 15.0])


def test_cif_per_atom_optin_matches_legacy():
    """mmCIF explicit per_atom opt-in keeps the original P0 item 6
    semantics: CA from conformer A (0.60), CB from conformer B (0.70)."""
    s = IngestionProcessor.run(CHIMERA_CIF, "x.cif", "t", "experimental",
                               altloc_policy=ALTLOC_POLICY_PER_ATOM)
    assert np.allclose(_pos(s, "CA"), [4.0, 5.0, 6.0])
    assert np.allclose(_pos(s, "CB"), [13.0, 14.0, 15.0])


# ── real 3ERD under the new policy ──────────────────────────────────────

def test_3erd_per_residue_no_duplicates_and_audit_runs():
    structure = IngestionProcessor.run(
        FIXTURE.read_bytes(), "3erd.pdb", "Audit", "experimental",
        altloc_policy=ALTLOC_POLICY_PER_RESIDUE)
    keys = [(a.chain_id, a.res_seq, a.insertion_code, a.atom_name)
            for a in structure.atoms]
    assert len(keys) == len(set(keys)), "duplicate atom records survived"
    results, coverage, _ = Tier1Measurements.run_full_audit(structure)
    assert results["LAW-120"]["status"] in ("PASS", "VETO")
    assert results["LAW-120"]["sample"] > 0


def test_3erd_default_per_residue_asp321_still_conformer_a():
    """3ERD ASP A 321 carries conformers A and B at 0.50/0.50 for the
    backbone atoms, so the per-residue TOTALS tie (A = B) and the
    per-residue rule keeps the FIRST-SEEN conformer — A, exactly as the
    pre-N8 per-atom tie rule did. The P0-era pinned assertion therefore
    still holds under the v2 default, and the outcome change is limited
    to residues whose per-atom occupancies disagree (chimera inputs)."""
    structure = IngestionProcessor.run(FIXTURE.read_bytes(), "3erd.pdb",
                                       "Audit", "experimental")
    n_pos = get_atom_pos(structure.atoms, "A", 321, "N", "")
    assert n_pos is not None
    assert np.allclose(n_pos, [30.044, -0.841, -5.710])  # the 'A' record
