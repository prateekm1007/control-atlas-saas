import hashlib as _hashlib

STATION_METADATA = {
    "version": "22.5.3",
    "institution": "TOSCANINI Structural Governance",
    "posture": "CONSERVATIVE",
    "min_coverage_pass": 70.0,
    "ml_sigma_multiplier": 1.25,
    "schema_version": "1.0.0"
}

# ═══════════════════════════════════════════════════════════════════════════════
# PIL-CHM-03: Chemical Canon (Unit Standardization)
#
# CALIBRATION NOTES (v22.5.2):
#
# LAW-100 (Bond Integrity):
#   Previous: 4.0σ — too strict for macromolecular crystallography.
#   At 1.7Å resolution, coordinate uncertainty is ~0.1-0.2Å, which translates
#   to z-scores of 5-15σ against Engh-Huber ideals (σ = 0.01-0.02Å).
#   MolProbity flags bonds at 4σ as *noteworthy* but does not reject structures
#   unless >5% of bonds exceed 4σ. We now measure the PERCENTAGE of outlier
#   bonds (z > 4σ) against total bonds, with threshold at 5%.
#   Reference: Engh & Huber (1991), Acta Cryst. A47, 392-400.
#
# LAW-110 (Backbone Gap):
#   Previous: 1.5Å — catches normal peptide bonds at medium resolution.
#   Inter-residue C(i)-N(i+1) peptide bond: ideal 1.329Å, but coordinate
#   uncertainty at 1.5-2.5Å resolution means observed values of 1.3-1.7Å
#   are routine. True backbone gaps (missing residues) show C-N > 3.0Å.
#   Recalibrated to 2.0Å: catches real gaps, ignores coordinate imprecision.
#   Reference: Word et al. (1999), J. Mol. Biol. 285, 1735-1747.
#
# LAW-160 (Chain Integrity):
#   Previous: 4.2Å — too tight for normal Cα spacing variation.
#   Trans-peptide Cα(i)-Cα(i+1): ideal 3.8Å, observed range 3.6-4.3Å
#   in well-refined structures. Cis-proline: ~2.9Å. True chain breaks
#   show CA-CA > 5.0Å. Recalibrated to 4.5Å.
#   Reference: Ramachandran & Sasisekharan (1968), Adv. Protein Chem. 23.
#
# LAW-182 (Hydrophobic Burial):
#   Previous: 0.3 ratio with 12Å radius — too dependent on protein size.
#   Multi-chain complexes (e.g., 4HHB tetramer) spread the centroid,
#   making legitimate hydrophobic cores appear unburied. Adjusted to
#   per-chain evaluation with 15Å radius.
#   NOTE: This remains heuristic (LAW-182 is classified heuristic).
# ═══════════════════════════════════════════════════════════════════════════════

LAW_CANON = {
    "LAW-100": {"title": "Bond Integrity", "threshold": 5.0, "unit": "%", "operator": "<=", "scope": "core residues", "type": "Percentage"},
    "LAW-105": {"title": "Reliability Coverage", "threshold": 70.0, "unit": "%", "operator": ">=", "scope": "total residues", "type": "Percentage"},
    "LAW-110": {"title": "Backbone Gap", "threshold": 2.0, "unit": "Å", "operator": "<=", "scope": "sequential pairs", "type": "Scalar"},
    "LAW-120": {"title": "Bond Angle", "threshold": 18.0, "unit": "°", "operator": "<=", "scope": "core residues", "type": "RMSD"},
    "LAW-125": {"title": "Ramachandran", "threshold": 5.0, "unit": "%", "operator": "<=", "scope": "core residues", "type": "Percentage"},
    "LAW-130": {"title": "Clashscore", "threshold": 20.0, "unit": "clashes/1000 atoms", "operator": "<=", "scope": "core atoms", "type": "Rate"},
    "LAW-135": {"title": "Omega Planarity", "threshold": 3.0, "unit": "%", "operator": "<=", "scope": "core residues", "type": "Percentage"},
    "LAW-145": {"title": "Chirality", "threshold": 0.0, "unit": "count", "operator": "=", "scope": "core residues", "type": "Count"},
    "LAW-150": {"title": "Rotamer Audit", "threshold": 20.0, "unit": "%", "operator": "<=", "scope": "core residues", "type": "Percentage"},
    "LAW-160": {"title": "Chain Integrity", "threshold": 4.5, "unit": "Å", "operator": "<=", "scope": "sequential C-alphas", "type": "Scalar"},
    "LAW-170": {"title": "Residue Identity", "threshold": 0.0, "unit": "count", "operator": "=", "scope": "total residues", "type": "Count"},
    "LAW-182": {"title": "Hydrophobic Burial", "threshold": 0.3, "unit": "ratio", "operator": ">=", "scope": "hydrophobic core", "type": "Rate"},
    "LAW-195": {"title": "Disulfide Geometry", "threshold": 0.20, "unit": "Å", "operator": "<=", "scope": "cys pairs", "type": "Scalar"},
    "LAW-155": {"title": "Voxel Occupancy", "threshold": 2.0, "unit": "V", "operator": ">=", "scope": "core atomic grid", "type": "Scalar"},
    "LAW-200": {"title": "Packing Quality", "threshold": 300.0, "unit": "Å³/atom", "operator": "<=", "scope": "core bounding box", "type": "Rate"}
}

LAW_METHOD_CLASSIFICATIONS = {lid: "deterministic" if lid not in ["LAW-155", "LAW-182", "LAW-200"] else "heuristic" for lid in LAW_CANON}
DETERMINISTIC_LAWS = sorted([l for l, m in LAW_METHOD_CLASSIFICATIONS.items() if m == "deterministic"])
HEURISTIC_LAWS = sorted([l for l, m in LAW_METHOD_CLASSIFICATIONS.items() if m == "heuristic"])
TOTAL_LAWS, DETERMINISTIC_COUNT, HEURISTIC_COUNT = len(LAW_CANON), len(DETERMINISTIC_LAWS), len(HEURISTIC_LAWS)

LAW_CATEGORIES = {lid: "geometric" if lid not in ["LAW-145", "LAW-150", "LAW-170", "LAW-182", "LAW-195"] else "chemical" for lid in LAW_CANON}
CATEGORY_DISPLAY = {"geometric": "GEOMETRIC INVARIANTS", "chemical": "CHEMICAL INVARIANTS"}

def get_laws_by_category():
    g = {}
    for l, c in LAW_CATEGORIES.items(): g.setdefault(c, []).append(l)
    return {k: sorted(v) for k, v in g.items()}

SCORE_DEFINITIONS = {
    "DETERMINISTIC_SCORE": {"title": "Deterministic Integrity", "explanation": "Physical compliance relative to crystallographic ideals."},
    "ADVISORY_SCORE": {"title": "Heuristic Advisory", "explanation": "Statistical proxies for structural plausibility."},
    "ML_CONFIDENCE": {"title": "ML Confidence", "explanation": "Model-reported mean pLDDT."},
    "STRATEGIC_SCORE": {"title": "EPI Priority Index", "explanation": "Non-probabilistic prioritization metric for lead selection."}
}

BAYESIAN_FORMULA = "P = clamp(S6, 0, 1) x W_arch x (1 - M_S8)"

ARCHITECTURE_WEIGHTS = {
    "NONE": 1.0, "LINEAR": 0.95, "MULTIVALENT": 0.85, "ENGAGEMENT": 0.80,
    "MOTIFS": 0.90, "LINKERS": 0.88, "CONFORMATIONAL": 0.75, "METAL": 0.70
}

SIGMA_TABLE = {"N-CA": 0.010, "CA-C": 0.021, "C-O": 0.020, "C-N": 0.014, "CA-CB": 0.020, "C-S": 0.025, "S-S": 0.016, "N-CA-C": 2.5, "CA-C-N": 1.5, "O-C-N": 1.6}
IDEAL_TABLE = {"N-CA": 1.470, "CA-C": 1.525, "C-O": 1.231, "C-N": 1.329, "CA-CB": 1.530, "C-S": 1.803, "S-S": 2.033, "N-CA-C": 111.0, "CA-C-N": 116.6, "O-C-N": 123.0}

STANDARD_RESIDUES = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
HYDROPHOBIC_RESIDUES = {"ALA","VAL","ILE","LEU","MET","PHE","TRP","PRO"}

# ═══════════════════════════════════════════════════════════════════════════
# WP-SCI-1: Measurement-status vocabulary (Constitution §4; directive §2)
#
# Every law row carries exactly one of these statuses:
#   PASS               measured and within threshold
#   VETO               measured deterministic violation
#   FAIL               measured failure of a non-deterministic check
#   NOT_MEASURED       empty sample — nothing was measurable (never a zero,
#                      never a favorable PASS)
#   NOT_APPLICABLE     the law's scope does not apply to this structure
#                      (e.g. LAW-195 with no disulfide pair) — recorded,
#                      not counted as passed, not a violation
#   MEASUREMENT_FAILED the measurement attempt itself failed (acquisition,
#                      tooling) — an explicit non-success state
# A missing mandatory measurement (NOT_MEASURED / MEASUREMENT_FAILED on a
# deterministic law) prevents an authoritative PASS; the adjudicator yields
# INDETERMINATE unless an independent measured VETO exists.
# ═══════════════════════════════════════════════════════════════════════════
STATUS_MEASURED_PASS        = "PASS"
STATUS_MEASURED_VETO        = "VETO"
STATUS_MEASURED_FAIL        = "FAIL"
STATUS_NOT_MEASURED         = "NOT_MEASURED"
STATUS_NOT_APPLICABLE       = "NOT_APPLICABLE"
STATUS_MEASUREMENT_FAILED   = "MEASUREMENT_FAILED"

BLOCKING_STATUSES = {STATUS_MEASURED_VETO, STATUS_MEASURED_FAIL}
# Non-success states that block an authoritative PASS without claiming a
# physical violation:
INDETERMINATE_CAUSE_STATUSES = {STATUS_NOT_MEASURED, STATUS_MEASUREMENT_FAILED}
# States that neither pass nor veto:
NON_PARTICIPATING_STATUSES = {STATUS_NOT_APPLICABLE}

# ═══════════════════════════════════════════════════════════════════════════
# WP-SCI-1: Law containment registry (directive §1)
#
# A contained law is one whose measurement definition or threshold is under
# repair and has not been re-validated by an independent calibration and
# validation protocol. A contained law CANNOT independently issue an
# authoritative VETO: if its failure would be the only blocking evidence,
# the adjudicator returns INDETERMINATE with explicit containment
# provenance instead of PASS or VETO. Containment is narrow, explicit,
# recorded, and lifted only by a reviewed policy change once the lift
# condition is met with evidence.
#
# LAW_CANON above remains byte-identical (canon hash 6a9cd4b4349b81de is
# unchanged); containment is a separate, auditable policy layer.
# ═══════════════════════════════════════════════════════════════════════════
LAW_CONTAINMENT = {
    "LAW-120": {
        "status": "contained",
        "reason": "measurement-definition-repair",
        "detail": (
            "Peptide O-C-N angle previously measured O(i)-C(i)-N(i) using the "
            "same residue's amide nitrogen; corrected to O(i)-C(i)-N(i+1) over "
            "valid sequential pairs (WP-SCI-1). The pinned 18.0 deg threshold "
            "was selected from the uncorrected metric on the v23 cohort and is "
            "NOT validated for the corrected metric."
        ),
        "recorded": "2026-10-19",
        "lift_condition": (
            "v24 diagnostic calibration recorded AND independent validation "
            "protocol (backend/calibration/protocol/evaluation_protocol_v1.md) "
            "acceptance criteria met on the held-out set AND reviewed policy "
            "update lifting containment"
        ),
        "effect": (
            "a measured LAW-120 VETO/FAIL alone cannot produce an authoritative "
            "VETO; verdict becomes INDETERMINATE with containment provenance"
        ),
    },
}

CONTAINED_LAWS = {lid for lid, c in LAW_CONTAINMENT.items() if c.get("status") == "contained"}


def compute_canon_hash():
    s = "".join([f"{k}|{LAW_METHOD_CLASSIFICATIONS[k]}" for k in sorted(LAW_CANON.keys())])
    return _hashlib.sha256(s.encode()).hexdigest()[:16]

LAW_CANON_HASH = compute_canon_hash()


# ═══════════════════════════════════════════════════════════════════════════
# WP-SCI-1: Semantic law hash (directive §6)
#
# The legacy canon hash above binds only sorted law IDs and method
# classifications — a threshold or measurement-definition change does not
# move it. compute_law_semantic_hash() binds MEANING: per-law IDs, titles,
# thresholds, operators, units, scopes, types, method classifications,
# containment records, the pinned ideal/sigma geometry tables, and the raw
# bytes digest of the modality-specific enforcement matrix (PIL-CAL-03).
#
# CI pins this digest (governance/law_semantic.sha256). Any meaning or
# threshold drift changes the digest and fails the gate until a reviewed
# update re-pins it. Historical results remain bound to their originally
# recorded policy hashes and are never reinterpreted under new semantics.
# ═══════════════════════════════════════════════════════════════════════════
def compute_law_semantic_hash() -> str:
    import json as _json
    from tos.governance.modality_matrix import compute_matrix_hash as _mmh

    fingerprint = {
        "schema": "law-semantic-fingerprint/1.0.0",
        "laws": [
            {
                "id": lid,
                "title": LAW_CANON[lid]["title"],
                "threshold": LAW_CANON[lid]["threshold"],
                "operator": LAW_CANON[lid]["operator"],
                "unit": LAW_CANON[lid]["unit"],
                "scope": LAW_CANON[lid]["scope"],
                "type": LAW_CANON[lid]["type"],
                "method_class": LAW_METHOD_CLASSIFICATIONS[lid],
                "category": LAW_CATEGORIES[lid],
                "containment": LAW_CONTAINMENT.get(lid),
            }
            for lid in sorted(LAW_CANON.keys())
        ],
        "ideal_table": {k: IDEAL_TABLE[k] for k in sorted(IDEAL_TABLE)},
        "sigma_table": {k: SIGMA_TABLE[k] for k in sorted(SIGMA_TABLE)},
        "modality_matrix_sha256": _mmh(),
        "status_vocabulary": sorted([
            STATUS_MEASURED_PASS, STATUS_MEASURED_VETO, STATUS_MEASURED_FAIL,
            STATUS_NOT_MEASURED, STATUS_NOT_APPLICABLE, STATUS_MEASUREMENT_FAILED,
        ]),
    }
    canonical = _json.dumps(fingerprint, sort_keys=True, separators=(",", ":"))
    return _hashlib.sha256(canonical.encode("utf-8")).hexdigest()
