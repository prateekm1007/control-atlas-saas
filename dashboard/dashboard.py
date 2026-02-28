"""
dashboard/dashboard.py
Toscanini Dashboard v2 — Forensic Structural Validation UI

Replaces legacy dashboard.py (Arena/Warhead edition).
All backend contracts verified against live API before writing.

LOCKED INVARIANTS:
  - /ingest      POST  Form(mode, candidate_id, t3_category) + File(optional)
  - /v1/batch    POST  Form(mode, t3_category) + File(zip, required)
  - /search      POST  Form(query)
  - /laws        GET   -> {"laws": {law_id: {title,threshold,unit,operator,scope,type}}}
  - /nkg         GET   -> {"vetoes": [], "successes": [], "total_vetoes": int, "total_successes": int}
  - /definitions GET   -> {score_key: {title, explanation}}

RESPONSE SHAPE (ToscaniniResponse):
  verdict.binary                    -> "PASS" | "VETO" | "INDETERMINATE"
  verdict.deterministic_score       -> int 0-100  [HERO METRIC]
  verdict.advisory_score            -> int 0-100
  verdict.physical_score            -> int 0-100
  verdict.confidence_score          -> float
  verdict.det_passed / det_total    -> int
  verdict.heur_passed / heur_total  -> int
  verdict.coverage_pct              -> float
  verdict.suppression_reason        -> str | None
  governance.audit_id               -> str
  governance.station_version        -> str
  governance.timestamp_utc          -> str
  governance.governance_fingerprint.canon_hash
  governance.governance_fingerprint.matrix_hash
  tier1.laws[]                      -> law_id, title, status, method, observed,
                                       threshold, operator, units, deviation,
                                       sample_size, scope, principle
  tier3.probability                 -> int 0-100
  characterization.total_atoms      -> int
  characterization.total_residues   -> int
  characterization.source_type      -> str
  characterization.resolution       -> float | None
  witness_reports.executive         -> str
  witness_reports.deep_dive         -> str
  witness_reports.recommendation    -> str
  ai_model_used                     -> str | None
  pdf_b64                           -> str
  pdb_b64                           -> str

BATCH SHAPE (BatchResponse):
  summary.total / passed / vetoed / indeterminate / errors -> int
  summary.mean_deterministic_score  -> float
  summary.common_failing_laws[]     -> {law_id: str, count: int}
  results[]                         -> filename, candidate_id, success,
                                       response (ToscaniniResponse|None), error

STYLE:
  apply_arena_theme()               -> injects CSS
  get_active_theme()                -> dict of color tokens
  get_active_theme_name()           -> str
  score_color(value, max_val=100)   -> hex color str
  render_lifecycle_header(state)    -> renders chip strip
  THEMES                            -> {"Biotech Noir": {...}, "Clinical White": {...}}

DEPENDENCIES (requirements.txt unchanged):
  streamlit>=1.29.0
  requests>=2.31.0
  py3Dmol>=2.0.4
"""

import streamlit as st
import requests
from remediation_generator import generate_remediation_zip, should_show_remediation
from comparison_engine import compare_audits
import base64
import os
import hashlib
import zipfile
import io
from datetime import datetime, timezone

from style_utils import (
    apply_arena_theme,
    get_active_theme,
    get_active_theme_name,
    score_color,
    render_lifecycle_header,
    THEMES,
)

# ===============================================================
# CONFIGURATION
# ===============================================================

BACKEND = os.getenv("BACKEND_URL", "http://brain:8000")
API_KEY = os.getenv("TOSCANINI_API_KEY", "")

T3_CATEGORIES = ["NONE"]  # Simplified — refinement studio replaces weighting

BENCHMARK_CRYSTALS = [
    {"pdb_id": "3NIR", "name": "Endothiapepsin (Ultra-High)", "resolution": 0.48, "method": "X-ray",   "af_id": "AF-P11838-F1"},
    {"pdb_id": "2VB1", "name": "Insulin (Atomic)",            "resolution": 0.65, "method": "X-ray",   "af_id": "AF-P01308-F1"},
    {"pdb_id": "1CRN", "name": "Crambin",                     "resolution": 1.50, "method": "X-ray",   "af_id": "AF-P01542-F1"},
    {"pdb_id": "4HHB", "name": "Hemoglobin (Deoxy)",          "resolution": 1.74, "method": "X-ray",   "af_id": "AF-P69905-F1"},
    {"pdb_id": "1UBQ", "name": "Ubiquitin",                   "resolution": 1.80, "method": "X-ray",   "af_id": "AF-P0CG47-F1"},
    {"pdb_id": "6LZG", "name": "SARS-CoV-2 RBD",              "resolution": 2.50, "method": "X-ray",   "af_id": "AF-P0DTC2-F1"},
    {"pdb_id": "7BV2", "name": "RdRp Complex (Cryo-EM)",      "resolution": 2.50, "method": "Cryo-EM", "af_id": "AF-P0DTD1-F1"},
    {"pdb_id": "1G03", "name": "Protein G (NMR)",             "resolution": 0.00, "method": "NMR",     "af_id": "AF-P06654-F1"},
]

# ===============================================================
# PAGE CONFIG — must be first Streamlit call
# ===============================================================

st.set_page_config(
    page_title="Toscanini · Structural Validation",
    layout="wide",
    page_icon="🔬",
)

# ===============================================================
# THEME BOOTSTRAP
# noir_toggle must exist in session_state before apply_arena_theme()
# ===============================================================

if "noir_toggle" not in st.session_state:
    st.session_state.noir_toggle = True

apply_arena_theme()

# ===============================================================
# SESSION STATE INITIALISATION
# ===============================================================

if "audit_result" not in st.session_state:
    st.session_state.audit_result = None
if "baseline_audit" not in st.session_state:
    st.session_state.baseline_audit = None
if "candidates" not in st.session_state:
    st.session_state.candidates = []
if "benchmark_results" not in st.session_state:
    st.session_state.benchmark_results = {}


# ===============================================================
# §0  THE SINGLE API CONDUIT
# Every outbound call passes through here.
# Injects X-API-Key. Returns parsed JSON or None on failure.
# No physics. No scoring. Transport only.
# ===============================================================

def api(method: str, path: str, **kwargs):
    headers = kwargs.pop("headers", {})
    if API_KEY:
        headers["X-API-Key"] = API_KEY
    try:
        r = requests.request(
            method,
            f"{BACKEND}{path}",
            timeout=180,
            headers=headers,
            **kwargs,
        )
        r.raise_for_status()
        return r.json()
    except requests.exceptions.ConnectionError:
        st.error(f"**Connection refused** -> `{BACKEND}{path}`\n\nIs the backend container running?")
        return None
    except requests.exceptions.HTTPError as exc:
        code = exc.response.status_code
        try:
            detail = exc.response.json().get("detail", exc.response.text[:300])
        except Exception:
            detail = exc.response.text[:300]
        st.error(f"**API {code}** on `{path}`: {detail}")
        return None
    except Exception as exc:
        st.error(f"**Unexpected error** on `{path}`: {exc}")
        return None


# ===============================================================
# §1  SHARED RENDERING HELPERS
# Display only. No math, no logic, no scoring.
# ===============================================================

def _binary_color(binary: str) -> str:
    return {"PASS": "#00CC66", "VETO": "#FF4444"}.get(binary, "#FFB833")


def render_verdict_banner(binary: str, det_score: int, det_passed: int, det_total: int) -> None:
    t = get_active_theme()
    color = _binary_color(binary)
    label = {"PASS": "PASS", "VETO": "VETO"}.get(binary, binary)
    icon  = {"PASS": "✅", "VETO": "🛑"}.get(binary, "⚠️")
    st.markdown(
        f"""
        <div style="
            background:linear-gradient(135deg,{color}1A,{color}08);
            border-left:4px solid {color};
            padding:1rem 1.5rem;
            border-radius:4px;
            margin-bottom:1rem;
            font-family:'Courier New',monospace;
        ">
            <span style="font-size:1.5rem;font-weight:800;color:{color};">
                {icon} {label}
            </span>
            <span style="margin-left:2rem;font-size:1rem;color:{t['text']};">
                Deterministic: <strong>{det_score}%</strong>
            </span>
            <span style="margin-left:1.5rem;font-size:0.85rem;color:{t['muted']};">
                {det_passed}/{det_total} laws satisfied
            </span>
        </div>
        """,
        unsafe_allow_html=True,
    )


def render_score_row(verdict: dict) -> None:
    det  = verdict.get("deterministic_score", 0)
    adv  = verdict.get("advisory_score", 0)
    phys = verdict.get("physical_score", 0)
    conf = verdict.get("confidence_score", 0.0)
    t    = get_active_theme()
    det_color = score_color(det)

    st.markdown(
        f"""
        <div style="
            background:{det_color}18;
            border:1px solid {det_color}55;
            border-radius:6px;
            padding:0.8rem 1.2rem;
            margin-bottom:0.75rem;
            font-family:'Courier New',monospace;
        ">
            <div style="font-size:0.7rem;color:{t['muted']};letter-spacing:0.08em;">
                DETERMINISTIC INTEGRITY
            </div>
            <div style="font-size:2.2rem;font-weight:800;color:{det_color};line-height:1.1;">
                {det}%
            </div>
            <div style="font-size:0.72rem;color:{t['muted']};margin-top:2px;">
                Physical law compliance relative to crystallographic ideals
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    c1, c2, c3 = st.columns(3)
    c1.metric("Advisory Score", f"{adv}%",        help="Statistical proxies for structural plausibility.")
    c2.metric("Physical Score", f"{phys}%",        help="Same basis as deterministic; present for traceability.")
    c3.metric("ML Confidence",  f"{round(conf,1)}", help="Model-reported mean pLDDT.")


def render_law_ledger(laws: list) -> None:
    if not laws:
        st.info("No law evaluations returned.")
        return

    METHOD_META = {
        "deterministic":         ("🟢", "Deterministic",           "Core physical invariant. Failure triggers a VETO."),
        "advisory_experimental": ("🟡", "Advisory (Experimental)", "Reported for transparency. Excluded from deterministic score."),
        "heuristic":             ("🔵", "Heuristic",               "Statistical proxy. Does not trigger VETO."),
    }

    for law in laws:
        status    = law.get("status", "UNKNOWN")
        method    = law.get("method", "deterministic")
        law_id    = law.get("law_id", "?")
        title     = law.get("title", law_id)
        observed  = law.get("observed", "—")
        threshold = law.get("threshold", "—")
        operator  = law.get("operator", "")
        units     = law.get("units", "")
        deviation = law.get("deviation", "—")
        scope     = law.get("scope", "—")
        size      = law.get("sample_size", "—")
        principle = law.get("principle", "N/A")

        status_icon           = "✅" if status == "PASS" else "🛑"
        m_icon, m_label, m_help = METHOD_META.get(method, ("⚪", method, ""))

        with st.expander(
            f"{status_icon} {law_id}: {title} — {status} [{m_icon} {m_label}]",
            expanded=False,
        ):
            if m_help:
                st.caption(f"ℹ️ {m_help}")
            col_a, col_b = st.columns(2)
            col_a.markdown(f"**Observed:** `{observed} {units}`")
            col_a.markdown(f"**Threshold:** `{operator} {threshold} {units}`")
            col_a.markdown(f"**Deviation:** {deviation}")
            col_b.markdown(f"**Scope:** {scope}")
            col_b.markdown(f"**Sample Size:** {size}")
            if principle and principle != "N/A":
                st.caption(f"Principle: {principle}")


def render_governance_fingerprint(governance: dict) -> None:
    with st.expander("🛡️ Governance Fingerprint", expanded=False):
        fp = governance.get("governance_fingerprint", {})
        c1, c2 = st.columns(2)
        c1.code(f"Canon Hash  : {fp.get('canon_hash', 'n/a')}")
        matrix_hash = fp.get('matrix_hash', 'n/a')
        c2.code(f"Matrix Hash : {matrix_hash[:32]}...")
        st.caption(
            f"Audit ID: `{governance.get('audit_id', 'n/a')}` · "
            f"Engine: v{governance.get('station_version', '?')} · "
            f"{governance.get('timestamp_utc', '')[:19].replace('T', ' ')} UTC"
        )
        st.caption(
            f"Schema: {fp.get('matrix_schema_version', '?')} · "
            f"Policy: {fp.get('policy_ref', '?')}"
        )


def render_characterization(char: dict) -> None:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total Atoms",    char.get("total_atoms",    "—"))
    c2.metric("Total Residues", char.get("total_residues", "—"))
    c3.metric("Source Type",    char.get("source_type",    "—"))
    res = char.get("resolution")
    c4.metric("Resolution", f"{res:.2f} Å" if res else "N/A (NMR/Predicted)")



def _extract_plddt_from_pdb(pdb_b64_str):
    """Extract per-residue B-factors (pLDDT) from base64 PDB CA atoms."""
    import base64 as _b64
    confidences = []
    try:
        pdb_text = _b64.b64decode(pdb_b64_str).decode("utf-8", errors="ignore")
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and len(line) >= 66:
                if line[12:16].strip() == "CA":
                    try:
                        confidences.append(float(line[60:66]))
                    except (ValueError, IndexError):
                        confidences.append(50.0)
    except Exception:
        pass
    return confidences


def _render_sidebar_confidence_strip(confidences):
    """Render a per-residue colored strip as pure HTML in the sidebar."""
    n = len(confidences)
    if n < 3:
        return
    colors = []
    for c in confidences:
        if c >= 90: colors.append("#0066CC")
        elif c >= 70: colors.append("#4DA6FF")
        elif c >= 50: colors.append("#FFB84D")
        else: colors.append("#CC0000")
    stops = ",".join(["%s %.1f%% %.1f%%" % (colors[i], i/n*100, (i+1)/n*100) for i in range(n)])
    n_vh = sum(1 for c in confidences if c >= 90)
    n_h = sum(1 for c in confidences if 70 <= c < 90)
    n_m = sum(1 for c in confidences if 50 <= c < 70)
    n_l = sum(1 for c in confidences if c < 50)
    # Find contiguous low regions
    low_spans = []
    in_r, start = False, 0
    for i, c in enumerate(confidences):
        if c < 50:
            if not in_r: start, in_r = i, True
        elif in_r:
            if i - start >= 3: low_spans.append((start+1, i))
            in_r = False
    if in_r and n - start >= 3: low_spans.append((start+1, n))
    low_html = ""
    if low_spans:
        spans_str = ", ".join(["%d-%d" % (s, e) for s, e in low_spans])
        low_html = '<div style="font-size:10px;color:#CC0000;margin-top:4px;">Disordered: residues %s</div>' % spans_str
    html = (
        '<div style="margin:8px 0;">'
        '<div style="font-size:11px;font-weight:600;color:#333;margin-bottom:4px;">'
        'Per-Residue Confidence (%d residues)</div>'
        '<div style="width:100%%;height:14px;border-radius:3px;'
        'background:linear-gradient(to right,%s);border:1px solid #ddd;"></div>'
        '<div style="display:flex;justify-content:space-between;font-size:9px;color:#999;margin-top:2px;">'
        '<span>1</span><span>%d</span></div>'
        '<div style="display:flex;gap:8px;margin-top:4px;flex-wrap:wrap;">'
        '<span style="font-size:9px;"><span style="color:#0066CC;">&#9632;</span> VH:%d</span>'
        '<span style="font-size:9px;"><span style="color:#4DA6FF;">&#9632;</span> H:%d</span>'
        '<span style="font-size:9px;"><span style="color:#FFB84D;">&#9632;</span> M:%d</span>'
        '<span style="font-size:9px;"><span style="color:#CC0000;">&#9632;</span> L:%d</span>'
        '</div>%s</div>'
    ) % (n, stops, n, n_vh, n_h, n_m, n_l, low_html)
    st.markdown(html, unsafe_allow_html=True)


def render_3d_viewer(pdb_b64: str, binary: str, height: int = 500, width: int = 700) -> None:
    try:
        pdb_str = base64.decodebytes(pdb_b64.encode()).decode("utf-8", errors="ignore")
        if len(pdb_str.strip()) < 50:
            st.warning("3D viewer: PDB payload too small to render.")
            return
        import py3Dmol
        import streamlit.components.v1 as components
        t = get_active_theme()
        view = py3Dmol.view(width=width, height=height)
        view.addModel(pdb_str, "pdb")
        if binary == "PASS":
            view.setStyle({"cartoon": {"color": "spectrum"}})
        else:
            view.setStyle({
                "stick":  {"colorscheme": "orangeCarbon"},
                "sphere": {"radius": 0.5},
            })
        view.setBackgroundColor(t["viz_bg"])
        view.zoomTo()
        components.html(view._make_html(), height=height + 20, width=width + 20, scrolling=False)
    except ImportError:
        st.warning("py3Dmol not installed — 3D viewer unavailable.")
    except Exception as exc:
        st.warning(f"Visualizer unavailable: {exc}")


def render_witness_reports(witness: dict, ai_model: str | None) -> None:
    with st.expander("🤖 AI Witness Analysis", expanded=False):
        executive = witness.get("executive", "")
        if executive:
            st.info(executive)
        else:
            st.caption("Executive summary unavailable.")
        col_a, col_b = st.columns(2)
        with col_a:
            st.markdown("**Deep Dive Analysis**")
            st.write(witness.get("deep_dive", "Analysis pending."))
        with col_b:
            st.markdown("**Governance Recommendation**")
            st.write(witness.get("recommendation", "Awaiting human review."))
        if ai_model:
            st.caption(f"Narrative Model: {ai_model}")


def render_pdf_download(pdf_b64: str, audit_id: str) -> None:
    if pdf_b64:
        pdf_bytes = base64.b64decode(pdf_b64)
        st.download_button(
            label="📄 Download Forensic Dossier (PDF)",
            data=pdf_bytes,
            file_name=f"TOSCANINI_DOSSIER_{audit_id}.pdf",
            mime="application/pdf",
            use_container_width=True,
        )
    else:
        st.error("PIL-PAR-14 VIOLATION: PDF byte-stream not found in payload.")


def render_full_audit_result(res: dict) -> None:
    v       = res.get("verdict", {})
    binary  = v.get("binary", "UNKNOWN")
    gov     = res.get("governance", {})
    char    = res.get("characterization", {})
    witness = res.get("witness_reports") or {}
    laws    = res.get("tier1", {}).get("laws", [])

    render_verdict_banner(
        binary,
        v.get("deterministic_score", 0),
        v.get("det_passed", 0),
        v.get("det_total", 0),
    )
    render_score_row(v)

    st.divider()
    st.subheader("Structure Characterization")
    render_characterization(char)

    st.divider()
    st.subheader("📋 Tier-1 Diagnostic Ledger")
    st.caption(
        "🟢 Deterministic = governs PASS/VETO  |  "
        "🟡 Advisory = reported only, not scored  |  "
        "🔵 Heuristic = statistical proxy"
    )
    render_law_ledger(laws)

    st.divider()
    if witness:
        render_witness_reports(witness, res.get("ai_model_used"))

    render_governance_fingerprint(gov)

    st.divider()
    st.subheader("🧬 Forensic 3D Reconstruction")
    pdb_b64 = res.get("pdb_b64", "")
    if pdb_b64:
        render_3d_viewer(pdb_b64, binary, height=550, width=780)
    else:
        st.info("No PDB payload returned.")

    st.divider()
    render_pdf_download(res.get("pdf_b64", ""), gov.get("audit_id", "ARTIFACT"))


def render_law_comparison_table(crystal_laws: list, af_laws: list) -> None:
    t     = get_active_theme()
    c_map = {l["law_id"]: l for l in (crystal_laws or [])}
    a_map = {l["law_id"]: l for l in (af_laws or [])}
    all_ids = sorted(set(list(c_map.keys()) + list(a_map.keys())))

    if not all_ids:
        st.info("No laws to compare.")
        return

    header = (
        f'<tr style="background:{t["surface"]};border-bottom:2px solid {t["border"]};">'
        f'<th style="padding:8px;color:{t["text"]};text-align:left;">Law</th>'
        f'<th style="padding:8px;color:{t["text"]};text-align:left;">Title</th>'
        f'<th style="padding:8px;color:{t["text"]};text-align:center;">Crystal</th>'
        f'<th style="padding:8px;color:{t["text"]};text-align:center;">AlphaFold</th>'
        f'<th style="padding:8px;color:{t["text"]};text-align:right;">Threshold</th>'
        f'</tr>'
    )

    def _cell(law_row):
        if not law_row:
            return f'<span style="color:{t["muted"]};">—</span>'
        s     = law_row.get("status", "?")
        val   = law_row.get("observed", "?")
        units = law_row.get("units", "")
        color = "#00CC66" if s == "PASS" else "#FF4444"
        badge = ""
        if law_row.get("method") == "advisory_experimental":
            badge = (
                ' <span style="font-size:0.65rem;color:#D4A017;background:#3a3000;'
                'padding:1px 4px;border-radius:3px;">ADV</span>'
            )
        return f'<span style="color:{color};font-weight:bold;">{val} {units}</span> ({s}){badge}'

    rows = ""
    for lid in all_ids:
        c        = c_map.get(lid)
        a        = a_map.get(lid)
        row_data = c or a or {}
        title    = row_data.get("title", lid)
        thresh   = row_data.get("threshold", "—")
        op       = row_data.get("operator", "")
        units    = row_data.get("units", "")
        rows += (
            f'<tr style="border-bottom:1px solid {t["border"]};">'
            f'<td style="padding:6px;color:{t["text"]};font-family:monospace;">{lid}</td>'
            f'<td style="padding:6px;color:{t["text"]};">{title}</td>'
            f'<td style="padding:6px;text-align:center;">{_cell(c)}</td>'
            f'<td style="padding:6px;text-align:center;">{_cell(a)}</td>'
            f'<td style="padding:6px;text-align:right;color:{t["muted"]};font-family:monospace;">'
            f'{op} {thresh} {units}</td>'
            f'</tr>'
        )

    st.markdown(
        f'<table style="width:100%;border-collapse:collapse;">{header}{rows}</table>',
        unsafe_allow_html=True,
    )


# ===============================================================
# §2  SIDEBAR
# ===============================================================

with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/1042/1042307.png", width=72)
    st.title("TOSCANINI")

    st.toggle("Dark Mode", key="noir_toggle")
    apply_arena_theme()

    # t3_category hardcoded — no user selection
    t3_category = "NONE"

    st.divider()

    # ── ENGINE STATUS ──
    try:
        health = requests.get(f"{BACKEND}/health", timeout=3).json()
        engine_ver = health.get("version", "?")
        engine_status = health.get("status", "?").upper()
        if engine_status == "OPERATIONAL":
            st.success(f"Engine v{engine_ver} -- {engine_status}")
        else:
            st.warning(f"Engine v{engine_ver} -- {engine_status}")
    except Exception:
        st.error("Engine unreachable")
        engine_ver = "?"

    st.divider()

    # ═══════════════════════════════════════════
    # STRUCTURE REFINEMENT STUDIO
    # ═══════════════════════════════════════════
    _ar = st.session_state.get("audit_result")

    if _ar is not None:
        _v = _ar.get("verdict", {})
        _binary = _v.get("binary", "UNKNOWN")
        _cov = _v.get("coverage_pct", 0)
        _det_passed = _v.get("det_passed", 0)
        _det_total = _v.get("det_total", 12)
        _det_score = _v.get("deterministic_score", 0)
        _laws = _ar.get("tier1", {}).get("laws", [])
        _failing = [l for l in _laws if l.get("status") != "PASS" and l.get("method") == "deterministic"]
        _char = _ar.get("characterization", {})
        _total_res = _char.get("total_residues", "N/A")
        _audit_id = _ar.get("governance", {}).get("audit_id", "ARTIFACT")

        st.markdown("### REFINEMENT STUDIO")

        # ── VERDICT STATUS INDICATOR ──
        if _binary == "PASS" and _cov >= 70:
            st.success("Structure is high-confidence -- ready for drug discovery")
        elif _binary == "INDETERMINATE":
            st.warning("INDETERMINATE -- refinement recommended")
        elif _binary == "VETO":
            st.error("VETO -- critical physics violations detected")
        else:
            st.info(f"Status: {_binary}")

        # ── DIAGNOSIS SUMMARY ──
        st.markdown("**Diagnosis**")
        _diag_cols = st.columns(2)
        _diag_cols[0].metric("Coverage", f"{_cov:.1f}%",
                             delta="OK" if _cov >= 70 else "LOW",
                             delta_color="normal" if _cov >= 70 else "inverse")
        _diag_cols[1].metric("Det. Laws", f"{_det_passed}/{_det_total}",
                             delta="OK" if _det_passed == _det_total else f"{_det_total - _det_passed} FAIL",
                             delta_color="normal" if _det_passed == _det_total else "inverse")

        # Show failing laws
        if _failing:
            st.markdown("**Failing Laws:**")
            for _fl in _failing:
                st.markdown(
                    f"- **{_fl['law_id']}** {_fl.get('title', '')}: "
                    f"{_fl.get('observed', 'N/A')} "
                    f"(threshold: {_fl.get('operator', '')} {_fl.get('threshold', 'N/A')})"
                )

        st.divider()

        # ── PER-RESIDUE CONFIDENCE MAP ──
        _pdb_b64 = _ar.get("pdb_b64", "")
        if _pdb_b64:
            _residue_conf = _extract_plddt_from_pdb(_pdb_b64)
            if len(_residue_conf) >= 3:
                _render_sidebar_confidence_strip(_residue_conf)

        st.divider()

        # ── RECOMMENDED REFINEMENT ACTIONS ──
        _needs_refinement = (_binary != "PASS" or _cov < 70)

        if _needs_refinement:
            st.markdown("**Recommended Refinement Actions**")
            st.caption("Physics-first methods to fix structural issues")

            # Priority table
            _refine_data = [
                {"Priority": "1", "Method": "Rosetta Fast Relax",
                 "Expected Improvement": "Fix rotamer outliers and clashes",
                 "Est. Time": "30-90 sec", "Recommendation": "Strongly recommended"},
                {"Priority": "2", "Method": "Short MD Equilibration",
                 "Expected Improvement": "Improve hydrophobic burial",
                 "Est. Time": "3-8 min", "Recommendation": "Recommended"},
                {"Priority": "3", "Method": "Targeted Loop Modeling",
                 "Expected Improvement": "Rebuild disordered regions",
                 "Est. Time": "2-5 min", "Recommendation": "If loops critical"},
                {"Priority": "4", "Method": "Pocket Refinement",
                 "Expected Improvement": "Optimize binding site geometry",
                 "Est. Time": "8-15 min", "Recommendation": "For drug design"},
            ]

            # Determine which methods are most relevant
            _has_rotamer = any("rotam" in l.get("title", "").lower() or "155" in l.get("law_id", "") for l in _failing)
            _has_burial = any("burial" in l.get("title", "").lower() or "hydrophob" in l.get("title", "").lower() or "182" in l.get("law_id", "") for l in _failing)
            _low_coverage = _cov < 50

            for _rd in _refine_data:
                _method = _rd["Method"]
                _rec = _rd["Recommendation"]

                # Highlight most relevant
                if _method == "Rosetta Fast Relax" and _has_rotamer:
                    _rec = "CRITICAL -- rotamer issues detected"
                elif _method == "Short MD Equilibration" and _has_burial:
                    _rec = "CRITICAL -- burial issues detected"
                elif _method == "Targeted Loop Modeling" and _low_coverage:
                    _rec = "CRITICAL -- low coverage regions"

                _rd["Recommendation"] = _rec

            st.dataframe(_refine_data, use_container_width=True, hide_index=True)

            st.divider()

            # ── REMEDIATION PACKAGE DOWNLOAD ──
            st.markdown("**Remediation Package**")
            st.caption(
                "Download a parameterized remediation kit: Rosetta FastRelax protocol, "
                "OpenMM equilibration script, and per-law prescription report. "
                "Execute in your validated environment, then re-upload for re-certification."
            )
            try:
                _zip_bytes = generate_remediation_zip(_ar)
                _zip_name = f"TOSCANINI_REMEDIATION_{_audit_id}.zip"
                st.download_button(
                    label="Download Remediation Package (ZIP)",
                    data=_zip_bytes,
                    file_name=_zip_name,
                    mime="application/zip",
                    use_container_width=True,
                    help="Contains: Rosetta XML, flags, OpenMM script, per-law report, README",
                )
                st.caption(
                    "Package includes: rosetta_relax.xml, rosetta.flags, "
                    "openmm_equilibrate.py, loop_modeling.xml, remediation_report.json, README.txt"
                )
            except Exception as _rem_err:
                st.error(f"Remediation package generation failed: {_rem_err}")

            st.divider()

            # ── EXPECTED OUTCOME ──
            st.markdown("**Expected After Refinement**")
            if _cov < 70:
                st.caption(
                    f"Current coverage: {_cov:.1f}%. "
                    "After Rosetta relax + short MD, expected improvement to 65-85% "
                    "depending on initial structure quality."
                )
            if _failing:
                _fail_ids = ", ".join([f['law_id'] for f in _failing])
                st.caption(f"Failing: {_fail_ids}. Fast Relax resolves rotamer/clash issues in ~80% of cases.")
            st.caption("Re-upload refined structure to Toscanini for re-validation.")

            st.divider()

            # ── B1: CALLBACK TOKEN ──────────────────────────────────────
            st.markdown("**🎫 Managed Refinement (B1 — Callback)**")
            st.caption("Refine locally, upload result, get automatic before/after comparison.")

            _email_input = st.text_input(
                "Email (optional)",
                key="b1_email",
                placeholder="you@institution.edu"
            )

            if st.button("Generate Callback Token", key="b1_token_btn"):
                try:
                    import sys
                    sys.path.insert(0, '/app')
                    from tos.security.tokens import create_refinement_token
                    _token = create_refinement_token(
                        _audit_id,
                        _email_input if _email_input else None
                    )
                    # Store token in session state for BYOC notebook generator
                    st.session_state["last_callback_token"] = _token
                    _zip_managed = generate_remediation_zip(_ar, callback_token=_token)
                    st.success("✅ Token generated and embedded in ZIP")
                    st.download_button(
                        label="📦 Download Managed Remediation Package",
                        data=_zip_managed,
                        file_name=f"TOSCANINI_REMEDIATION_{_audit_id}_MANAGED.zip",
                        mime="application/zip",
                        use_container_width=True,
                        key="b1_managed_zip"
                    )
                    with st.expander("View Token", expanded=False):
                        st.code(_token, language="text")
                        st.caption("Token valid 7 days. Single-use. See README.txt in ZIP.")
                except Exception as _te:
                    st.error(f"Token generation failed: {_te}")

            st.divider()

            # ── B2: RUN WITH TOSCANINI ───────────────────────────────────
            st.markdown("**⚡ Run with Toscanini (B2 — Managed GPU)**")
            st.caption("We execute refinement on our GPU. Auto re-audit on completion.")

            _protocol_choice = st.selectbox(
                "Protocol",
                options=["auto", "openmm", "rosetta"],
                key="b2_protocol",
                help="Auto = Toscanini selects based on violations"
            )

            _b2_email = st.text_input(
                "Email for notification",
                key="b2_email",
                placeholder="you@institution.edu"
            )

            if st.button("⚡ Run with Toscanini (Beta)", key="b2_run_btn", type="primary"):
                try:
                    import requests as _req

                    # Re-upload the original PDB for GPU execution
                    if "uploaded_file_bytes" in st.session_state and st.session_state.uploaded_file_bytes:
                        _pdb_bytes = st.session_state.uploaded_file_bytes
                        _pdb_name  = st.session_state.get("uploaded_file_name", "structure.pdb")

                        _resp = _req.post(
                            "http://brain:8000/refinement/submit",
                            data={
                                "audit_id":   _audit_id,
                                "protocol":   _protocol_choice,
                                "user_email": _b2_email if _b2_email else ""
                            },
                            files={"file": (_pdb_name, _pdb_bytes, "application/octet-stream")},
                            timeout=30
                        )

                        if _resp.status_code == 200:
                            _job = _resp.json()
                            st.success(f"✅ Job submitted: `{_job['job_id']}`")
                            st.info(
                                "Protocol: **" + str(_job.get("protocol","?")) + "**  \n"
                                + "Est. time: **" + str(_job.get("estimated_minutes","?")) + " min**  \n"
                                + "Poll: /refinement/status/" + str(_job.get("job_id","?"))
                            )
                            st.session_state["b2_job_id"] = _job["job_id"]
                        else:
                            st.warning(_resp.json().get("message", "GPU worker not yet deployed."))
                    else:
                        st.warning("Original file not in session. Re-upload to use managed execution.")

                except Exception as _b2e:
                    st.warning(f"B2 GPU worker not deployed yet: {_b2e}")

            # ── B2: JOB STATUS POLLING ───────────────────────────────────
            if "b2_job_id" in st.session_state and st.session_state.b2_job_id:
                _job_id = st.session_state.b2_job_id
                st.markdown(f"**Job Status: `{_job_id}`**")

                if st.button("🔄 Refresh Status", key="b2_refresh"):
                    try:
                        import requests as _req2
                        _status_resp = _req2.get(
                            f"http://brain:8000/refinement/status/{_job_id}",
                            timeout=10
                        )
                        _status = _status_resp.json()
                        _state  = _status.get("state", "unknown")

                        _state_colors = {
                            "queued":    "🟡",
                            "running":   "🔵",
                            "success":   "🟢",
                            "failed":    "🔴",
                            "timeout":   "🟠",
                            "not_found": "⚪"
                        }
                        _icon = _state_colors.get(_state, "⚪")

                        st.markdown(f"{_icon} **State:** {_state.upper()}")

                        if _status.get("logs"):
                            st.caption(f"Logs: {_status['logs']}")

                        if _state == "success" and _status.get("comparison_url"):
                            st.success("✅ Refinement complete!")
                            st.markdown(
                                f"[View Before/After Comparison]({_status['comparison_url']})"
                            )
                            st.session_state["b2_job_id"] = None

                        if _state == "failed":
                            st.error(f"❌ Job failed: {_status.get('error', 'Unknown error')}")
                            st.session_state["b2_job_id"] = None

                    except Exception as _se:
                        st.warning(f"Status check failed: {_se}")

        # ── PDF DOWNLOAD (always visible) ──
        st.divider()
        _pdf_b64 = _ar.get("pdf_b64", "")
        if _pdf_b64:
            import base64 as _b64
            _pdf_bytes = _b64.b64decode(_pdf_b64)
            st.download_button(
                label="Download Forensic Dossier (PDF)",
                data=_pdf_bytes,
                file_name=f"TOSCANINI_DOSSIER_{_audit_id}.pdf",
                mime="application/pdf",
                use_container_width=True,
            )

    else:
        # No audit result — minimal sidebar
        st.markdown("### REFINEMENT STUDIO")
        st.caption("Upload a structure to activate diagnosis and refinement tools.")
        st.caption("The studio auto-detects issues and recommends physics-first fixes.")

    st.divider()
    st.caption(f"Theme: {get_active_theme_name()}")


# ===============================================================
# §3  HEADER + LIFECYCLE STRIP
# ===============================================================

t = get_active_theme()
st.markdown(
    f'<h1 style="color:{t["accent"]};text-align:center;font-family:Courier New;'
    f'margin-bottom:0px;">🔬 TOSCANINI</h1>',
    unsafe_allow_html=True,
)
st.markdown(
    f'<p style="text-align:center;color:{t["muted"]};font-size:0.8rem;'
    f'margin-top:0;font-family:Courier New;">Forensic Structural Validation Engine</p>',
    unsafe_allow_html=True,
)

if st.session_state.audit_result:
    _lifecycle_state = "AUDITED"
elif st.session_state.candidates:
    _lifecycle_state = "ACQUIRING"
else:
    _lifecycle_state = "INSTANTIATED"

render_lifecycle_header(_lifecycle_state)


# ===============================================================
# §4  TAB LAYOUT
# ===============================================================

tab_audit, tab_batch, tab_refine, tab_bench, tab_ref, tab_history = st.tabs(
    ["⚡ Audit", "📦 Batch", "⚡ Refinement", "🔬 Benchmark", "📜 Reference", "📊 History"]
)


# ===============================================================
# TAB 1 — AUDIT (Single Structure)
# ===============================================================

with tab_audit:
    if st.session_state.audit_result:
        res = st.session_state.audit_result

        if st.button("🔄 Return to Command Center", use_container_width=True):
            st.session_state.audit_result = None
            st.session_state.candidates   = []
            st.rerun()

        st.divider()
        render_full_audit_result(res)

        # ── BEFORE/AFTER COMPARISON ──────────────────────────────────────
        st.divider()
        col_baseline, col_compare = st.columns(2)
        
        with col_baseline:
            if st.button("📌 Set as Baseline for Comparison", use_container_width=True,
                         help="Store this audit as the 'before' state for refinement comparison"):
                st.session_state.baseline_audit = res
                st.success(f"Baseline stored: {res.get('governance', {}).get('audit_id', 'N/A')}")
        
        with col_compare:
            if st.session_state.baseline_audit:
                baseline_id = st.session_state.baseline_audit.get("governance", {}).get("audit_id", "N/A")
                st.info(f"Baseline: {baseline_id}")
            else:
                st.caption("No baseline set")
        
        # If baseline exists and current != baseline, show comparison
        if st.session_state.baseline_audit:
            baseline_id = st.session_state.baseline_audit.get("governance", {}).get("audit_id", "")
            current_id = res.get("governance", {}).get("audit_id", "")
            
            if baseline_id != current_id:
                st.divider()
                st.subheader("🔬 Before/After Comparison")
                
                try:
                    comparison = compare_audits(st.session_state.baseline_audit, res)
                    
                    # Summary metrics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        delta_cov = comparison["coverage_delta"]
                        st.metric("Coverage Change", f"{delta_cov:+.1f}%",
                                  delta=f"{delta_cov:.1f}%",
                                  delta_color="normal" if delta_cov >= 0 else "inverse")
                    with col2:
                        delta_viol = comparison["violation_count_delta"]
                        st.metric("Violations", f"{comparison['violation_count_after']}",
                                  delta=f"{delta_viol:+d}",
                                  delta_color="inverse" if delta_viol < 0 else "normal")
                    with col3:
                        st.metric("Verdict Change", comparison["verdict_change"],
                                  delta="Improved" if comparison["verdict_improved"] else "No improvement",
                                  delta_color="normal" if comparison["verdict_improved"] else "off")
                    
                    st.divider()
                    
                    # Regressions warning
                    if comparison["regressions"]:
                        st.error(f"⚠️ **REGRESSIONS DETECTED:** {len(comparison['regressions'])} deterministic laws worsened")
                        st.caption(", ".join(comparison["regressions"]))
                    
                    # Improvements
                    if comparison["improvements"]:
                        st.success(f"✅ **IMPROVEMENTS:** {len(comparison['improvements'])} violations resolved")
                        st.caption(", ".join(comparison["improvements"]))
                    
                    # Law-by-law table
                    if comparison["law_changes"]:
                        st.markdown("**Detailed Law Changes**")
                        change_df = []
                        for lc in comparison["law_changes"]:
                            status_emoji = {
                                "RESOLVED": "✅",
                                "REGRESSED": "❌",
                                "CHANGED": "🔄",
                                "UNCHANGED": "➖",
                                "STRUCTURAL_CHANGE": "⚠️",
                            }.get(lc["change_type"], "")
                            
                            change_df.append({
                                "Law": lc["law_id"],
                                "Title": lc["title"],
                                "Before": lc["before_status"],
                                "After": lc["after_status"],
                                "Status": f"{status_emoji} {lc['change_type']}",
                            })
                        st.dataframe(change_df, use_container_width=True, hide_index=True)
                    
                    st.caption(
                        f"Baseline: {comparison['baseline_audit_id']} | "
                        f"Refined: {comparison['refined_audit_id']}"
                    )
                
                except Exception as e:
                    st.error(f"Comparison failed: {e}")

    else:
        st.subheader("Upload Structure for Audit")
        uploaded_file = st.file_uploader(
            "Upload a PDB or CIF file",
            type=["pdb", "cif"],
            key="audit_file_upload",
            help="File is validated deterministically through the full 15-law canon.",
        )

        if uploaded_file is not None:
            raw_name = uploaded_file.name
            stem = raw_name.rsplit(".", 1)[0] if "." in raw_name else raw_name
            if stem.strip():
                candidate_id = stem.strip()
            else:
                candidate_id = "UPLOAD_" + hashlib.md5(uploaded_file.getvalue()).hexdigest()[:8].upper()

            st.caption(f"Candidate ID: `{candidate_id}` · Weighting: `{t3_category}`")

            if st.button("⚡ Run Audit", type="primary", use_container_width=True):
                with st.spinner(f"Validating `{candidate_id}` through 15-law canon…"):
                    res = api(
                        "POST", "/ingest",
                        data={"mode": "Audit", "candidate_id": candidate_id, "t3_category": t3_category},
                        files={"file": (uploaded_file.name, uploaded_file.getvalue(), "application/octet-stream")},
                    )
                if res:
                    st.session_state.audit_result = res
                    # Store file bytes for B2 managed execution
                    st.session_state.uploaded_file_bytes = uploaded_file.getvalue()
                    st.session_state.uploaded_file_name  = uploaded_file.name
                    st.rerun()

        st.divider()

        with st.expander("🔍 Or fetch by protein name / UniProt / AlphaFold ID", expanded=False):
            q = st.text_input(
                "Search query",
                placeholder="e.g. insulin, EGFR, p53, AF-P01308-F1",
                label_visibility="collapsed",
            )
            if st.button("Search Global Repositories", type="secondary", use_container_width=True):
                if q.strip():
                    with st.spinner("Searching…"):
                        resp = api("POST", "/search", data={"query": q.strip()})
                    if resp and resp.get("results"):
                        st.session_state.candidates = resp["results"]
                    else:
                        st.warning("No results returned.")
                        st.session_state.candidates = []

            if st.session_state.candidates:
                st.divider()
                for i, c in enumerate(st.session_state.candidates):
                    col_label, col_btn = st.columns([5, 1])
                    col_label.markdown(f"**{c.get('id', '?')}** — {c.get('label', '')}")
                    if col_btn.button("Audit", key=f"audit_candidate_{i}"):
                        with st.spinner(f"Fetching and auditing `{c['id']}`…"):
                            res = api(
                                "POST", "/ingest",
                                data={
                                    "mode":         "Discovery",
                                    "candidate_id": c["id"],
                                    "t3_category":  t3_category,
                                },
                            )
                        if res:
                            st.session_state.audit_result = res
                            st.rerun()


# ===============================================================
# TAB 2 — BATCH (Industrial Validation)
# ===============================================================

with tab_batch:
    st.subheader("Batch Industrial Validation")
    st.caption(
        "Upload a ZIP archive containing PDB files. "
        "Each structure is validated independently through the full deterministic canon. "
        "Maximum 100 structures per batch."
    )

    zip_file = st.file_uploader(
        "Upload ZIP archive",
        type=["zip"],
        key="batch_zip_upload",
        help="Each file inside the ZIP must be a valid .pdb file.",
    )

    if zip_file is not None:
        try:
            with zipfile.ZipFile(io.BytesIO(zip_file.getvalue())) as zf:
                pdb_names = [
                    n for n in zf.namelist()
                    if n.lower().endswith(".pdb") and not n.startswith("__")
                ]
            struct_count = len(pdb_names)
            if struct_count > 100:
                st.warning(
                    f"⚠️ Archive contains {struct_count} PDB files. "
                    "Backend enforces a maximum of 100. The first 100 will be processed."
                )
            else:
                st.info(f"Detected **{struct_count}** PDB file(s) in archive.")
        except Exception:
            st.warning("Could not inspect ZIP contents client-side.")

        if st.button("🚀 Run Batch Validation", type="primary", use_container_width=True):
            with st.spinner("Processing batch — this may take several minutes…"):
                batch_res = api(
                    "POST", "/v1/batch",
                    data={"mode": "Audit", "t3_category": t3_category},
                    files={"file": (zip_file.name, zip_file.getvalue(), "application/zip")},
                )

            if batch_res:
                summary = batch_res.get("summary", {})
                results = batch_res.get("results", [])

                st.divider()
                st.subheader("Batch Summary")
                s1, s2, s3, s4, s5 = st.columns(5)
                s1.metric("Total",             summary.get("total", "—"))
                s2.metric("✅ Passed",         summary.get("passed", "—"))
                s3.metric("🛑 Vetoed",         summary.get("vetoed", "—"))
                s4.metric("⚠️ Indeterminate",  summary.get("indeterminate", "—"))
                s5.metric("❌ Errors",         summary.get("errors", "—"))

                mean_score = summary.get("mean_deterministic_score", 0)
                mean_color = score_color(mean_score)
                t_ = get_active_theme()
                st.markdown(
                    f'<div style="margin:0.5rem 0;font-family:Courier New;">'
                    f'Mean Deterministic Score: '
                    f'<span style="color:{mean_color};font-weight:700;font-size:1.1rem;">'
                    f'{mean_score:.1f}%</span></div>',
                    unsafe_allow_html=True,
                )

                failing = summary.get("common_failing_laws", [])
                if failing:
                    st.markdown("**Most Common Failing Laws:**")
                    total_structs = max(summary.get("total", 1), 1)
                    for fl in failing[:5]:
                        law_id  = fl.get("law_id", "?")
                        count   = fl.get("count", 0)
                        bar_pct = int(count / total_structs * 100)
                        st.markdown(
                            f'<div style="margin:3px 0;font-family:monospace;font-size:0.85rem;">'
                            f'<span style="color:{t_["text"]};">{law_id}</span>'
                            f'<div style="background:{t_["border"]};border-radius:3px;height:8px;margin-top:2px;">'
                            f'<div style="background:#FF4444;width:{bar_pct}%;height:8px;border-radius:3px;"></div>'
                            f'</div>'
                            f'<span style="color:{t_["muted"]};font-size:0.75rem;">{count} failures</span>'
                            f'</div>',
                            unsafe_allow_html=True,
                        )

                if results:
                    st.divider()
                    st.subheader("Per-Structure Results")

                    table_rows = []
                    for r in results:
                        if r.get("success") and r.get("response"):
                            v      = r["response"].get("verdict", {})
                            binary = v.get("binary", "ERROR")
                            det    = v.get("deterministic_score", 0)
                            cov    = v.get("coverage_pct", 0)
                        else:
                            binary = "ERROR"
                            det    = 0
                            cov    = 0

                        table_rows.append({
                            "Filename":       r.get("filename", "?"),
                            "Candidate ID":   r.get("candidate_id", "?"),
                            "Verdict":        binary,
                            "Det. Score (%)": det,
                            "Coverage (%)":   round(cov, 1),
                            "Error":          r.get("error") or "",
                        })

                    st.dataframe(table_rows, use_container_width=True)

                    keys = ["Filename", "Candidate ID", "Verdict", "Det. Score (%)", "Coverage (%)", "Error"]
                    csv_lines = [",".join(keys)]
                    for row in table_rows:
                        csv_lines.append(",".join(str(row.get(k, "")) for k in keys))
                    csv_bytes = "\n".join(csv_lines).encode("utf-8")

                    st.download_button(
                        "⬇️ Download Results as CSV",
                        data=csv_bytes,
                        file_name=f"toscanini_batch_{datetime.now(timezone.utc):%Y%m%d_%H%M%S}.csv",
                        mime="text/csv",
                    )


# ===============================================================
# TAB 3 — BENCHMARK (Experimental vs Predicted)
# ===============================================================

def _run_benchmark_audit(candidate_id: str) -> dict | None:
    return api(
        "POST", "/ingest",
        data={"mode": "Benchmark", "candidate_id": candidate_id, "t3_category": "NONE"},
    )



# ===============================================================
# TAB 3 — REFINEMENT (Upload Refined Structure)
# ===============================================================

with tab_refine:
    st.subheader("Upload Refined Structure")
    st.caption(
        "If you have refined your structure locally (Rosetta/OpenMM), "
        "upload it here with your callback token for automatic re-audit and before/after comparison."
    )

    st.markdown("---")

    col1, col2 = st.columns([3, 1])
    with col1:
        _ref_token = st.text_input(
            "Callback Token",
            key="refine_token",
            placeholder="Paste token from README.txt in your remediation package",
            help="Token is valid for 7 days, single-use"
        )
    with col2:
        st.write("")
        st.write("")
        _ref_file = st.file_uploader(
            "Refined PDB",
            type=["pdb", "cif"],
            key="refine_file"
        )

    if st.button("🚀 Submit Refined Structure", type="primary", key="refine_submit"):
        if not _ref_token:
            st.error("⚠️ Please provide a callback token")
        elif not _ref_file:
            st.error("⚠️ Please upload a refined PDB file")
        else:
            with st.spinner("Uploading and re-auditing refined structure..."):
                try:
                    _ref_resp = api(
                        "POST", "/refinement/callback",
                        data={"token": _ref_token},
                        files={"file": (_ref_file.name, _ref_file.getvalue(), "application/octet-stream")}
                    )
                    if _ref_resp:
                        _ref_verdict  = _ref_resp.get("verdict", "UNKNOWN")
                        _ref_coverage = _ref_resp.get("coverage_pct", 0)
                        _ref_score    = _ref_resp.get("deterministic_score", 0)
                        _orig_id      = _ref_resp.get("original_audit_id", "?")
                        _new_id       = _ref_resp.get("refined_audit_id", "?")

                        st.success("✅ Refined structure audited successfully!")

                        c1, c2, c3 = st.columns(3)
                        c1.metric("Verdict",   _ref_verdict)
                        c2.metric("Coverage",  f"{_ref_coverage:.1f}%")
                        c3.metric("Det. Score",f"{_ref_score}/100")

                        st.info(
                            "**Comparison ready:**  \n"
                            + f"Original: `{_orig_id}`  \n"
                            + f"Refined:  `{_new_id}`"
                        )

                        st.session_state["comparison_baseline"] = _orig_id
                        st.session_state["comparison_refined"]   = _new_id
                        st.balloons()
                    else:
                        st.error("Upload failed. Check token validity and file format.")
                except Exception as _re:
                    st.error(f"Upload error: {_re}")

    st.divider()

    # ── BYOC NOTEBOOK GENERATOR ──────────────────────────────────────
    # BYOC_PLATFORM_SELECTOR
    st.markdown("**🚀 Run on Free GPU — Kaggle / Colab / Paperspace**")
    st.caption(
        "No GPU? No problem. Generate a pre-filled notebook and run it on "
        "Kaggle (30 hrs/week free), Google Colab, or Paperspace. "
        "Your refined structure posts back automatically."
    )

    if st.session_state.get("audit_result"):
        _byoc_audit = st.session_state["audit_result"]

        _byoc_col1, _byoc_col2 = st.columns(2)

        with _byoc_col1:
            _byoc_platform = st.selectbox(
                "Platform",
                options=["kaggle", "colab", "paperspace"],
                format_func=lambda x: {
                    "kaggle":     "🟦 Kaggle (30 hrs/week free)",
                    "colab":      "🟠 Google Colab (free T4 GPU)",
                    "paperspace": "🟣 Paperspace Gradient (free tier)",
                }[x],
                key="byoc_platform"
            )

        with _byoc_col2:
            try:
                from notebook_generator import select_protocol_for_audit
                _auto_proto = select_protocol_for_audit(_byoc_audit)
            except Exception:
                _auto_proto = "openmm"

            _byoc_protocol = st.selectbox(
                "Protocol",
                options=["openmm", "rosetta"],
                index=0 if _auto_proto == "openmm" else 1,
                format_func=lambda x: {
                    "openmm":  "⚗️ OpenMM (no license needed)",
                    "rosetta": "🔬 Rosetta FastRelax (academic license)",
                }[x],
                key="byoc_protocol",
                help="Auto-selected based on your failing laws. Override if needed."
            )

        # Platform info card
        try:
            from notebook_generator import get_platform_info
            _pinfo = get_platform_info(_byoc_platform)
            st.info(
                f"**{_pinfo.get('name')}** · "
                f"{_pinfo.get('gpu_type')} · "
                f"{_pinfo.get('free_hours')} · "
                f"{_pinfo.get('notes')}"
            )
        except Exception:
            pass

        # Generate + download notebook
        _byoc_token = st.session_state.get("last_callback_token", "")
        if not _byoc_token:
            st.caption(
                "💡 Generate a callback token first via the remediation package "
                "(Download Remediation ZIP above) to enable automatic re-upload."
            )

        if st.button("📓 Generate Notebook", key="byoc_generate_btn", type="primary"):
            with st.spinner("Generating notebook..."):
                try:
                    # ── B5: Credit check before generation ───────────────
                    # byoc_credit_check
                    import sys as _nbsys
                    _nbsys.path.insert(0, "/app/gpu_worker")
                    try:
                        from worker.credits import check_credits, deduct_credits
                        _nb_identifier = (
                            st.session_state.get("b1_email") or
                            st.session_state.get("b2_email") or
                            "anonymous"
                        )
                        _nb_credit_check = check_credits(_nb_identifier, "notebook")
                        if not _nb_credit_check["allowed"]:
                            st.error(
                                f"No GPU credits remaining for `{_nb_identifier}`. "
                                f"Each notebook export costs 1 credit. "
                                f"Credits remaining: {_nb_credit_check['credits_remaining']}"
                            )
                            st.stop()
                    except ImportError:
                        # Credits module not available in this environment — allow
                        _nb_identifier  = "anonymous"
                        _nb_credit_check = {"allowed": True}

                    from notebook_generator import generate_notebook_bytes, get_platform_redirect_url
                    _nb_bytes = generate_notebook_bytes(
                        _byoc_audit,
                        _byoc_token or "PASTE_YOUR_CALLBACK_TOKEN_HERE",
                        _byoc_platform,
                        _byoc_protocol
                    )
                    _nb_filename = f"toscanini_{_byoc_audit.get('governance',{}).get('audit_id','AUDIT')}_{_byoc_protocol}.ipynb"

                    st.download_button(
                        label=f"⬇️ Download {_nb_filename}",
                        data=_nb_bytes,
                        file_name=_nb_filename,
                        mime="application/json",
                        key="byoc_download_btn"
                    )

                    # ── B5: Deduct credit after successful generation ──────
                    try:
                        _nb_audit_id = _byoc_audit.get("governance", {}).get("audit_id", "UNKNOWN")
                        deduct_credits(_nb_identifier, "notebook", f"NB_{_nb_audit_id}")
                        _remaining = _nb_credit_check["credits_remaining"] - 1
                        st.caption(f"✓ 1 credit deducted · {_remaining} remaining")
                    except Exception:
                        pass  # Non-fatal — notebook already generated

                    # Platform-specific instructions
                    if _byoc_platform == "kaggle":
                        st.success(
                            "**Kaggle steps:**  \n"
                            "1. Download the notebook above  \n"
                            "2. Go to [kaggle.com/kernels](https://www.kaggle.com/kernels)  \n"
                            "3. Click **New Notebook** → **File** → **Import Notebook**  \n"
                            "4. Upload the .ipynb file  \n"
                            "5. Enable GPU: **Settings** → **Accelerator** → T4 GPU  \n"
                            "6. Click **Run All** — results post back automatically"
                        )
                    elif _byoc_platform == "colab":
                        st.success(
                            "**Colab steps:**  \n"
                            "1. Download the notebook above  \n"
                            "2. Go to [colab.research.google.com](https://colab.research.google.com)  \n"
                            "3. **File** → **Upload notebook**  \n"
                            "4. **Runtime** → **Change runtime type** → T4 GPU  \n"
                            "5. **Runtime** → **Run all** — results post back automatically"
                        )
                    elif _byoc_platform == "paperspace":
                        st.success(
                            "**Paperspace steps:**  \n"
                            "1. Download the notebook above  \n"
                            "2. Go to [console.paperspace.com](https://console.paperspace.com)  \n"
                            "3. Create a new Notebook project  \n"
                            "4. Upload the .ipynb file  \n"
                            "5. Select free GPU runtime and run all cells"
                        )

                    if _byoc_token:
                        st.caption(f"Token embedded: `{_byoc_token[:20]}...` (single-use, 7 days)")
                    else:
                        st.warning(
                            "No callback token embedded. After running the notebook, "
                            "manually download the refined PDB and upload it above."
                        )

                except Exception as _be:
                    st.error(f"Notebook generation failed: {_be}")
    else:
        st.info("Run an audit first to enable notebook generation.")

    st.divider()

    # Job status polling for B2 managed jobs
    st.markdown("**⚡ B2 Job Status**")
    st.caption("If you submitted a managed GPU job, check its status here.")

    _poll_job_id = st.text_input(
        "Job ID",
        key="poll_job_id",
        placeholder="e.g. A1B2C3D4",
        value=st.session_state.get("b2_job_id", "")
    )

    if st.button("🔄 Check Job Status", key="poll_status_btn"):
        if not _poll_job_id:
            st.warning("Enter a Job ID")
        else:
            try:
                _st_resp = api("GET", f"/refinement/status/{_poll_job_id}")
                if _st_resp:
                    _state = _st_resp.get("state", "unknown")
                    _icons = {
                        "queued":    "🟡 QUEUED",
                        "running":   "🔵 RUNNING",
                        "success":   "🟢 SUCCESS",
                        "failed":    "🔴 FAILED",
                        "timeout":   "🟠 TIMEOUT",
                        "not_found": "⚪ NOT FOUND"
                    }
                    st.markdown(f"**Status:** {_icons.get(_state, _state.upper())}")

                    if _st_resp.get("logs"):
                        st.caption(f"Logs: {_st_resp['logs']}")
                    if _st_resp.get("error"):
                        st.error(f"Error: {_st_resp['error']}")
                    if _state == "success" and _st_resp.get("comparison_url"):
                        st.success("Refinement complete!")
                        st.markdown(f"Comparison URL: `{_st_resp['comparison_url']}`")
                    if _st_resp.get("started_at"):
                        st.caption(f"Started: {_st_resp['started_at']}")
                    if _st_resp.get("completed_at"):
                        st.caption(f"Completed: {_st_resp['completed_at']}")
            except Exception as _pe:
                st.error(f"Status check failed: {_pe}")

    # ── AUTO-COMPARISON DISPLAY ────────────────────────────────────
    st.divider()
    st.markdown("**📊 Before/After Comparison**")

    if "comparison_baseline" in st.session_state and "comparison_refined" in st.session_state:
        _base_id = st.session_state.get("comparison_baseline")
        _refn_id = st.session_state.get("comparison_refined")

        if _base_id and _refn_id:
            st.info(f"Comparing: `{_base_id}` → `{_refn_id}`")

            try:
                _cmp_resp = api(
                    "POST", "/compare",
                    data={"baseline_id": _base_id, "refined_id": _refn_id}
                )

                if _cmp_resp and _cmp_resp.get("status") == "success":
                    _cmp = _cmp_resp["comparison"]

                    # High-level metrics
                    c1, c2, c3 = st.columns(3)
                    c1.metric(
                        "Verdict",
                        _cmp.get("verdict_change", "?"),
                        delta="Improved" if _cmp.get("verdict_improved") else "No change"
                    )
                    c2.metric(
                        "Coverage",
                        f"+{_cmp.get('coverage_delta', 0):.1f}%" if _cmp.get("coverage_delta", 0) > 0 else f"{_cmp.get('coverage_delta', 0):.1f}%"
                    )
                    c3.metric(
                        "Violations",
                        f"{_cmp.get('violation_count_after', '?')}",
                        delta=str(_cmp.get("violation_count_delta", 0))
                    )

                    # Improvements
                    if _cmp.get("improvements"):
                        st.success("✅ Laws resolved: " + ", ".join(_cmp["improvements"]))
                    if _cmp.get("regressions"):
                        st.error("⚠️ Regressions: " + ", ".join(_cmp["regressions"]))

                    # Law-by-law breakdown
                    _changes = _cmp.get("law_changes", [])
                    if _changes:
                        with st.expander("🔬 Law-Level Changes", expanded=False):
                            for _ch in _changes:
                                _icon = "✅" if _ch.get("change_type") == "RESOLVED" else "⚠️"
                                st.markdown(
                                    f"{_icon} **{_ch['law_id']}** ({_ch.get('title','')}): "
                                    f"{_ch['before_status']} → {_ch['after_status']} "
                                    f"(observed: {_ch.get('before_observed','?')} → {_ch.get('after_observed','?')})"
                                )

                                # Residue improvements
                                _ri = _ch.get("residue_improvements")
                                if _ri and _ri.get("fixed_count", 0) > 0:
                                    st.caption(
                                        f"  Fixed residues: {', '.join(_ri['fixed_residues'][:10])}"
                                        + (f" (+ {_ri['fixed_count'] - 10} more)" if _ri['fixed_count'] > 10 else "")
                                    )
                                if _ri and _ri.get("still_count", 0) > 0:
                                    st.caption(
                                        f"  Still failing: {', '.join(_ri['still_failing'][:10])}"
                                    )
                else:
                    st.caption("Comparison not available yet. Both audits must be stored.")

            except Exception as _ce:
                st.caption(f"Comparison engine: {_ce}")
    else:
        st.caption("Submit a refined structure above to see before/after comparison.")

with tab_bench:
    st.subheader("Benchmark — Experimental vs. Predicted")
    st.markdown(
        "Each structure passes through the **identical** deterministic pipeline. "
        "Crystal structures at high resolution must PASS to confirm engine calibration."
    )
    st.divider()

    for entry in BENCHMARK_CRYSTALS:
        pdb_id = entry["pdb_id"]
        af_id  = entry["af_id"]

        with st.expander(
            f"**{entry['name']}** — {pdb_id} "
            f"({'%.2f' % entry['resolution']} Å, {entry['method']})",
            expanded=False,
        ):
            cached = st.session_state.benchmark_results.get(pdb_id, {})

            if not cached.get("crystal"):
                if st.button(f"⚡ Audit {pdb_id} + {af_id}", key=f"bench_{pdb_id}"):
                    with st.spinner(f"Auditing {pdb_id} (crystal)…"):
                        crystal_res = _run_benchmark_audit(pdb_id)

                    af_res = None
                    if crystal_res:
                        with st.spinner(f"Auditing {af_id} (AlphaFold)…"):
                            af_res = _run_benchmark_audit(af_id)
                        if af_res and "verdict" not in af_res:
                            af_res = None

                    st.session_state.benchmark_results[pdb_id] = {
                        "crystal":   crystal_res,
                        "alphafold": af_res,
                    }
                    st.rerun()

            else:
                cr       = cached["crystal"]
                ar       = cached.get("alphafold")
                cv       = cr.get("verdict", {})
                binary_c = cv.get("binary", "ERROR")

                col_c, col_a = st.columns(2)

                with col_c:
                    st.markdown(f"**🔵 Crystal: {pdb_id}**")
                    det_laws_c = [l for l in cr.get("tier1", {}).get("laws", []) if l.get("method") == "deterministic"]
                    adv_laws_c = [l for l in cr.get("tier1", {}).get("laws", []) if "advisory" in l.get("method", "")]
                    det_pass_c = sum(1 for l in det_laws_c if l.get("status") == "PASS")
                    if binary_c == "PASS":
                        st.success(
                            f"✅ {binary_c} — {det_pass_c}/{len(det_laws_c)} Det"
                            + (f" + {len(adv_laws_c)} Advisory" if adv_laws_c else "")
                        )
                    elif binary_c == "VETO":
                        st.error(f"🛑 {binary_c}")
                    else:
                        st.warning(f"⚠️ {binary_c}")
                    st.metric("Deterministic Score", f"{cv.get('deterministic_score', 0)}%")
                    st.metric("Coverage",            f"{cv.get('coverage_pct', 0):.1f}%")
                    st.caption(f"Method: {entry['method']} | Resolution: {entry['resolution']} Å")

                with col_a:
                    st.markdown(f"**🟠 AlphaFold: {af_id}**")
                    if ar:
                        av       = ar.get("verdict", {})
                        binary_a = av.get("binary", "ERROR")
                        if binary_a == "PASS":
                            st.success(f"✅ {binary_a} — {av.get('det_passed','?')}/{av.get('det_total','?')} Det")
                        elif binary_a == "VETO":
                            st.error(f"🛑 {binary_a}")
                        else:
                            st.warning(f"⚠️ {binary_a}")
                        st.metric("Deterministic Score", f"{av.get('deterministic_score', 0)}%")
                        st.metric("pLDDT (Core)",        f"{round(av.get('confidence_score', 0), 1)}")
                    else:
                        st.warning("⚠️ AlphaFold model unavailable for this target.")

                if cr.get("pdb_b64"):
                    st.divider()
                    st.markdown("**3D Structural Comparison**")
                    v3d_l, v3d_r = st.columns(2)
                    with v3d_l:
                        st.caption(f"Crystal: {pdb_id}")
                        render_3d_viewer(cr["pdb_b64"], binary_c, height=380, width=420)
                    with v3d_r:
                        if ar and ar.get("pdb_b64"):
                            binary_a_3d = ar.get("verdict", {}).get("binary", "ERROR")
                            st.caption(f"AlphaFold: {af_id}")
                            render_3d_viewer(ar["pdb_b64"], binary_a_3d, height=380, width=420)
                        else:
                            st.caption(f"AlphaFold: {af_id}")
                            st.info("Predicted structure not available for 3D comparison.")

                st.divider()
                st.markdown("**Deterministic Law Comparison**")
                render_law_comparison_table(
                    cr.get("tier1", {}).get("laws", []),
                    ar.get("tier1", {}).get("laws", []) if ar else [],
                )

                if cr.get("pdf_b64"):
                    st.divider()
                    audit_id = cr.get("governance", {}).get("audit_id", pdb_id)
                    render_pdf_download(cr["pdf_b64"], audit_id)


# ===============================================================
# TAB 4 — REFERENCE (Unified Source of Truth)
# ===============================================================

with tab_ref:
    st.subheader("Reference — Structural Law Canon & Knowledge Graph")

    ref_laws, ref_nkg, ref_defs = st.tabs(
        ["📐 Law Canon", "🧠 Knowledge Graph", "📖 Definitions"]
    )

    with ref_laws:
        with st.spinner("Loading Law Canon…"):
            canon = api("GET", "/laws")

        if canon:
            laws_dict = canon.get("laws", {})
            c1, c2 = st.columns(2)
            c1.metric("Total Laws",          canon.get("total_laws", len(laws_dict)))
            c2.metric("Deterministic Count", canon.get("deterministic_count", "?"))
            st.divider()
            for lid, info in laws_dict.items():
                with st.expander(f"**{lid}** — {info.get('title', '?')}", expanded=False):
                    col_a, col_b, col_c = st.columns(3)
                    col_a.metric(
                        "Threshold",
                        f"{info.get('operator','')} {info.get('threshold','?')} {info.get('unit','')}",
                    )
                    col_b.metric("Type",  info.get("type",  "?"))
                    col_c.metric("Scope", info.get("scope", "?"))

    with ref_nkg:
        with st.spinner("Loading Knowledge Graph…"):
            nkg = api("GET", "/nkg")

        if nkg:
            c1, c2 = st.columns(2)
            c1.metric("Total Vetoes",    nkg.get("total_vetoes",    0))
            c2.metric("Total Successes", nkg.get("total_successes", 0))
            st.divider()
            vetoes    = nkg.get("vetoes",    [])
            successes = nkg.get("successes", [])
            if vetoes:
                with st.expander(f"🛑 Veto Records ({len(vetoes)})", expanded=False):
                    st.json(vetoes)
            else:
                st.info("No veto records. Knowledge graph is empty.")
            if successes:
                with st.expander(f"✅ Success Records ({len(successes)})", expanded=False):
                    st.json(successes)

    with ref_defs:
        with st.spinner("Loading Definitions…"):
            defs = api("GET", "/definitions")

        if defs:
            for key, info in defs.items():
                with st.expander(f"**{info.get('title', key)}**", expanded=False):
                    st.write(info.get("explanation", "No explanation."))
                    st.caption(f"Key: `{key}`")
        else:
            st.info("Definitions unavailable.")


# ===============================================================
# TAB 6 — INSTITUTIONAL HISTORY
# ===============================================================

with tab_history:
    st.subheader("📊 Refinement History")
    st.caption(
        "Track your refinement journey. Every structure you have submitted "
        "for refinement appears here with before/after metrics."
    )

    st.markdown("---")

    # ── Identifier input ──────────────────────────────────────────────────────
    _hist_col1, _hist_col2 = st.columns([3, 1])
    with _hist_col1:
        _hist_identifier = st.text_input(
            "Email or IP",
            key="hist_identifier",
            placeholder="you@institution.edu  or  your IP address",
            help="Use the same email you provided when submitting refinement jobs"
        )
    with _hist_col2:
        st.write("")
        st.write("")
        _hist_use_ip = st.checkbox("Use my IP", key="hist_use_ip",
                                   help="Look up history by IP address instead of email")

    if st.button("🔍 Load History", key="hist_load_btn", type="primary"):
        if not _hist_identifier and not _hist_use_ip:
            st.warning("Enter an email address or check 'Use my IP'")
        else:
            with st.spinner("Loading refinement history..."):
                try:
                    if _hist_use_ip:
                        # Get client IP via a simple echo endpoint or use placeholder
                        _hist_resp = api("GET", "/history/ip/client")
                    else:
                        import urllib.parse
                        _encoded = urllib.parse.quote(_hist_identifier, safe="")
                        _hist_resp = api("GET", f"/history/{_encoded}")

                    if _hist_resp and _hist_resp.get("status") == "success":
                        st.session_state["history_data"] = _hist_resp
                    else:
                        st.error("Could not load history. Check identifier.")
                except Exception as _he:
                    st.error(f"History load failed: {_he}")

    # ── Display history ───────────────────────────────────────────────────────
    if st.session_state.get("history_data"):
        _hdata = st.session_state["history_data"]
        _hcount = _hdata.get("refinement_count", 0)
        _hcomps = _hdata.get("comparisons", [])

        # Summary metrics
        st.markdown("---")
        _hm1, _hm2, _hm3, _hm4 = st.columns(4)
        _hm1.metric("Total Refinements", _hcount)

        # Calculate aggregate stats
        _improved = sum(
            1 for c in _hcomps
            if c.get("refined_verdict") == "PASS" and c.get("original_verdict") != "PASS"
        )
        _still_veto = sum(
            1 for c in _hcomps
            if c.get("refined_verdict") == "VETO"
        )
        _protocols = {}
        for c in _hcomps:
            proto = c.get("protocol", "unknown")
            _protocols[proto] = _protocols.get(proto, 0) + 1

        _hm2.metric("Resolved to PASS", _improved,
                    delta=f"{_improved/_hcount*100:.0f}%" if _hcount > 0 else "0%")
        _hm3.metric("Still VETO", _still_veto)
        _hm4.metric("Most Used Protocol",
                    max(_protocols, key=_protocols.get) if _protocols else "—")

        st.markdown("---")

        if _hcount == 0:
            st.info(
                "No refinement history found for this identifier. "
                "Submit a structure for refinement to start tracking your journey."
            )
        else:
            st.markdown(f"**{_hcount} refinement(s) found — newest first**")

            for _i, _comp in enumerate(_hcomps):
                _orig_id   = _comp.get("original_audit_id", "?")
                _ref_id    = _comp.get("refined_audit_id", "?")
                _proto     = _comp.get("protocol", "unknown")
                _method    = _comp.get("refinement_method", "unknown")
                _created   = _comp.get("created_at", "")[:19].replace("T", " ")
                _orig_v    = _comp.get("original_verdict", "?")
                _ref_v     = _comp.get("refined_verdict", "?")
                _orig_s    = _comp.get("original_score", "?")
                _ref_s     = _comp.get("refined_score", "?")
                _orig_cov  = _comp.get("original_coverage", "?")
                _ref_cov   = _comp.get("refined_coverage", "?")

                # Verdict change indicator
                _verdict_icon = "✅" if _ref_v == "PASS" and _orig_v != "PASS" else (
                    "🔴" if _ref_v == "VETO" else "🟡"
                )

                with st.expander(
                    f"{_verdict_icon} [{_created}] {_orig_id} → {_ref_id} "
                    f"({_orig_v} → {_ref_v}) · {_proto}",
                    expanded=(_i == 0)  # expand most recent
                ):
                    _dc1, _dc2, _dc3 = st.columns(3)

                    # Verdict delta
                    _dc1.metric(
                        "Verdict",
                        _ref_v,
                        delta=f"{_orig_v} → {_ref_v}",
                        delta_color="normal"
                    )

                    # Score delta
                    if isinstance(_orig_s, (int, float)) and isinstance(_ref_s, (int, float)):
                        _dc2.metric(
                            "Det. Score",
                            f"{_ref_s}/100",
                            delta=f"{_ref_s - _orig_s:+d} pts"
                        )
                    else:
                        _dc2.metric("Det. Score", str(_ref_s))

                    # Coverage delta
                    if isinstance(_orig_cov, (int, float)) and isinstance(_ref_cov, (int, float)):
                        _dc3.metric(
                            "Coverage",
                            f"{_ref_cov:.1f}%",
                            delta=f"{_ref_cov - _orig_cov:+.1f}%"
                        )
                    else:
                        _dc3.metric("Coverage", str(_ref_cov))

                    # Metadata row
                    st.caption(
                        f"Protocol: `{_proto}` · "
                        f"Method: `{_method}` · "
                        f"Original: `{_orig_id}` · "
                        f"Refined: `{_ref_id}`"
                    )

                    # Load comparison button
                    if st.button(
                        "📊 Load Comparison",
                        key=f"hist_cmp_{_i}_{_orig_id}",
                        help="Load this comparison into the Refinement tab"
                    ):
                        st.session_state["comparison_baseline"] = _orig_id
                        st.session_state["comparison_refined"]  = _ref_id
                        st.success(
                            "Comparison loaded. Switch to the **Refinement** tab "
                            "to see the full before/after analysis."
                        )
    else:
        st.info("Enter your email or IP above and click **Load History** to see your refinement journey.")

    # ── Retention metrics (admin view) ────────────────────────────────────────
    with st.expander("📈 Platform Metrics (Admin)", expanded=False):
        st.caption("Aggregate usage across all users. Metadata only — no coordinates stored.")
        try:
            _cost_resp = api("GET", "/costs/summary")
            if _cost_resp:
                _cm1, _cm2, _cm3 = st.columns(3)
                _cm1.metric("Total Jobs",       _cost_resp.get("total_jobs", "—"))
                _cm2.metric("Total GPU Minutes", f"{_cost_resp.get('total_gpu_minutes', 0):.1f}")
                _cm3.metric("Total Cost (USD)",  f"${_cost_resp.get('total_usd', 0):.4f}")
            else:
                st.caption("Cost summary endpoint not available.")
        except Exception as _me:
            st.caption("Metrics unavailable.")
