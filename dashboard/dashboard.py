import streamlit as st
import requests, base64, json, os
from style_utils import (
    apply_arena_theme, render_lifecycle_header, score_color,
    THEMES, get_active_theme, get_active_theme_name,
)

st.set_page_config(page_title="TOSCANINI OS", layout="wide", page_icon="🛡️")

if "noir_toggle" in st.session_state:
    st.session_state.toscanini_theme = (
        "Biotech Noir" if st.session_state.noir_toggle else "Clinical White"
    )

apply_arena_theme()

BACKEND = os.getenv("BACKEND_URL", "http://brain:8000")
API_KEY = os.getenv("TOSCANINI_API_KEY", "")

WARHEADS = [
    "NONE", "LINEAR", "MULTIVALENT", "ENGAGEMENT",
    "MOTIFS", "LINKERS", "CONFORMATIONAL", "METAL",
]

BENCHMARK_CRYSTALS = [
    {"pdb_id": "3NIR", "name": "Endothiapepsin (Ultra-High)", "resolution": 0.48, "method": "X-ray", "af_id": None},
    {"pdb_id": "2VB1", "name": "Insulin (Atomic)", "resolution": 0.65, "method": "X-ray", "af_id": "AF-P01308-F1"},
    {"pdb_id": "1CRN", "name": "Crambin", "resolution": 1.50, "method": "X-ray", "af_id": "AF-P01542-F1"},
    {"pdb_id": "4HHB", "name": "Hemoglobin (Deoxy)", "resolution": 1.74, "method": "X-ray", "af_id": "AF-P69905-F1"},
    {"pdb_id": "1UBQ", "name": "Ubiquitin", "resolution": 1.80, "method": "X-ray", "af_id": "AF-P0CG47-F1"},
    {"pdb_id": "6LZG", "name": "SARS-CoV-2 Spike RBD", "resolution": 2.50, "method": "X-ray", "af_id": "AF-P0DTC2-F1"},
]


def api(method, path, **kwargs):
    headers = kwargs.pop("headers", {})
    if API_KEY:
        headers["X-API-Key"] = API_KEY
    try:
        r = requests.request(
            method, f"{BACKEND}{path}", timeout=180, headers=headers, **kwargs
        )
        r.raise_for_status()
        return r.json()
    except Exception as e:
        st.error(f"Backend Link Error ({path}): {e}")
        return None


def render_3d_viewer(pdb_b64, verdict, height=500, width=None):
    try:
        pdb_str = base64.decodebytes(pdb_b64.encode()).decode("utf-8", errors="ignore")
        if len(pdb_str.strip()) < 50:
            return
        import py3Dmol
        import streamlit.components.v1 as components

        t = get_active_theme()
        w = width or 700
        view = py3Dmol.view(width=w, height=height)
        view.addModel(pdb_str, "pdb")
        if verdict == "PASS":
            view.setStyle({"cartoon": {"color": "spectrum"}})
        else:
            view.setStyle(
                {"stick": {"colorscheme": "orangeCarbon"}, "sphere": {"radius": 0.5}}
            )
        view.setBackgroundColor(t["viz_bg"])
        view.zoomTo()
        components.html(view._make_html(), height=height + 20, width=w + 20, scrolling=False)
    except Exception as e:
        st.warning(f"Visualizer Unavailable: {e}")


def run_benchmark_audit(candidate_id):
    return api(
        "POST", "/ingest",
        data={"mode": "Benchmark", "candidate_id": candidate_id, "t3_category": "NONE"},
    )


def render_law_comparison_table(crystal_laws, af_laws):
    t = get_active_theme()
    header = (
        f'<tr style="background:{t["surface"]};border-bottom:2px solid {t["border"]};">'
        f'<th style="padding:8px;color:{t["text"]};text-align:left;">Law</th>'
        f'<th style="padding:8px;color:{t["text"]};text-align:left;">Title</th>'
        f'<th style="padding:8px;color:{t["text"]};text-align:center;">Crystal</th>'
        f'<th style="padding:8px;color:{t["text"]};text-align:center;">AlphaFold</th>'
        f'<th style="padding:8px;color:{t["text"]};text-align:right;">Threshold</th>'
        f'</tr>'
    )
    c_map = {l["law_id"]: l for l in crystal_laws} if crystal_laws else {}
    a_map = {l["law_id"]: l for l in af_laws} if af_laws else {}
    all_ids = sorted(set(list(c_map.keys()) + list(a_map.keys())))
    rows = ""
    for lid in all_ids:
        c = c_map.get(lid)
        a = a_map.get(lid)
        if c and c.get("method") not in ("deterministic", "advisory_experimental"):
            continue
        if a and a.get("method") not in ("deterministic", "advisory_experimental"):
            continue
        title = (c or a or {}).get("title", lid)
        thresh = (c or a or {}).get("threshold", "—")
        op = (c or a or {}).get("operator", "")
        units = (c or a or {}).get("units", "")

        def _cell(law_row):
            if not law_row:
                return f'<span style="color:{t["muted"]};">—</span>'
            s = law_row["status"]
            val = law_row["observed"]
            color = "#00CC66" if s == "PASS" else "#FF4444"
            method_badge = ""
            if law_row.get("method") == "advisory_experimental":
                method_badge = f' <span style="font-size:0.65rem;color:#D4A017;background:#3a3000;padding:1px 4px;border-radius:3px;">ADV</span>'
            return f'<span style="color:{color};font-weight:bold;">{val} {units}</span> ({s}){method_badge}'

        rows += (
            f'<tr style="border-bottom:1px solid {t["border"]};">'
            f'<td style="padding:6px;color:{t["text"]};font-family:monospace;">{lid}</td>'
            f'<td style="padding:6px;color:{t["text"]};">{title}</td>'
            f'<td style="padding:6px;text-align:center;">{_cell(c)}</td>'
            f'<td style="padding:6px;text-align:center;">{_cell(a)}</td>'
            f'<td style="padding:6px;text-align:right;color:{t["muted"]};font-family:monospace;">{op} {thresh} {units}</td>'
            f'</tr>'
        )
    st.markdown(
        f'<table style="width:100%;border-collapse:collapse;margin-top:10px;">{header}{rows}</table>',
        unsafe_allow_html=True,
    )


def render_law_badge(law):
    """Render a single law with method-aware badge (Task 1)."""
    method = law.get("method", "unknown")
    status = law["status"]

    if status == "PASS":
        icon = "✅"
    else:
        icon = "🛑"

    if method == "advisory_experimental":
        badge = "🟡 Advisory (Experimental)"
        badge_help = "Reported for transparency. Excluded from deterministic score."
    elif method == "deterministic":
        badge = "🟢 Deterministic"
        badge_help = "Counts toward deterministic governance score."
    elif method == "heuristic":
        badge = "🔵 Heuristic"
        badge_help = "Statistical proxy. Does not trigger VETO."
    else:
        badge = f"⚪ {method}"
        badge_help = ""

    return icon, badge, badge_help


# ── Session State ───────────────────────────────────────────
if "audit_result" not in st.session_state:
    st.session_state.audit_result = None
if "candidates" not in st.session_state:
    st.session_state.candidates = []
if "benchmark_results" not in st.session_state:
    st.session_state.benchmark_results = {}

# ═══════════════════════════════════════════════════════════════
# SIDEBAR
# ═══════════════════════════════════════════════════════════════
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/1042/1042307.png", width=80)
    st.title("STATION CONFIG")
    st.markdown(
        '<div style="display:flex;align-items:center;gap:6px;margin:4px 0 8px 0;">'
        '<span style="font-size:1.1rem;">☀️</span>'
        '</div>',
        unsafe_allow_html=True,
    )
    st.toggle("🌙", value=True, key="noir_toggle", label_visibility="visible")
    st.divider()
    warhead = st.radio(
        "🎯 ARCHITECTURE WARHEAD", WARHEADS, index=0,
        help="Selects therapeutic class weighting.",
    )
    st.divider()
    st.caption(f"Engine: v22.5.x | Theme: {get_active_theme_name()}")
    st.caption("Status: ONLINE")

# ── Header ──────────────────────────────────────────────────
t = get_active_theme()
st.markdown(
    f'<h1 style="color:{t["accent"]};text-align:center;font-family:Courier New;'
    f'margin-bottom:0px;">🛡️ TOSCANINI</h1>',
    unsafe_allow_html=True,
)
current_state = "INSTANTIATED"
if st.session_state.audit_result:
    current_state = "AUDITED"
elif st.session_state.candidates:
    current_state = "ACQUIRING"
render_lifecycle_header(current_state)

t_work, t_bench, t_valid, t_nkg, t_laws = st.tabs(
    ["⚡ Resolve", "🔬 Benchmark", "📊 Validation Dataset", "🧠 Knowledge Graph", "📜 Law Canon"]
)

# ═══════════════════════════════════════════════════════════════
# TAB 1: RESOLVE
# ═══════════════════════════════════════════════════════════════
with t_work:
    if not st.session_state.audit_result:
        st.subheader("🔍 Structure Discovery")
        q = st.text_input("UNIPROT / ALPHAFOLD RESOLUTION", placeholder="e.g. insulin, EGFR, p53...")
        if st.button("EXECUTE EXHAUSTIVE RESOLUTION", type="primary", use_container_width=True):
            if q.strip():
                with st.spinner("Searching Global Repositories..."):
                    resp = api("POST", "/search", data={"query": q})
                    if resp and "results" in resp:
                        st.session_state.candidates = resp["results"]
                    else:
                        st.warning("Resolution returned 0 isoforms.")
        if st.session_state.candidates:
            st.divider()
            for i, c in enumerate(st.session_state.candidates):
                col_label, col_btn = st.columns([5, 1])
                with col_label:
                    st.markdown(f"**{c.get('id')}** - {c.get('label')}")
                with col_btn:
                    if st.button("⚡ AUDIT", key=f"audit_{i}"):
                        with st.spinner("Executing Refusal Engine..."):
                            res = api("POST", "/ingest", data={"mode": "Discovery", "candidate_id": c["id"], "t3_category": warhead})
                            if res:
                                st.session_state.audit_result = res
                                st.rerun()
    else:
        res = st.session_state.audit_result
        v = res.get("verdict", {})
        binary = v.get("binary", "ERROR")
        if binary == "PASS":
            st.success("✅ STRUCTURE PASSED: ALL DETERMINISTIC INVARIANTS SATISFIED")
        elif binary == "VETO":
            st.error(f"🛑 DETERMINISTIC VETO: {', '.join(v.get('killer_laws', []))}")
        else:
            st.warning(f"⚠️ VERDICT: {binary}")

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("DETERMINISTIC", f"{v.get('deterministic_score', 0)}%")
        m2.metric("ADVISORY", f"{v.get('advisory_score', 0)}%")
        m3.metric("pLDDT (CORE)", f"{round(v.get('confidence_score', 0), 1)}")
        m4.metric("EPI PRIORITY", f"{res.get('tier3', {}).get('probability', 0)}%")

        st.divider()
        st.subheader("🤖 PhD Witness Report")
        witness = res.get("witness_reports", {})
        st.info(witness.get("executive", "Executive Summary Unreachable."))
        col_a, col_b = st.columns(2)
        with col_a:
            with st.expander("🔬 Deep Dive Analysis", expanded=True):
                st.write(witness.get("deep_dive", "Analysis pending."))
        with col_b:
            with st.expander("💡 Governance Recommendation", expanded=True):
                st.write(witness.get("recommendation", "Awaiting human review."))
        st.caption(f"Narrative Model: {res.get('ai_model_used', 'Internal Fallback')}")

        st.divider()
        st.subheader("🧬 Forensic 3D Reconstruction")
        if res.get("pdb_b64"):
            render_3d_viewer(res["pdb_b64"], binary, height=600, width=800)

        # ── TASK 1: Advisory-aware Diagnostic Ledger ────────────
        st.divider()
        st.subheader("📋 Tier-1 Diagnostic Ledger")
        st.caption("🟢 Deterministic = counts toward score | 🟡 Advisory = reported only, not scored | 🔵 Heuristic = statistical proxy")

        laws = res.get("tier1", {}).get("laws", [])
        for law in laws:
            icon, badge, badge_help = render_law_badge(law)
            with st.expander(f"{icon} {law['law_id']}: {law['title']} — {law['status']}  [{badge}]"):
                if badge_help:
                    st.caption(f"ℹ️ {badge_help}")
                st.write(f"**Observed:** {law['observed']} {law['units']}")
                st.write(f"**Threshold:** {law['operator']} {law['threshold']} {law['units']}")
                st.write(f"**Deviation:** {law['deviation']}")
                st.write(f"**Evaluation Scope:** {law['sample_size']} {law['scope']}")
                st.write(f"**Classification:** {badge}")

        st.divider()
        if res.get("pdf_b64"):
            pdf_bytes = base64.b64decode(res["pdf_b64"])
            st.download_button(
                label="📄 DOWNLOAD FORENSIC DOSSIER (PDF)", data=pdf_bytes,
                file_name=f"TOSCANINI_DOSSIER_{res.get('governance', {}).get('audit_id', 'ARTIFACT')}.pdf",
                mime="application/pdf", use_container_width=True,
            )
        else:
            st.error("PIL-PAR-14 VIOLATION: PDF byte-stream not found in payload.")

        if st.button("🔄 RETURN TO COMMAND CENTER", use_container_width=True):
            st.session_state.audit_result = None
            st.rerun()

# ═══════════════════════════════════════════════════════════════
# TAB 2: BENCHMARK (Task 3 — expanded dataset)
# ═══════════════════════════════════════════════════════════════
with t_bench:
    st.subheader("🔬 Benchmark Mode — Experimental vs Predicted")
    st.markdown(
        "Each structure is audited through the **identical** deterministic pipeline. "
        "Crystal structures **must** PASS to prove the engine is calibrated to biological reality."
    )
    st.divider()

    for entry in BENCHMARK_CRYSTALS:
        pdb_id = entry["pdb_id"]
        af_id = entry.get("af_id")
        key = pdb_id

        with st.expander(
            f"**{entry['name']}** — {pdb_id} ({entry['resolution']} Å, {entry['method']})",
            expanded=False,
        ):
            cached = st.session_state.benchmark_results.get(key, {})
            has_crystal = cached.get("crystal") is not None
            has_af = cached.get("alphafold") is not None

            if not has_crystal:
                if st.button(f"⚡ Audit {pdb_id}" + (f" + {af_id}" if af_id else ""), key=f"bench_{key}"):
                    with st.spinner(f"Auditing {pdb_id}..."):
                        crystal_res = run_benchmark_audit(pdb_id)
                    af_res = None
                    if af_id:
                        with st.spinner(f"Auditing {af_id}..."):
                            af_res = run_benchmark_audit(af_id)
                    st.session_state.benchmark_results[key] = {"crystal": crystal_res, "alphafold": af_res}
                    st.rerun()
            else:
                cr = cached["crystal"]
                ar = cached.get("alphafold")
                cv = cr.get("verdict", {})

                c1, c2 = st.columns(2) if ar else (st.container(), None)
                with c1:
                    st.markdown(f"**🔵 Crystal: {pdb_id}**")
                    binary_c = cv.get("binary", "ERROR")
                    det_laws_c = [l for l in cr.get("tier1", {}).get("laws", []) if l["method"] == "deterministic"]
                    det_pass_c = sum(1 for l in det_laws_c if l["status"] == "PASS")
                    if binary_c == "PASS":
                        st.success(f"✅ {binary_c} — {det_pass_c}/{len(det_laws_c)} Det + 1 Advisory")
                    else:
                        st.error(f"🛑 {binary_c}")
                    st.metric("Coverage", f"{cv.get('coverage_pct', 0)}%")

                if c2 and ar:
                    av = ar.get("verdict", {})
                    with c2:
                        st.markdown(f"**🟠 AlphaFold: {af_id}**")
                        binary_a = av.get("binary", "ERROR")
                        if binary_a == "PASS":
                            st.success(f"✅ {binary_a} — {av.get('det_passed', '?')}/{av.get('det_total', '?')} Det")
                        elif binary_a == "VETO":
                            st.error(f"🛑 {binary_a}")
                        else:
                            st.warning(f"⚠️ {binary_a}")
                        st.metric("pLDDT", f"{round(av.get('confidence_score', 0), 1)}")

                # Side-by-side 3D
                if cr.get("pdb_b64"):
                    st.divider()
                    if ar and ar.get("pdb_b64"):
                        st.markdown("**3D Structural Comparison**")
                        v3d_l, v3d_r = st.columns(2)
                        with v3d_l:
                            st.caption(f"Crystal: {pdb_id}")
                            render_3d_viewer(cr["pdb_b64"], cv.get("binary", "ERROR"), height=400, width=420)
                        with v3d_r:
                            st.caption(f"AlphaFold: {af_id}")
                            render_3d_viewer(ar["pdb_b64"], av.get("binary", "ERROR"), height=400, width=420)
                    else:
                        st.caption(f"Crystal: {pdb_id}")
                        render_3d_viewer(cr["pdb_b64"], cv.get("binary", "ERROR"), height=400, width=600)

                # Comparison table
                if ar:
                    st.divider()
                    st.markdown("**Deterministic Law Comparison**")
                    render_law_comparison_table(
                        cr.get("tier1", {}).get("laws", []),
                        ar.get("tier1", {}).get("laws", []),
                    )

# ═══════════════════════════════════════════════════════════════
# TAB 3: VALIDATION DATASET
# ═══════════════════════════════════════════════════════════════
with t_valid:
    st.subheader("📊 Validation Dataset — Statistical Trust Summary")
    st.markdown(
        "Aggregate calibration metrics across all benchmarked experimental structures. "
        "A false veto rate >0% against high-resolution crystals indicates engine miscalibration."
    )
    st.divider()

    bench = st.session_state.benchmark_results
    if not bench:
        st.info("No benchmark data yet. Run comparative audits in the **Benchmark** tab first.")
    else:
        n_tested = len(bench)
        crystal_pass = 0
        crystal_total = 0
        all_clash = []
        resolution_range = []

        for key, data in bench.items():
            cr = data.get("crystal")
            if not cr:
                continue
            crystal_total += 1
            cv = cr.get("verdict", {})
            if cv.get("binary") == "PASS":
                crystal_pass += 1
            for entry in BENCHMARK_CRYSTALS:
                if entry["pdb_id"] == key:
                    resolution_range.append(entry["resolution"])
                    break
            for law in cr.get("tier1", {}).get("laws", []):
                if law["law_id"] == "LAW-130":
                    try:
                        all_clash.append(float(law["observed"]))
                    except (ValueError, TypeError):
                        pass

        false_veto_rate = round((crystal_total - crystal_pass) / crystal_total * 100, 1) if crystal_total > 0 else 0.0

        s1, s2, s3, s4 = st.columns(4)
        s1.metric("Structures Tested", n_tested)
        s2.metric("Crystal Pass Rate", f"{round(crystal_pass / max(crystal_total, 1) * 100, 1)}%")
        s3.metric("False Veto Rate", f"{false_veto_rate}%")
        s4.metric("Resolution Range", f"{min(resolution_range, default=0):.2f}–{max(resolution_range, default=0):.2f} Å")

        if all_clash:
            import numpy as np
            st.divider()
            st.metric("Mean Clashscore (Crystal)", f"{round(float(np.mean(all_clash)), 2)} clashes/1000 atoms")

        st.divider()
        if false_veto_rate == 0.0 and crystal_total > 0:
            st.success(
                f"✅ Toscanini deterministic invariants benchmarked against {crystal_total} "
                f"experimental crystal structures (resolution: {min(resolution_range):.2f}–{max(resolution_range):.2f} Å). "
                f"**No false veto events observed** (n={crystal_total})."
            )
        elif crystal_total > 0:
            st.error(f"🛑 CALIBRATION ALERT: {crystal_total - crystal_pass} false veto(es).")

with t_nkg:
    st.subheader("🧠 Knowledge Graph (NKG)")
    nkg_data = api("GET", "/nkg")
    if nkg_data:
        st.metric("Total System Vetoes", nkg_data.get("total_vetoes", 0))
        if nkg_data.get("vetoes"):
            st.json(nkg_data["vetoes"])
    else:
        st.info("Knowledge Graph initialized. Moat is empty.")

with t_laws:
    st.subheader("📜 Structural Law Canon")
    laws_data = api("GET", "/laws")
    if laws_data:
        for lid, info in laws_data.get("laws", {}).items():
            with st.expander(f"{lid}: {info['title']}"):
                st.write(f"**Principle:** {info.get('principle', 'N/A')}")
                st.write(f"**Metric Type:** {info.get('type')}")
                st.write(f"**Threshold:** {info.get('operator')} {info.get('threshold')} {info.get('unit')}")
