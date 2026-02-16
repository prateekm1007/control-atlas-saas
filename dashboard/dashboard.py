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

# ═══════════════════════════════════════════════════════════════════════════════
# BENCHMARK DATASET — Canonical experimental structures
# Must PASS 12/12 deterministic or engine is miscalibrated
# ═══════════════════════════════════════════════════════════════════════════════
BENCHMARK_CRYSTALS = [
    {
        "pdb_id": "4HHB",
        "name": "Hemoglobin (Deoxy)",
        "resolution": 1.74,
        "method": "X-ray",
        "uniprot": "P69905",
        "af_id": "AF-P69905-F1",
    },
    {
        "pdb_id": "1CRN",
        "name": "Crambin",
        "resolution": 0.945,
        "method": "X-ray",
        "uniprot": "P01542",
        "af_id": "AF-P01542-F1",
    },
    {
        "pdb_id": "2PTC",
        "name": "Trypsin (Bovine)",
        "resolution": 1.55,
        "method": "X-ray",
        "uniprot": "P00760",
        "af_id": "AF-P00760-F1",
    },
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
    """PIL-VIS-01: Theme-aware 3D viewer. Works for single or side-by-side."""
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
    """Run a single benchmark audit through the standard /ingest pipeline."""
    return api(
        "POST",
        "/ingest",
        data={
            "mode": "Benchmark",
            "candidate_id": candidate_id,
            "t3_category": "NONE",
        },
    )


def render_law_comparison_table(crystal_laws, af_laws):
    """Deterministic-only side-by-side comparison table."""
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
        if c and c.get("method") != "deterministic":
            continue
        if a and a.get("method") != "deterministic":
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
            return f'<span style="color:{color};font-weight:bold;">{val} {units}</span> ({s})'

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


# ── Session State ───────────────────────────────────────────
if "audit_result" not in st.session_state:
    st.session_state.audit_result = None
if "candidates" not in st.session_state:
    st.session_state.candidates = []
if "benchmark_results" not in st.session_state:
    st.session_state.benchmark_results = {}

# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ═══════════════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/1042/1042307.png", width=80)
    st.title("STATION CONFIG")

    # ── Compact ☀️ / 🌙 toggle ─────────────────────────────
    st.markdown(
        '<div style="display:flex;align-items:center;gap:6px;margin:4px 0 8px 0;">'
        '<span style="font-size:1.1rem;">☀️</span>'
        '</div>',
        unsafe_allow_html=True,
    )
    st.toggle("🌙", value=True, key="noir_toggle", label_visibility="visible")

    st.divider()
    warhead = st.radio(
        "🎯 ARCHITECTURE WARHEAD",
        WARHEADS,
        index=0,
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

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN TABS
# ═══════════════════════════════════════════════════════════════════════════════
t_work, t_bench, t_valid, t_nkg, t_laws = st.tabs(
    ["⚡ Resolve", "🔬 Benchmark", "📊 Validation Dataset", "🧠 Knowledge Graph", "📜 Law Canon"]
)

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 1: RESOLVE (existing Command Center)
# ═══════════════════════════════════════════════════════════════════════════════
with t_work:
    if not st.session_state.audit_result:
        st.subheader("🔍 Structure Discovery")
        q = st.text_input(
            "UNIPROT / ALPHAFOLD RESOLUTION",
            placeholder="e.g. insulin, EGFR, p53...",
        )

        if st.button(
            "EXECUTE EXHAUSTIVE RESOLUTION",
            type="primary",
            use_container_width=True,
        ):
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
                            res = api(
                                "POST",
                                "/ingest",
                                data={
                                    "mode": "Discovery",
                                    "candidate_id": c["id"],
                                    "t3_category": warhead,
                                },
                            )
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

        col_narr_a, col_narr_b = st.columns(2)
        with col_narr_a:
            with st.expander("🔬 Deep Dive Analysis", expanded=True):
                st.write(witness.get("deep_dive", "Analysis pending."))
        with col_narr_b:
            with st.expander("💡 Governance Recommendation", expanded=True):
                st.write(witness.get("recommendation", "Awaiting human review."))

        st.caption(
            f"Narrative Model: {res.get('ai_model_used', 'Internal Fallback')}"
        )

        st.divider()
        st.subheader("🧬 Forensic 3D Reconstruction")
        if res.get("pdb_b64"):
            render_3d_viewer(res["pdb_b64"], binary, height=600, width=800)

        st.divider()
        st.subheader("📋 Tier-1 Diagnostic Ledger")
        laws = res.get("tier1", {}).get("laws", [])
        for law in laws:
            icon = "✅" if law["status"] == "PASS" else "🛑"
            with st.expander(
                f"{icon} {law['law_id']}: {law['title']} — {law['status']}"
            ):
                st.write(f"**Observed:** {law['observed']} {law['units']}")
                st.write(
                    f"**Threshold:** {law['operator']} {law['threshold']} {law['units']}"
                )
                st.write(f"**Deviation:** {law['deviation']}")
                st.write(
                    f"**Evaluation Scope:** {law['sample_size']} {law['scope']}"
                )
                st.write(f"**Scientific Principle:** {law['principle']}")

        st.divider()
        if res.get("pdf_b64"):
            pdf_bytes = base64.b64decode(res["pdf_b64"])
            st.download_button(
                label="📄 DOWNLOAD FORENSIC DOSSIER (PDF)",
                data=pdf_bytes,
                file_name=(
                    f"TOSCANINI_DOSSIER_"
                    f"{res.get('governance', {}).get('audit_id', 'ARTIFACT')}.pdf"
                ),
                mime="application/pdf",
                use_container_width=True,
            )
        else:
            st.error("PIL-PAR-14 VIOLATION: PDF byte-stream not found in payload.")

        if st.button("🔄 RETURN TO COMMAND CENTER", use_container_width=True):
            st.session_state.audit_result = None
            st.rerun()

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 2: BENCHMARK — Side-by-Side Crystal vs AlphaFold
# ═══════════════════════════════════════════════════════════════════════════════
with t_bench:
    st.subheader("🔬 Benchmark Mode — Experimental vs Predicted")
    st.markdown(
        "Each structure is audited through the **identical** deterministic pipeline. "
        "Crystal structures **must** PASS 12/12 or the engine is miscalibrated."
    )
    st.divider()

    for entry in BENCHMARK_CRYSTALS:
        pdb_id = entry["pdb_id"]
        af_id = entry["af_id"]
        key = pdb_id

        with st.expander(
            f"**{entry['name']}** — {pdb_id} ({entry['resolution']} Å, {entry['method']})",
            expanded=False,
        ):
            cached = st.session_state.benchmark_results.get(key, {})

            if not cached.get("crystal") or not cached.get("alphafold"):
                if st.button(f"⚡ Run Comparative Audit: {pdb_id} vs {af_id}", key=f"bench_{key}"):
                    with st.spinner(f"Auditing {pdb_id} (Crystal)..."):
                        crystal_res = run_benchmark_audit(pdb_id)
                    with st.spinner(f"Auditing {af_id} (AlphaFold)..."):
                        af_res = run_benchmark_audit(af_id)

                    st.session_state.benchmark_results[key] = {
                        "crystal": crystal_res,
                        "alphafold": af_res,
                    }
                    st.rerun()
            else:
                cr = cached["crystal"]
                ar = cached["alphafold"]

                # ── Executive Comparison Strip ──────────────
                cv = cr.get("verdict", {})
                av = ar.get("verdict", {})

                c1, c2 = st.columns(2)
                with c1:
                    st.markdown(f"**🔵 Crystal: {pdb_id}**")
                    binary_c = cv.get("binary", "ERROR")
                    if binary_c == "PASS":
                        st.success(f"✅ {binary_c} — {cv.get('det_passed', '?')}/{cv.get('det_total', '?')} Deterministic")
                    else:
                        st.error(f"🛑 {binary_c}")
                    st.metric("Deterministic", f"{cv.get('deterministic_score', 0)}%")
                    st.metric("Coverage", f"{cv.get('coverage_pct', 0)}%")

                with c2:
                    st.markdown(f"**🟠 AlphaFold: {af_id}**")
                    binary_a = av.get("binary", "ERROR")
                    if binary_a == "PASS":
                        st.success(f"✅ {binary_a} — {av.get('det_passed', '?')}/{av.get('det_total', '?')} Deterministic")
                    elif binary_a == "VETO":
                        st.error(f"🛑 {binary_a}")
                    else:
                        st.warning(f"⚠️ {binary_a}")
                    st.metric("Deterministic", f"{av.get('deterministic_score', 0)}%")
                    st.metric("pLDDT", f"{round(av.get('confidence_score', 0), 1)}")

                # ── Side-by-Side 3D ─────────────────────────
                st.divider()
                st.markdown("**3D Structural Comparison**")
                v3d_l, v3d_r = st.columns(2)
                with v3d_l:
                    st.caption(f"Crystal: {pdb_id}")
                    if cr.get("pdb_b64"):
                        render_3d_viewer(cr["pdb_b64"], cv.get("binary", "ERROR"), height=400, width=420)
                with v3d_r:
                    st.caption(f"AlphaFold: {af_id}")
                    if ar.get("pdb_b64"):
                        render_3d_viewer(ar["pdb_b64"], av.get("binary", "ERROR"), height=400, width=420)

                # ── Quantitative Comparison Table ───────────
                st.divider()
                st.markdown("**Deterministic Law Comparison**")
                c_laws = cr.get("tier1", {}).get("laws", [])
                a_laws = ar.get("tier1", {}).get("laws", [])
                render_law_comparison_table(c_laws, a_laws)

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 3: VALIDATION DATASET — Statistical Trust
# ═══════════════════════════════════════════════════════════════════════════════
with t_valid:
    st.subheader("📊 Validation Dataset — Statistical Trust Summary")
    st.markdown(
        "Aggregate calibration metrics across all benchmarked experimental structures. "
        "A false veto rate >0% against high-resolution crystals indicates engine miscalibration."
    )
    st.divider()

    bench = st.session_state.benchmark_results
    if not bench:
        st.info(
            "No benchmark data yet. Run comparative audits in the **Benchmark** tab first."
        )
    else:
        n_tested = len(bench)
        crystal_pass = 0
        crystal_total = 0
        all_clash = []
        all_bond = []
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

            c_laws = cr.get("tier1", {}).get("laws", [])
            for law in c_laws:
                if law["law_id"] == "LAW-130":
                    try:
                        all_clash.append(float(law["observed"]))
                    except (ValueError, TypeError):
                        pass
                if law["law_id"] == "LAW-100":
                    try:
                        all_bond.append(float(law["observed"]))
                    except (ValueError, TypeError):
                        pass

        false_veto_rate = (
            round((crystal_total - crystal_pass) / crystal_total * 100, 1)
            if crystal_total > 0
            else 0.0
        )

        s1, s2, s3, s4 = st.columns(4)
        s1.metric("Structures Tested", n_tested)
        s2.metric("Crystal Pass Rate", f"{round(crystal_pass / max(crystal_total, 1) * 100, 1)}%")
        s3.metric("False Veto Rate", f"{false_veto_rate}%")
        s4.metric(
            "Resolution Range",
            f"{min(resolution_range, default=0):.2f}–{max(resolution_range, default=0):.2f} Å",
        )

        if all_clash:
            st.divider()
            import numpy as np

            mc1, mc2 = st.columns(2)
            mc1.metric("Mean Clashscore (Crystal)", f"{round(float(np.mean(all_clash)), 2)} clashes/1000 atoms")
            if all_bond:
                mc2.metric("Mean Bond Outliers (Crystal)", f"{round(float(np.mean(all_bond)), 2)}")

        st.divider()
        if false_veto_rate == 0.0 and crystal_total > 0:
            st.success(
                f"✅ Toscanini deterministic invariants were benchmarked against {crystal_total} "
                f"high-resolution experimental crystal structures. "
                f"**No false veto events observed** within validation dataset (n={crystal_total})."
            )
        elif crystal_total > 0:
            st.error(
                f"🛑 CALIBRATION ALERT: {crystal_total - crystal_pass} crystal structure(s) "
                f"received false vetoes. Engine requires threshold review."
            )

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 4: KNOWLEDGE GRAPH
# ═══════════════════════════════════════════════════════════════════════════════
with t_nkg:
    st.subheader("🧠 Knowledge Graph (NKG)")
    nkg_data = api("GET", "/nkg")
    if nkg_data:
        st.metric("Total System Vetoes", nkg_data.get("total_vetoes", 0))
        if nkg_data.get("vetoes"):
            st.json(nkg_data["vetoes"])
    else:
        st.info("Knowledge Graph initialized. Moat is empty.")

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 5: LAW CANON
# ═══════════════════════════════════════════════════════════════════════════════
with t_laws:
    st.subheader("📜 Structural Law Canon")
    laws_data = api("GET", "/laws")
    if laws_data:
        for lid, info in laws_data.get("laws", {}).items():
            with st.expander(f"{lid}: {info['title']}"):
                st.write(f"**Principle:** {info.get('principle', 'N/A')}")
                st.write(f"**Metric Type:** {info.get('type')}")
                st.write(
                    f"**Threshold:** {info.get('operator')} "
                    f"{info.get('threshold')} {info.get('unit')}"
                )
