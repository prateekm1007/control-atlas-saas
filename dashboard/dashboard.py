import streamlit as st
import requests, base64, json, os
from style_utils import apply_arena_theme, render_lifecycle_header, score_color

st.set_page_config(page_title="TOSCANINI OS", layout="wide", page_icon="🛡️")
apply_arena_theme()

BACKEND = os.getenv("BACKEND_URL", "http://brain:8000")
API_KEY = os.getenv("TOSCANINI_API_KEY", "")
WARHEADS = ["NONE", "LINEAR", "MULTIVALENT", "ENGAGEMENT", "MOTIFS", "LINKERS", "CONFORMATIONAL", "METAL"]


def api(method, path, **kwargs):
    headers = kwargs.pop("headers", {})
    if API_KEY:
        headers["X-API-Key"] = API_KEY
    try:
        r = requests.request(method, f"{BACKEND}{path}", timeout=180, headers=headers, **kwargs)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        st.error(f"Backend Error: {e}")
        return None


def render_3d_viewer(pdb_b64):
    """Render 3D protein structure using py3Dmol HTML embed."""
    try:
        pdb_str = base64.b64decode(pdb_b64).decode("utf-8", errors="ignore")
        if len(pdb_str.strip()) < 50:
            st.warning("Structure too small for 3D rendering")
            return
        import py3Dmol
        import streamlit.components.v1 as components
        view = py3Dmol.view(width=700, height=500)
        view.addModel(pdb_str, "pdb")
        view.setStyle({"cartoon": {"color": "spectrum"}})
        view.setBackgroundColor("#0a0a0a")
        view.zoomTo()
        html_str = view._make_html()
        components.html(html_str, height=520, width=720, scrolling=False)
    except Exception as e:
        st.warning(f"3D rendering unavailable: {e}")

# ── Session State ──────────────────────────────────────────
for key, default in [
    ("audit_result", None),
    ("candidates", []),
    ("warhead", "NONE"),
]:
    if key not in st.session_state:
        st.session_state[key] = default

# ── Header ─────────────────────────────────────────────────
st.markdown(
    '<h1 style="color:#00D4FF;text-align:center;font-family:Courier New;'
    'margin-bottom:0px;">🛡️ TOSCANINI</h1>',
    unsafe_allow_html=True,
)
st.markdown(
    '<p style="text-align:center;color:#666;font-size:0.8rem;margin-top:0;">'
    'Structural Governance Station — Forensic Falsification Engine</p>',
    unsafe_allow_html=True,
)

if st.session_state.audit_result:
    current_state = "AUDITED"
elif st.session_state.candidates:
    current_state = "ACQUIRING"
else:
    current_state = "INSTANTIATED"
render_lifecycle_header(current_state)

# ── Tabs ───────────────────────────────────────────────────
t_work, t_nkg, t_laws = st.tabs(["⚡ Command Center", "🧠 Knowledge Graph", "📜 Law Canon"])

# ══════════════════════════════════════════════════════════
#  TAB 1: COMMAND CENTER
# ══════════════════════════════════════════════════════════
with t_work:

    # ──────────────────────────────────────────────────────
    #  PHASE 1: SEARCH + UPLOAD (no audit result yet)
    # ──────────────────────────────────────────────────────
    if st.session_state.audit_result is None:

        # ── Search Bar ──
        col_search, col_warhead = st.columns([3, 1])

        with col_search:
            q = st.text_input(
                "🔍 Search UniProt / AlphaFold DB",
                placeholder="e.g. insulin, EGFR, hemoglobin, P53, BRCA1...",
            )

        with col_warhead:
            warhead = st.selectbox(
                "🎯 Architecture",
                WARHEADS,
                index=WARHEADS.index(st.session_state.warhead),
                help="Architectural context for EPI prioritization",
            )
            st.session_state.warhead = warhead

        col_btn_search, col_btn_upload = st.columns([1, 1])

        with col_btn_search:
            search_clicked = st.button("🔎 Search Database", type="primary", use_container_width=True)

        with col_btn_upload:
            uploaded = st.file_uploader(
                "📁 Or upload PDB/CIF",
                type=["pdb", "cif"],
                label_visibility="collapsed",
            )

        # ── Execute Search ──
        if search_clicked and q.strip():
            with st.spinner(f"Searching AlphaFold DB for '{q}'..."):
                resp = api("POST", "/search", data={"query": q})
                if resp and isinstance(resp, dict):
                    st.session_state.candidates = resp.get("results", [])
                elif resp and isinstance(resp, list):
                    st.session_state.candidates = resp
                else:
                    st.session_state.candidates = []

        # ── Execute Upload Audit ──
        if uploaded is not None:
            if st.button("⚡ Audit Uploaded File", type="primary", use_container_width=True):
                with st.spinner(f"Auditing {uploaded.name}..."):
                    res = api(
                        "POST", "/ingest",
                        data={"mode": "Upload", "candidate_id": uploaded.name, "t3_category": warhead},
                        files={"file": (uploaded.name, uploaded.getvalue())},
                    )
                    if res:
                        st.session_state.audit_result = res
                        st.rerun()

        # ── Candidate Cards ──
        candidates = st.session_state.candidates
        if candidates:
            st.divider()
            st.markdown(
                f'<h3 style="color:#00D4FF;">📋 {len(candidates)} Candidates Found</h3>',
                unsafe_allow_html=True,
            )

            for i, c in enumerate(candidates):
                cid = c.get("id", "N/A")
                label = c.get("label", "Unknown")
                source = c.get("source", "")

                with st.container():
                    col_info, col_action = st.columns([5, 1])

                    with col_info:
                        st.markdown(
                            f'<div style="background:#111;border:1px solid #333;'
                            f'border-radius:8px;padding:10px 15px;margin:3px 0;">'
                            f'<span style="color:#00D4FF;font-weight:bold;font-size:0.9rem;">{cid}</span>'
                            f'<br/>'
                            f'<span style="color:#aaa;font-size:0.8rem;">{label}</span>'
                            f'</div>',
                            unsafe_allow_html=True,
                        )

                    with col_action:
                        if st.button(
                            f"⚡ Audit",
                            key=f"audit_{i}",
                            use_container_width=True,
                            type="primary",
                        ):
                            with st.spinner(f"Fetching & auditing {cid}..."):
                                res = api(
                                    "POST", "/ingest",
                                    data={
                                        "mode": "Discovery",
                                        "candidate_id": cid,
                                        "t3_category": st.session_state.warhead,
                                    },
                                )
                                if res:
                                    st.session_state.audit_result = res
                                    st.rerun()

        elif search_clicked:
            st.warning("No candidates found. Try a different query.")

    # ──────────────────────────────────────────────────────
    #  PHASE 2: AUDIT RESULTS
    # ──────────────────────────────────────────────────────
    else:
        res = st.session_state.audit_result
        v = res.get("verdict", {})
        binary = v.get("binary", "ERROR")
        prov = res.get("provenance", {})
        conf_meta = res.get("confidence_meta", {})

        # ── Verdict Banner ──
        if binary == "PASS":
            st.success("✅ STRUCTURE PASSED — All deterministic invariants satisfied")
        elif binary == "VETO":
            killers = v.get("killer_laws", [])
            msg = ", ".join(killers[:3]) if killers else "Deterministic violation"
            st.error(f"🛑 DETERMINISTIC VETO — {msg}")
        elif binary == "INDETERMINATE":
            st.warning("⚠️ INDETERMINATE — Insufficient coverage or fringe escalation")
        else:
            st.error(f"❌ ERROR: {binary}")

        # ── Source Info ──
        st.caption(
            f"Source: **{prov.get('source', 'N/A')}** | "
            f"Bytes: {prov.get('byte_count', 0):,} | "
            f"Version: {prov.get('station_version', '?')} | "
            f"SHA-256: `{prov.get('hash', 'N/A')[:16]}...`"
        )

        # ── 4-Metric Grid ──
        m1, m2, m3, m4 = st.columns(4)
        det_score = v.get("deterministic_score", 0)
        adv_score = v.get("advisory_score", 0)
        conf_score = v.get("confidence_score", 0)
        epi = res.get("tier3", {}).get("probability", 0)

        m1.metric(
            "DETERMINISTIC", f"{det_score}%",
            delta=f"{v.get('det_passed', 0)}/{v.get('det_total', 0)} laws",
            delta_color="normal" if det_score == 100 else "inverse",
        )
        m2.metric(
            "ADVISORY", f"{adv_score}%",
            delta=f"{v.get('heur_passed', 0)}/{v.get('heur_total', 0)} laws",
        )
        m3.metric(
            "pLDDT", f"{conf_score}",
            delta=conf_meta.get("source_type", ""),
        )
        m4.metric("EPI PRIORITY", f"{epi}%")

        # ── Coverage Bar ──
        coverage = v.get("coverage_pct", 0)
        cov_float = float(coverage) if isinstance(coverage, (int, float)) else 0
        st.progress(
            min(cov_float / 100, 1.0),
            text=f"Core Coverage: {coverage}% (threshold: 30%)",
        )

        # ── 3D Structure Viewer ──
        pdb_b64 = res.get("pdb_b64", "")
        if pdb_b64:
            st.divider()
            st.subheader("🧬 3D Structure")
            render_3d_viewer(pdb_b64)

        # ── Law Details ──
        st.divider()
        st.subheader("📋 Tier-1 Physical Invariant Audit")

        laws = res.get("tier1", {}).get("laws", [])
        det_laws = [l for l in laws if l.get("method") == "deterministic"]
        heur_laws = [l for l in laws if l.get("method") == "heuristic"]

        col_det, col_heur = st.columns([2, 1])

        with col_det:
            st.markdown(f"**Deterministic Laws** ({len(det_laws)} — VETO authority)")
            for l in det_laws:
                status = l.get("status", "?")
                icon = "✅" if status == "PASS" else "🛑" if status in ("FAIL", "VETO") else "⚠️"
                color = "#00CC66" if status == "PASS" else "#FF4444"
                with st.expander(f"{icon} {l['law_id']}: {l.get('title', '')} — **{status}**"):
                    st.write(f"**Measurement:** {l.get('measurement', 'N/A')}")
                    st.write(f"**Principle:** {l.get('principle', 'N/A')}")
                    st.write(f"**Method:** Deterministic (exact geometric computation)")

        with col_heur:
            st.markdown(f"**Heuristic Laws** ({len(heur_laws)} — Advisory)")
            for l in heur_laws:
                status = l.get("status", "?")
                icon = "✅" if status == "PASS" else "⚠️"
                with st.expander(f"{icon} {l['law_id']}: {l.get('title', '')} — **{status}**"):
                    st.write(f"**Measurement:** {l.get('measurement', 'N/A')}")
                    st.write(f"**Principle:** {l.get('principle', 'N/A')}")
                    st.write(f"**Method:** Heuristic (approximate algorithm)")

        # ── AI Narrative ──
        narratives = res.get("witness_reports", {})
        if narratives and any(narratives.values()):
            st.divider()
            st.subheader("🤖 AI Clinical Assessment")
            st.markdown(f"*Model: {res.get('ai_model_used', 'N/A')}*")

            exec_text = narratives.get("executive", "")
            if exec_text:
                st.info(exec_text)

            col_dd, col_rec = st.columns(2)
            with col_dd:
                with st.expander("🔬 Deep Dive", expanded=False):
                    st.write(narratives.get("deep_dive", "N/A"))
            with col_rec:
                with st.expander("💡 Recommendation", expanded=False):
                    st.write(narratives.get("recommendation", "N/A"))

        # ── Characterization ──
        char = res.get("characterization", {})
        if char:
            st.divider()
            st.subheader("📊 Structural Characterization")
            c1, c2, c3 = st.columns(3)
            c1.metric("α-Helix", f"{char.get('helix', 0)}%")
            c2.metric("β-Sheet", f"{char.get('sheet', 0)}%")
            c3.metric("Coil/Loop", f"{char.get('loop', 0)}%")
            st.caption(
                f"Total atoms: {char.get('total_atoms', 'N/A')} | "
                f"Residues: {char.get('total_residues', 'N/A')} | "
                f"Method: {char.get('method', 'N/A')}"
            )

        # ── PDF Download ──
        if pdb_b64:
            st.divider()
            col_pdf, col_pdb = st.columns(2)

            with col_pdf:
                pdf_b64 = res.get("pdf_b64", "")
                if pdf_b64:
                    try:
                        pdf_bytes = base64.b64decode(pdf_b64)
                        st.download_button(
                            "📄 Download Forensic Dossier (PDF)",
                            data=pdf_bytes,
                            file_name=f"toscanini_{prov.get('source', 'audit')}.pdf",
                            mime="application/pdf",
                            use_container_width=True,
                        )
                    except Exception:
                        st.warning("PDF decode failed")

            with col_pdb:
                try:
                    pdb_bytes = base64.b64decode(res.get("pdb_b64", ""))
                    st.download_button(
                        "🧬 Download Structure (PDB)",
                        data=pdb_bytes,
                        file_name=f"{prov.get('source', 'structure')}.pdb",
                        mime="chemical/x-pdb",
                        use_container_width=True,
                    )
                except Exception:
                    pass

        # ── Reset ──
        st.divider()
        if st.button("🔄 New Search", use_container_width=True, type="secondary"):
            st.session_state.audit_result = None
            st.session_state.candidates = []
            st.rerun()


# ══════════════════════════════════════════════════════════
#  TAB 2: KNOWLEDGE GRAPH
# ══════════════════════════════════════════════════════════
with t_nkg:
    st.subheader("🧠 Knowledge Graph")
    nkg_data = api("GET", "/nkg")
    if nkg_data:
        total_v = nkg_data.get("total_vetoes", 0)
        total_s = nkg_data.get("total_successes", 0)
        total = total_v + total_s

        col_tot, col_pass, col_veto = st.columns(3)
        col_tot.metric("Total Audits", total)
        col_pass.metric("Successes", total_s)
        col_veto.metric("Vetoes", total_v)

        if total == 0:
            st.info("Knowledge graph is empty. Audit structures to populate it.")
        else:
            if nkg_data.get("vetoes"):
                st.markdown("**🛑 Recent Vetoes**")
                for rec in reversed(nkg_data["vetoes"][-10:]):
                    with st.expander(
                        f"{rec.get('timestamp', '?')[:19]} | "
                        f"{rec.get('failure_class', '?')} | "
                        f"Phys: {rec.get('physical_score', 0)}%"
                    ):
                        st.json(rec)

            if nkg_data.get("successes"):
                st.markdown("**✅ Recent Successes**")
                for rec in reversed(nkg_data["successes"][-10:]):
                    with st.expander(
                        f"{rec.get('timestamp', '?')[:19]} | "
                        f"Phys: {rec.get('physical_score', 0)}% | "
                        f"EPI: {rec.get('epi_index', 0)}%"
                    ):
                        st.json(rec)
    else:
        st.warning("Could not connect to Knowledge Graph backend")

with t_laws:
    st.subheader("📜 Toscanini Law Canon")
    laws_data = api("GET", "/laws")
    if laws_data:
        st.caption(
            f"Total: {laws_data.get('total_laws')} laws | "
            f"Deterministic: {laws_data.get('deterministic_count')} (VETO gate) | "
            f"Heuristic: {laws_data.get('heuristic_count')} (Advisory)"
        )
        for lid in sorted(laws_data.get("laws", {}).keys()):
            info = laws_data["laws"][lid]
            sev = info.get("severity", "N/A")
            icon = "🔴" if sev == "FATAL" else "🟡"
            with st.expander(f"{icon} {lid}: {info['title']} [{sev}]"):
                st.write(f"**Principle:** {info.get('principle', 'N/A')}")
                st.write(f"**Threshold:** {info.get('threshold')} {info.get('unit')}")
                st.write(f"**Severity:** {sev}")

    defs = api("GET", "/definitions")
    if defs:
        st.divider()
        st.subheader("📊 Score Definitions")
        for key, d in defs.items():
            with st.expander(d.get("title", key)):
                st.write(d.get("explanation", ""))
