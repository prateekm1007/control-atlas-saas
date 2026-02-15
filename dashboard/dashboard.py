import streamlit as st
import requests, base64, json, os
from style_utils import apply_arena_theme, render_lifecycle_header, score_color

# 🛡️ PIL-HYG-21: Institutional Hygiene
st.set_page_config(page_title="TOSCANINI OS", layout="wide", page_icon="🛡️")
apply_arena_theme()

BACKEND = os.getenv("BACKEND_URL", "http://brain:8000")
API_KEY = os.getenv("TOSCANINI_API_KEY", "")

# 🛡️ PIL-ARC-13: Architecture Lab
WARHEADS = ["NONE", "LINEAR", "MULTIVALENT", "ENGAGEMENT", "MOTIFS", "LINKERS", "CONFORMATIONAL", "METAL"]

def api(method, path, **kwargs):
    headers = kwargs.pop("headers", {})
    if API_KEY: headers["X-API-Key"] = API_KEY
    try:
        r = requests.request(method, f"{BACKEND}{path}", timeout=180, headers=headers, **kwargs)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        st.error(f"Backend Link Error ({path}): {e}")
        return None

def render_dual_engine_visualizer(pdb_b64, verdict):
    """🛡️ PIL-VIS-01: Dual-Engine Visualizer"""
    try:
        pdb_str = base64.decodebytes(pdb_b64.encode()).decode("utf-8", errors="ignore")
        if len(pdb_str.strip()) < 50: return
        import py3Dmol
        import streamlit.components.v1 as components
        
        view = py3Dmol.view(width=800, height=600)
        view.addModel(pdb_str, "pdb")
        
        if verdict == "PASS":
            view.setStyle({"cartoon": {"color": "spectrum"}})
        else:
            view.setStyle({"stick": {"colorscheme": "orangeCarbon"}, "sphere": {"radius": 0.5}})
            
        view.setBackgroundColor("#050505")
        view.zoomTo()
        components.html(view._make_html(), height=620, width=820, scrolling=False)
    except Exception as e:
        st.warning(f"Forensic Visualizer Unavailable: {e}")

# ── Session State (PIL-LCY-19) ─────────────────────────────
if "audit_result" not in st.session_state: st.session_state.audit_result = None
if "candidates" not in st.session_state: st.session_state.candidates = []

# ── 🛡️ PIL-ARC-13: Sidebar Configuration ───────────────────
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/1042/1042307.png", width=80)
    st.title("STATION CONFIG")
    st.divider()
    warhead = st.radio("🎯 ARCHITECTURE WARHEAD", WARHEADS, index=0, help="Selects therapeutic class weighting.")
    st.divider()
    st.caption(f"Engine: v22.5.x")
    st.caption(f"Status: ONLINE")

# ── Header ─────────────────────────────────────────────────
st.markdown('<h1 style="color:#00D4FF;text-align:center;font-family:Courier New;margin-bottom:0px;">🛡️ TOSCANINI</h1>', unsafe_allow_html=True)
current_state = "INSTANTIATED"
if st.session_state.audit_result: current_state = "AUDITED"
elif st.session_state.candidates: current_state = "ACQUIRING"
render_lifecycle_header(current_state)

t_work, t_nkg, t_laws = st.tabs(["⚡ Command Center", "🧠 Knowledge Graph", "📜 Law Canon"])

with t_work:
    if not st.session_state.audit_result:
        # 🛡️ PIL-DSC-04: Discovery Console
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
                            res = api("POST", "/ingest", data={"mode":"Discovery", "candidate_id":c['id'], "t3_category": warhead})
                            if res:
                                st.session_state.audit_result = res
                                st.rerun()
    else:
        # ── 🛡️ RESULTS INTEGRITY DISPLAY ────────────────────
        res = st.session_state.audit_result
        v = res.get("verdict", {})
        binary = v.get("binary", "ERROR")
        
        # Verdict Banner
        if binary == "PASS": st.success("✅ STRUCTURE PASSED: ALL DETERMINISTIC INVARIANTS SATISFIED")
        elif binary == "VETO": st.error(f"🛑 DETERMINISTIC VETO: {', '.join(v.get('killer_laws', []))}")
        else: st.warning(f"⚠️ VERDICT: {binary}")

        # Metric Grid (PIL-EPI-07)
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("DETERMINISTIC", f"{v.get('deterministic_score',0)}%")
        m2.metric("ADVISORY", f"{v.get('advisory_score',0)}%")
        m3.metric("pLDDT (CORE)", f"{round(v.get('confidence_score',0), 1)}")
        m4.metric("EPI PRIORITY", f"{res.get('tier3',{}).get('probability',0)}%")

        # 🛡️ PIL-NAR-11: PhD Witness Report
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
        
        st.caption(f"Narrative Model: {res.get('ai_model_used', 'Internal Fallback')}")

        # 🛡️ PIL-VIS-01: Forensic Visualizer
        st.divider()
        st.subheader("🧬 Forensic 3D Reconstruction")
        if res.get("pdb_b64"): render_dual_engine_visualizer(res["pdb_b64"], binary)

        # 🛡️ PIL-LED-02: Diagnostic Ledger (FIXED SYNC)
        st.divider()
        st.subheader("📋 Tier-1 Diagnostic Ledger")
        laws = res.get("tier1", {}).get("laws", [])
        for l in laws:
            icon = "✅" if l['status'] == "PASS" else "🛑"
            with st.expander(f"{icon} {l['law_id']}: {l['title']} — {l['status']}"):
                # 🛡️ PIL-PAR-14: Numeric Parity (No measurement string)
                st.write(f"**Observed:** {l['observed']} {l['units']}")
                st.write(f"**Threshold:** {l['operator']} {l['threshold']} {l['units']}")
                st.write(f"**Deviation:** {l['deviation']}")
                st.write(f"**Evaluation Scope:** {l['sample_size']} {l['scope']}")
                st.write(f"**Scientific Principle:** {l['principle']}")

        # 🛡️ PIL-PAR-14: Forensic Dossier (PDF) Download
        st.divider()
        if res.get("pdf_b64"):
            pdf_bytes = base64.b64decode(res["pdf_b64"])
            st.download_button(
                label="📄 DOWNLOAD FORENSIC DOSSIER (PDF)",
                data=pdf_bytes,
                file_name=f"TOSCANINI_DOSSIER_{res.get('governance', {}).get('audit_id', 'ARTIFACT')}.pdf",
                mime="application/pdf",
                use_container_width=True
            )
        else:
            st.error("PIL-PAR-14 VIOLATION: PDF byte-stream not found in payload.")
        
        if st.button("🔄 RETURN TO COMMAND CENTER", use_container_width=True):
            st.session_state.audit_result = None
            st.rerun()

with t_nkg:
    st.subheader("🧠 Knowledge Graph (NKG)")
    nkg_data = api("GET", "/nkg")
    if nkg_data:
        st.metric("Total System Vetoes", nkg_data.get('total_vetoes', 0))
        if nkg_data.get("vetoes"): st.json(nkg_data["vetoes"])
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
