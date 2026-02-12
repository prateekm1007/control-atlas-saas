import streamlit as st
import requests, base64, json, os, html
from style_utils import apply_arena_theme, render_lifecycle_header
import streamlit.components.v1 as components

st.set_page_config(page_title="TOSCANINI OS", layout="wide", page_icon="🛡️")
apply_arena_theme()
backend_url = os.getenv("BACKEND_URL", "http://brain:8000")

def safe_req(method, url, **kwargs):
    try:
        r = requests.request(method, url, timeout=120, **kwargs)
        r.raise_for_status(); return r.json()
    except Exception as e: st.error(f"Fault: {e}"); return None

if "audit_result" not in st.session_state: st.session_state.audit_result = None
if "candidates" not in st.session_state: st.session_state.candidates = []

with st.sidebar:
    st.header("🔬 DESIGN LAB")
    WARHEADS = ["NONE", "LINEAR", "MULTIVALENT", "ENGAGEMENT", "MOTIFS", "LINKERS", "CONFORMATIONAL", "METAL"]
    t3_choice = st.selectbox("ARCHITECTURE INTENT", WARHEADS, index=0)
    
    if st.button("🚀 APPLY ARCHITECTURE", use_container_width=True, type="primary"):
        if st.session_state.get("last_id"):
            res = safe_req("POST", f"{backend_url}/ingest", data={"mode":"Discovery","candidate_id":st.session_state.last_id,"t3_category":t3_choice})
            if res: st.session_state.audit_result = res; st.rerun()

st.markdown('<h1 style="color:#00D4FF;text-align:center;font-family:Courier New;margin-bottom:0px;">🛡️ TOSCANINI</h1>', unsafe_allow_html=True)
current_res = st.session_state.audit_result or {}
render_lifecycle_header(current_res.get('governance', {}).get('state', 'INSTANTIATED'))

t_work, t_nkg, t_pillars = st.tabs(["⚡ Command Center", "🧠 Knowledge Graph", "🛡️ 21 Pillars"])

with t_work:
    if not st.session_state.audit_result:
        # DUAL MODE: Search or Upload
        m_search, m_upload = st.tabs(["🎯 Discovery Search", "📤 Institutional Upload"])
        
        with m_search:
            q = st.text_input("QUERY", placeholder="e.g. insulin, kinase...")
            if st.button("Fetch Exhaustive Results", type="primary"):
                resp = safe_req("POST", f"{backend_url}/search", data={"query": q})
                st.session_state.candidates = resp if isinstance(resp, list) else []
            for i, c in enumerate(st.session_state.candidates):
                if st.button(f"📥 {c['label']}", key=f"c_{i}"):
                    st.session_state.last_id = c['id']
                    res = safe_req("POST", f"{backend_url}/ingest", data={"mode":"Discovery","candidate_id":c['id'],"t3_category":t3_choice})
                    if res: st.session_state.audit_result = res; st.rerun()
        
        with m_upload:
            st.info("Manual PDB/mmCIF Ingestion")
            up_file = st.file_uploader("Select Structure File", type=["pdb", "cif", "mmcif"])
            if up_file and st.button("Execute Manual Audit"):
                files = {"file": (up_file.name, up_file.getvalue())}
                res = safe_req("POST", f"{backend_url}/ingest", data={"mode":"Upload","candidate_id":up_file.name,"t3_category":t3_choice}, files=files)
                if res: st.session_state.audit_result = res; st.rerun()

    else:
        # Results display...
        res = st.session_state.audit_result; v = res.get("verdict", {}); defs = res.get("definitions", {})
        m1, m2, m3 = st.columns(3)
        m1.metric("PHYSICAL INTEGRITY", f"{v.get('physical_score')}%", help=defs.get('PHYSICAL_SCORE',{}).get('explanation'))
        m2.metric("PRIORITY INDEX (EPI)", f"{res.get('tier3',{}).get('probability')}%", help=defs.get('STRATEGIC_SCORE',{}).get('explanation'))
        m3.metric("ML CONFIDENCE", f"{v.get('confidence_score')}%", help=defs.get('ML_CONFIDENCE',{}).get('explanation'))
        
        col_l, col_r = st.columns([1, 1.5])
        with col_l:
            st.subheader("📋 Audit Summary")
            reports = res.get("witness_reports", {})
            st.info(f"**EXECUTIVE SUMMARY**\n\n{reports.get('executive', 'Analyzing...')}")
            with st.expander("🔍 STRUCTURAL ASSESSMENT", expanded=True):
                st.write(reports.get("deep_dive", "Gathering data..."))
            st.success(f"**RESOURCE RECOMMENDATION:** {reports.get('recommendation', 'Evaluating...')}")
            if res.get("pdf_b64"): 
                st.download_button("📥 DOWNLOAD 6-PAGE DECISION DOSSIER", base64.b64decode(res["pdf_b64"]), file_name="toscanini.pdf", use_container_width=True, type="primary")
        with col_r:
            pdb_data = base64.b64decode(res.get("pdb_b64", "")).decode()
            pdb_escaped = html.escape(json.dumps(pdb_data))
            style = "{stick:{colorscheme:'orangeCarbon'}}" if v.get('binary')=="VETO" else "{cartoon:{color:'spectrum'}}"
            html_block = f"<div id='v' style='height:400px;background:#000;' data-pdb='{pdb_escaped}'></div><script src='https://3Dmol.org/build/3Dmol-min.js'></script><script>var el=document.getElementById('v'); var pdb=JSON.parse(el.dataset.pdb); var v=$3Dmol.createViewer(el); v.addModel(pdb,'pdb'); v.setStyle({{}},{style}); v.zoomTo(); v.render();</script>"
            components.html(html_block, height=420)
        if st.button("🔄 New Search", use_container_width=True):
            st.session_state.audit_result = None; st.session_state.candidates = []; st.rerun()

with t_nkg:
    nkg = safe_req("GET", f"{backend_url}/nkg/records")
    if nkg:
        c1, c2 = st.columns(2)
        with c1:
            st.error("VETO RECORDS")
            for r in nkg.get('vetoes', []): st.code(json.dumps(r, indent=2), language="json")
        with c2:
            st.success("SUCCESS RECORDS")
            for r in nkg.get('successes', []): st.code(json.dumps(r, indent=2), language="json")

with t_pillars:
    defs = safe_req("GET", f"{backend_url}/definitions")
    if defs:
        for title, desc in defs.get('PILLAR_MANIFEST', []):
            with st.expander(f"**{title}**"): st.write(desc)
