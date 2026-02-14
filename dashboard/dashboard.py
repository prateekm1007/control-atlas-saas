import streamlit as st
import requests, base64, json, os, html
from style_utils import apply_arena_theme, render_lifecycle_header, score_color
import streamlit.components.v1 as components

st.set_page_config(page_title="TOSCANINI OS", layout="wide", page_icon="🛡️")
apply_arena_theme()

BACKEND = os.getenv("BACKEND_URL", "http://brain:8000")
API_KEY = os.getenv("TOSCANINI_API_KEY", "")
WARHEADS = ["NONE", "LINEAR", "MULTIVALENT", "ENGAGEMENT", "MOTIFS", "LINKERS", "CONFORMATIONAL", "METAL"]

def api(method, path, **kwargs):
    headers = kwargs.pop("headers", {})
    if API_KEY: headers["X-API-Key"] = API_KEY
    try:
        r = requests.request(method, f"{BACKEND}{path}", timeout=180, headers=headers, **kwargs)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        st.error(f"Backend Link Error: {e}")
        return None

if "audit_result" not in st.session_state: st.session_state.audit_result = None
if "candidates" not in st.session_state: st.session_state.candidates = []

st.markdown('<h1 style="color:#00D4FF;text-align:center;font-family:Courier New;margin-bottom:0px;">🛡️ TOSCANINI</h1>', unsafe_allow_html=True)
render_lifecycle_header(st.session_state.audit_result.get('governance', {}).get('state', 'INSTANTIATED') if st.session_state.audit_result else "INSTANTIATED")

t_work, t_nkg, t_laws = st.tabs(["⚡ Command Center", "🧠 Knowledge Graph", "📜 Law Canon"])

with t_work:
    if not st.session_state.audit_result:
        q = st.text_input("QUERY", placeholder="e.g. insulin, kinase...")
        if st.button("Fetch Exhaustive Results", type="primary"):
            resp = api("POST", "/search", data={"query": q})
            st.session_state.candidates = resp if isinstance(resp, list) else []
        
        for i, c in enumerate(st.session_state.candidates):
            if st.button(f"📥 {c['label']}", key=f"c_{i}"):
                with st.spinner("Executing Forensic Audit..."):
                    res = api("POST", "/ingest", data={"mode":"Discovery","candidate_id":c['id']})
                    if res: 
                        st.session_state.audit_result = res
                        st.rerun()
    else:
        res = st.session_state.audit_result
        v = res.get("verdict", {})
        binary = v.get("binary", "ERROR")
        
        if binary == "VETO":
            killer_list = v.get('killer_laws', [])
            msg = ", ".join(killer_list) if killer_list else "Deterministic Violation"
            # --- FIXED LINE BELOW ---
            st.error(f"🛑 DETERMINISTIC VETO: {msg}")
        else:
            st.success("✅ ALL DETERMINISTIC INVARIANTS SATISFIED")

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("DETERMINISTIC", f"{v.get('deterministic_score',0)}%")
        m2.metric("ADVISORY", f"{v.get('advisory_score',0)}%")
        m3.metric("pLDDT", f"{v.get('confidence_score',0)}")
        m4.metric("EPI PRIORITY", f"{res.get('tier3',{}).get('probability',0)}%")

        if st.button("🔄 New Search", use_container_width=True):
            st.session_state.audit_result = None
            st.rerun()

with t_laws:
    laws_data = api("GET", "/laws")
    if laws_data:
        for lid, info in laws_data.get("laws", {}).items():
            with st.expander(f"{lid}: {info['title']}"):
                st.write(f"**Principle:** {info.get('principle', 'N/A')}")
                st.write(f"**Threshold:** {info.get('threshold')} {info.get('unit')}")
