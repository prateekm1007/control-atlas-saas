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


# ── Session State ──────────────────────────────────────────
if "audit_result" not in st.session_state:
    st.session_state.audit_result = None
if "candidates" not in st.session_state:
    st.session_state.candidates = []
if "selected_candidate" not in st.session_state:
    st.session_state.selected_candidate = None

# ── Header ─────────────────────────────────────────────────
st.markdown(
    '<h1 style="color:#00D4FF;text-align:center;font-family:Courier New;margin-bottom:0px;">'
    '🛡️ TOSCANINI</h1>',
    unsafe_allow_html=True,
)
current_state = "INSTANTIATED"
if st.session_state.audit_result:
    current_state = "AUDITED"
elif st.session_state.candidates:
    current_state = "ACQUIRING"
render_lifecycle_header(current_state)

# ── Tabs ───────────────────────────────────────────────────
t_work, t_nkg, t_laws = st.tabs(["⚡ Command Center", "🧠 Knowledge Graph", "📜 Law Canon"])

# ══════════════════════════════════════════════════════════
#  TAB 1: COMMAND CENTER
# ══════════════════════════════════════════════════════════
with t_work:
    if st.session_state.audit_result is None:
        # ── SEARCH PHASE ──
        st.subheader("🔍 Structure Discovery")

        col_search, col_upload = st.columns([2, 1])

        with col_search:
            q = st.text_input("Search UniProt / AlphaFold", placeholder="e.g. insulin, EGFR, P53...")
            if st.button("🔎 Search", type="primary", use_container_width=True):
                if q.strip():
                    with st.spinner("Searching AlphaFold DB..."):
                        resp = api("POST", "/search", data={"query": q})
                        if resp and isinstance(resp, dict):
                            st.session_state.candidates = resp.get("results", [])
                        elif resp and isinstance(resp, list):
                            st.session_state.candidates = resp
                        else:
                            st.session_state.candidates = []
                    if st.session_state.candidates:
                        st.success(f"Found {len(st.session_state.candidates)} candidates")
                    else:
                        st.warning("No results found")

        with col_upload:
            st.markdown("**Or upload PDB/CIF**")
            uploaded = st.file_uploader("Upload structure", type=["pdb", "cif"], label_visibility="collapsed")

        # ── WARHEAD SELECTOR ──
        st.divider()
        warhead = st.selectbox(
            "🎯 Architecture Warhead (T3 Category)",
            WARHEADS,
            help="Selects the architectural context for EPI prioritization scoring."
        )

        # ── CANDIDATE LIST ──
        if st.session_state.candidates:
            st.subheader(f"📋 Candidates ({len(st.session_state.candidates)})")
            for i, c in enumerate(st.session_state.candidates):
                col_label, col_btn = st.columns([4, 1])
                with col_label:
                    st.markdown(f"**{c.get('id', 'N/A')}** — {c.get('label', 'Unknown')}")
                with col_btn:
                    if st.button("⚡ Audit", key=f"audit_{i}", use_container_width=True):
                        with st.spinner(f"Auditing {c['id']}..."):
                            res = api("POST", "/ingest", data={
                                "mode": "Discovery",
                                "candidate_id": c["id"],
                                "t3_category": warhead,
                            })
                            if res:
                                st.session_state.audit_result = res
                                st.rerun()

        # ── DIRECT UPLOAD AUDIT ──
        if uploaded is not None:
            if st.button("⚡ Audit Uploaded File", type="primary", use_container_width=True):
                with st.spinner("Auditing uploaded structure..."):
                    res = api("POST", "/ingest",
                              data={"mode": "Upload", "candidate_id": uploaded.name, "t3_category": warhead},
                              files={"file": (uploaded.name, uploaded.getvalue())})
                    if res:
                        st.session_state.audit_result = res
                        st.rerun()

    else:
        # ══════════════════════════════════════════════════════
        #  RESULTS DISPLAY
        # ══════════════════════════════════════════════════════
        res = st.session_state.audit_result
        v = res.get("verdict", {})
        binary = v.get("binary", "ERROR")
        prov = res.get("provenance", {})

        # ── Verdict Banner ──
        if binary == "PASS":
            st.success("✅ STRUCTURE PASSED — All deterministic invariants satisfied")
        elif binary == "VETO":
            killers = v.get("killer_laws", [])
            msg = ", ".join(killers) if killers else "Deterministic violation"
            st.error(f"🛑 DETERMINISTIC VETO — {msg}")
        elif binary == "INDETERMINATE":
            st.warning("⚠️ INDETERMINATE — Insufficient coverage or fringe escalation")
        else:
            st.error(f"❌ ERROR: {binary}")

        # ── Source info ──
        st.caption(f"Source: {prov.get('source', 'N/A')} | Bytes: {prov.get('byte_count', 0):,} | Version: {prov.get('station_version', '?')}")

        # ── 4-Metric Grid ──
        m1, m2, m3, m4 = st.columns(4)
        det_score = v.get("deterministic_score", 0)
        adv_score = v.get("advisory_score", 0)
        conf_score = v.get("confidence_score", 0)
        epi = res.get("tier3", {}).get("probability", 0)

        m1.metric("DETERMINISTIC", f"{det_score}%",
                  delta=f"{v.get('det_passed', 0)}/{v.get('det_total', 0)} laws")
        m2.metric("ADVISORY", f"{adv_score}%",
                  delta=f"{v.get('heur_passed', 0)}/{v.get('heur_total', 0)} laws")
        m3.metric("pLDDT", f"{conf_score}",
                  delta=res.get("confidence_meta", {}).get("source_type", ""))
        m4.metric("EPI PRIORITY", f"{epi}%")

        # ── Coverage ──
        coverage = v.get("coverage_pct", "N/A")
        st.progress(min(float(coverage) / 100 if isinstance(coverage, (int, float)) else 0, 1.0),
                    text=f"Core Coverage: {coverage}%")

        # ── Law Details ──
        st.divider()
        st.subheader("📋 Tier-1 Law Results")

        laws = res.get("tier1", {}).get("laws", [])
        det_laws = [l for l in laws if l.get("method") == "deterministic"]
        heur_laws = [l for l in laws if l.get("method") == "heuristic"]

        col_det, col_heur = st.columns([2, 1])

        with col_det:
            st.markdown("**Deterministic Laws** (VETO authority)")
            for l in det_laws:
                status = l.get("status", "?")
                icon = "✅" if status == "PASS" else "🛑" if status in ("FAIL", "VETO") else "⚠️"
                with st.expander(f"{icon} {l['law_id']}: {l.get('title', '')} — {status}"):
                    st.write(f"**Measurement:** {l.get('measurement', 'N/A')}")
                    st.write(f"**Principle:** {l.get('principle', 'N/A')}")

        with col_heur:
            st.markdown("**Heuristic Laws** (Advisory only)")
            for l in heur_laws:
                status = l.get("status", "?")
                icon = "✅" if status == "PASS" else "⚠️"
                with st.expander(f"{icon} {l['law_id']}: {l.get('title', '')} — {status}"):
                    st.write(f"**Measurement:** {l.get('measurement', 'N/A')}")

        # ── AI Narrative ──
        narratives = res.get("witness_reports", {})
        if narratives:
            st.divider()
            st.subheader("🤖 AI Assessment")
            st.info(narratives.get("executive", ""))
            with st.expander("Deep Dive"):
                st.write(narratives.get("deep_dive", ""))
            with st.expander("Recommendation"):
                st.write(narratives.get("recommendation", ""))
            st.caption(f"Model: {res.get('ai_model_used', 'N/A')}")

        # ── PDF Download ──
        pdf_b64 = res.get("pdf_b64", "")
        if pdf_b64:
            st.divider()
            try:
                pdf_bytes = base64.b64decode(pdf_b64)
                st.download_button(
                    "📄 Download Forensic Dossier (PDF)",
                    data=pdf_bytes,
                    file_name=f"toscanini_audit_{prov.get('source', 'structure')}.pdf",
                    mime="application/pdf",
                    use_container_width=True,
                )
            except Exception:
                st.warning("PDF generation failed for this audit")

        # ── Reset ──
        st.divider()
        if st.button("🔄 New Search", use_container_width=True, type="secondary"):
            st.session_state.audit_result = None
            st.session_state.candidates = []
            st.rerun()


# ══════════════════════════════════════════════════════════
#  TAB 2: KNOWLEDGE GRAPH (placeholder)
# ══════════════════════════════════════════════════════════
with t_nkg:
    st.subheader("🧠 Negative Knowledge Graph")
    st.info("NKG records accumulate as structures are audited. Historical failure patterns will appear here.")
    st.caption("Database: /app/nkg/piu_moat.jsonl")


# ══════════════════════════════════════════════════════════
#  TAB 3: LAW CANON
# ══════════════════════════════════════════════════════════
with t_laws:
    st.subheader("📜 Law Canon")
    laws_data = api("GET", "/laws")
    if laws_data:
        st.caption(f"Total: {laws_data.get('total_laws')} laws | "
                   f"Deterministic: {laws_data.get('deterministic_count')} | "
                   f"Heuristic: {laws_data.get('heuristic_count')}")
        for lid in sorted(laws_data.get("laws", {}).keys()):
            info = laws_data["laws"][lid]
            with st.expander(f"{lid}: {info['title']}"):
                st.write(f"**Principle:** {info.get('principle', 'N/A')}")
                st.write(f"**Threshold:** {info.get('threshold')} {info.get('unit')}")
                st.write(f"**Severity:** {info.get('severity', 'N/A')}")

    defs = api("GET", "/definitions")
    if defs:
        st.divider()
        st.subheader("📊 Score Definitions")
        for key, d in defs.items():
            with st.expander(d.get("title", key)):
                st.write(d.get("explanation", ""))
