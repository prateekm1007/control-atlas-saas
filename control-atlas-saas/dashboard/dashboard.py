import base64
import os
from typing import Dict, Optional

import py3Dmol
import requests
import streamlit as st
from stmol import showmol


st.set_page_config(page_title="Toscanini Pipeline Console", layout="wide")

BACKEND_URL = os.getenv("BACKEND_URL", "http://localhost:8000")

if "audit_result" not in st.session_state:
    st.session_state.audit_result = None
if "program_cache" not in st.session_state:
    st.session_state.program_cache = {}
if "selected_anchor" not in st.session_state:
    st.session_state.selected_anchor = None


def fetch_stats() -> int:
    try:
        resp = requests.get(f"{BACKEND_URL}/stats", timeout=10)
        resp.raise_for_status()
        return int(resp.json().get("unique_pius", 0))
    except Exception:
        return 0


def ingest_structure(payload: Dict, files: Optional[Dict] = None) -> Optional[Dict]:
    try:
        response = requests.post(f"{BACKEND_URL}/ingest", data=payload, files=files, timeout=180)
        response.raise_for_status()
        return response.json()
    except Exception as exc:
        st.error(f"⚠️ Forensic ingest failed: {exc}")
        return None


def fetch_dossier(program_id: str) -> Optional[bytes]:
    try:
        response = requests.get(f"{BACKEND_URL}/programs/{program_id}/dossier", timeout=30)
        response.raise_for_status()
        return response.content
    except Exception as exc:
        st.error(f"⚠️ Dossier retrieval failed: {exc}")
        return None


def program_state_panel(program_data: Dict) -> None:
    st.subheader("Program Lifecycle (S0 → S7)")
    state = program_data.get("state")
    hash_chain = program_data.get("hash_chain") or []
    steps = [
        "S0_INGESTED",
        "S1_ACQUISITION_COMPLETE",
        "S2_FOLDING_COMPLETE",
        "S3_TIER1_COMPLETE",
        "S4_ARCHITECTURE_DERIVED",
        "S5_TIER3_COMPLETE",
        "S6_FORENSIC_AGGREGATION",
        "S7_DOSSIER_READY",
    ]
    cols = st.columns(4)
    for idx, step in enumerate(steps):
        col = cols[idx % 4]
        marker = "✅" if state and steps.index(step) <= steps.index(state) else "⬜"
        col.write(f"{marker} {step}")
    integrity = "OK" if hash_chain else "UNVERIFIED"
    st.caption(f"Hash-chain entries: {len(hash_chain)} | Integrity: {integrity}")


def render_viewer(pdb_bytes: bytes, fmt: str, verdict: str, anchor: Optional[Dict]) -> None:
    pdb_str = pdb_bytes.decode("utf-8", errors="ignore")
    view = py3Dmol.view(width=700, height=520)
    view.addModel(pdb_str, fmt)
    if verdict == "VETO":
        view.setStyle({}, {"stick": {"radius": 0.2}, "sphere": {"radius": 0.6, "color": "red"}})
    else:
        view.setStyle({}, {"cartoon": {"color": "spectrum"}})
    if anchor:
        pos = anchor.get("pos")
        if pos:
            view.addSphere({"center": {"x": pos[0], "y": pos[1], "z": pos[2]}, "radius": 3, "color": "red", "wireframe": True})
            view.zoomTo({"center": {"x": pos[0], "y": pos[1], "z": pos[2]}})
    view.zoomTo()
    view.setBackgroundColor("0x0b0f17")
    showmol(view, height=520, width=700)


def ledger_section(laws: list) -> None:
    st.subheader("Diagnostic Ledger")
    for law in laws:
        status = law.get("status", "PASS")
        icon = "✅" if status == "PASS" else "❌"
        with st.expander(f"{icon} {law.get('law_id')} — {law.get('title')}"):
            st.write(law.get("measurement"))
            if law.get("anchor"):
                if st.button("🔍 Forensic Zoom", key=f"anchor_{law.get('law_id')}"):
                    st.session_state.selected_anchor = law.get("anchor")
                    st.experimental_rerun()
            st.caption(f"Principle: {law.get('principle')}")


st.title("Toscanini — Forensic Falsification-as-a-Service")
st.caption("The filter is more valuable than the generator. Schrödinger-Class credibility, every time.")

with st.sidebar:
    st.header("📉 NKG Intelligence")
    nkg_count = fetch_stats()
    st.metric("Failure Moat Entries", nkg_count)
    st.caption("Live count from piu_moat.jsonl")

    st.divider()
    st.header("Sovereign Ingestion")
    existing_program = st.text_input("Load Existing Program ID")
    if existing_program:
        cached = st.session_state.program_cache.get(existing_program)
        if cached:
            st.session_state.audit_result = cached
            st.success(f"Loaded program {existing_program} from session cache.")
        else:
            st.warning("Program not in session cache. Ingest again or keep session active.")

    afdb_id = st.text_input("AlphaFold DB ID", placeholder="AF-P12345-F1")
    uploaded_file = st.file_uploader("Upload PDB/CIF", type=["pdb", "cif"])
    run_audit = st.button("Ingest & Run Full Audit")

    if run_audit:
        payload = {"mode": "Discovery", "candidate_id": afdb_id or "manual"}
        files = None
        if uploaded_file:
            payload["mode"] = "Upload"
            files = {"file": (uploaded_file.name, uploaded_file.getvalue())}
        result = ingest_structure(payload, files=files)
        if result:
            st.session_state.audit_result = result
            program_id = result.get("program", {}).get("program_id")
            if program_id:
                st.session_state.program_cache[program_id] = result


if st.session_state.audit_result:
    result = st.session_state.audit_result
    if result.get("verdict") == "ERROR":
        st.error(result.get("details"))
    else:
        program = result.get("program", {})
        verdict = result.get("verdict")
        prob = program.get("tier3", {}).get("probability", 0.0)

        st.subheader("Program Continuity")
        st.write(f"Program ID: `{program.get('program_id')}`")

        program_state_panel(program)

        banner_color = "🟩" if verdict == "PASS" else "🟥"
        st.markdown(
            f"### {banner_color} FINAL VERDICT: **{verdict}** — Schrödinger-Class Badge",
        )
        st.metric("Tier-3 Probability", f"{prob * 100:.2f}%")

        col_left, col_right = st.columns([1, 1])
        with col_left:
            st.subheader("Architecture Category")
            arch = program.get("architecture") or {}
            st.write(f"**{arch.get('category', 'UNASSIGNED')}**")
            st.caption(arch.get("rationale", "No rationale recorded."))

            st.subheader("Tier-3 Probability Breakdown")
            breakdown = program.get("tier3") or {}
            st.progress(min(max(prob, 0.0), 1.0))
            st.json(
                {
                    "P_base": breakdown.get("p_base"),
                    "P_phys": breakdown.get("p_phys"),
                    "P_ml": breakdown.get("p_ml"),
                    "Tax": breakdown.get("tax"),
                    "Weight_Derivation": breakdown.get("weight_derivation"),
                    "Final_P": breakdown.get("probability"),
                }
            )

        with col_right:
            st.subheader("Dual-Engine Visualizer")
            if result.get("pdb_b64"):
                pdb_bytes = base64.b64decode(result.get("pdb_b64"))
                render_viewer(pdb_bytes, result.get("ext", "pdb"), verdict, st.session_state.selected_anchor)
            else:
                st.warning("No structure data available for visualization.")

        ledger_section(result.get("laws", []))

        st.subheader("Dossier Delivery")
        dossier_bytes = None
        if program.get("program_id"):
            if st.button("Download Nature-Grade Dossier"):
                dossier_bytes = fetch_dossier(program["program_id"])
            if dossier_bytes:
                st.download_button(
                    "Save PDF",
                    dossier_bytes,
                    file_name=f"{program['program_id']}_dossier.pdf",
                    mime="application/pdf",
                )
                pdf_b64 = base64.b64encode(dossier_bytes).decode("utf-8")
                st.components.v1.html(
                    f"<iframe src='data:application/pdf;base64,{pdf_b64}' width='100%' height='500'></iframe>",
                    height=520,
                )
