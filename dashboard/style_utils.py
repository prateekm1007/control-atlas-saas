import streamlit as st

def apply_arena_theme():
    st.markdown("""
    <style>
        .stApp { background-color: #000000; color: #E0E0E0; }
        .stMetric { background-color: #111; border-left: 3px solid #00D4FF; padding: 10px; border-radius: 5px; }
        .stMetric label { color: #888; }
        .stMetric [data-testid="stMetricValue"] { color: #00D4FF; font-size: 1.8rem; }
        .stExpander { border: 1px solid #333; background-color: #050505; }
        .stTabs [data-baseweb="tab"] { color: #888; }
        .stTabs [aria-selected="true"] { color: #00D4FF !important; border-bottom-color: #00D4FF !important; }
    </style>
    """, unsafe_allow_html=True)

def render_lifecycle_header(state: str):
    states = ["INSTANTIATED", "ACQUIRING", "AUDITED", "DECIDED", "SEALED"]
    chips = ""
    for s in states:
        color = "#00D4FF" if s == state else "#222"
        txt_color = "#000" if s == state else "#555"
        weight = "bold" if s == state else "normal"
        chips += f'<span style="background:{color};color:{txt_color};padding:4px 12px;margin:2px;border-radius:12px;font-size:0.75rem;font-weight:{weight};">{s}</span>'
    st.markdown(f'<div style="text-align:center;margin-bottom:20px;">{chips}</div>', unsafe_allow_html=True)
