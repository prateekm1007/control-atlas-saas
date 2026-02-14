import streamlit as st


def score_color(value: float, max_val: float = 100) -> str:
    """Return hex color based on score value. Green=good, Red=bad."""
    ratio = value / max_val if max_val > 0 else 0
    if ratio >= 0.95:
        return "#00CC66"   # Green — passing
    elif ratio >= 0.70:
        return "#FFB833"   # Amber — concerning
    else:
        return "#FF4444"   # Red — failing


def apply_arena_theme():
    st.markdown("""
    <style>
    .stApp { background-color: #000000; color: #E0E0E0; }
    .stExpander { border: 1px solid #333; background-color: #050505; }
    .stTabs [data-baseweb="tab"] { color: #888; }
    .stTabs [aria-selected="true"] { color: #00D4FF !important; border-bottom-color: #00D4FF !important; }
    div[data-testid="stMetricValue"] { font-size: 1.5rem; }
    .stDivider { border-color: #222; }
    </style>
    """, unsafe_allow_html=True)


def render_lifecycle_header(state: str):
    states = ["INSTANTIATED", "ACQUIRING", "AUDITED", "DECIDED", "SEALED"]
    chips = ""
    for i, s in enumerate(states):
        s_idx = states.index(state) if state in states else 0
        if i < s_idx:
            color, txt = "#006644", "#fff"        # Completed — dark green
        elif i == s_idx:
            color, txt = "#00D4FF", "#000"        # Current — cyan
        else:
            color, txt = "#1a1a1a", "#444"        # Future — dim
        weight = "bold" if i == s_idx else "normal"
        chips += (f'<span style="background:{color};color:{txt};padding:4px 12px;margin:2px;'
                  f'border-radius:12px;font-size:0.7rem;font-weight:{weight};'
                  f'display:inline-block;">{s}</span>')
    st.markdown(f'<div style="text-align:center;margin-bottom:15px;">{chips}</div>', unsafe_allow_html=True)
