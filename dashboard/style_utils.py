import streamlit as st

# ═══════════════════════════════════════════════════════════════════════════════
# PIL-SOL-01: Solar Flare Theme Registry
# ═══════════════════════════════════════════════════════════════════════════════

THEMES = {
    "Biotech Noir": {
        "bg": "#000000",
        "text": "#E0E0E0",
        "accent": "#00D4FF",
        "accent_dark": "#008FB3",
        "surface": "#050505",
        "border": "#333333",
        "muted": "#888888",
        "chip_done_bg": "#006644",
        "chip_done_text": "#ffffff",
        "chip_active_text": "#000000",
        "chip_future_bg": "#1a1a1a",
        "chip_future_text": "#444444",
        "viz_bg": "#050505",
        "metric_label": "#aaaaaa",
        "divider": "#222222",
        "sidebar_bg": "#0a0a0a",
        "input_bg": "#111111",
        "input_border": "#333333",
        "toggle_track": "#333333",
        "toggle_thumb": "#00D4FF",
    },
    "Clinical White": {
        "bg": "#FFFFFF",
        "text": "#1e293b",
        "accent": "#059669",
        "accent_dark": "#047857",
        "surface": "#f8fafc",
        "border": "#e2e8f0",
        "muted": "#64748b",
        "chip_done_bg": "#d1fae5",
        "chip_done_text": "#065f46",
        "chip_active_text": "#ffffff",
        "chip_future_bg": "#f1f5f9",
        "chip_future_text": "#94a3b8",
        "viz_bg": "#f0f4f8",
        "metric_label": "#475569",
        "divider": "#e2e8f0",
        "sidebar_bg": "#f1f5f9",
        "input_bg": "#ffffff",
        "input_border": "#cbd5e1",
        "toggle_track": "#cbd5e1",
        "toggle_thumb": "#059669",
    },
}


def get_active_theme() -> dict:
    if "noir_toggle" in st.session_state:
        name = "Biotech Noir" if st.session_state.noir_toggle else "Clinical White"
    else:
        name = "Biotech Noir"
    return THEMES.get(name, THEMES["Biotech Noir"])


def get_active_theme_name() -> str:
    if "noir_toggle" in st.session_state:
        return "Biotech Noir" if st.session_state.noir_toggle else "Clinical White"
    return "Biotech Noir"


def score_color(value: float, max_val: float = 100) -> str:
    ratio = value / max_val if max_val > 0 else 0
    if ratio >= 0.95:
        return "#00CC66"
    elif ratio >= 0.70:
        return "#FFB833"
    else:
        return "#FF4444"


def apply_arena_theme():
    t = get_active_theme()
    st.markdown(f"""
    <style>
    /* ── Base ───────────────────────────────────── */
    .stApp {{ background-color: {t["bg"]}; color: {t["text"]}; }}

    /* ── Toggle Force-Visibility ────────────────── */
    div[data-testid="stToggle"] label > div[role="checkbox"] {{
        background-color: {t["toggle_track"]} !important;
        border: 2px solid {t["border"]} !important;
    }}
    div[data-testid="stToggle"] label > div[role="checkbox"][aria-checked="true"] {{
        background-color: {t["accent"]} !important;
        border-color: {t["accent"]} !important;
    }}
    div[data-testid="stToggle"] label > div[role="checkbox"] > div {{
        background-color: {t["text"]} !important;
    }}

    /* ── Sidebar — FULL OVERRIDE ────────────────── */
    section[data-testid="stSidebar"] {{
        background-color: {t["sidebar_bg"]} !important;
    }}
    section[data-testid="stSidebar"] p,
    section[data-testid="stSidebar"] span,
    section[data-testid="stSidebar"] label,
    section[data-testid="stSidebar"] div,
    section[data-testid="stSidebar"] h1,
    section[data-testid="stSidebar"] h2,
    section[data-testid="stSidebar"] h3,
    section[data-testid="stSidebar"] li,
    section[data-testid="stSidebar"] .stMarkdown,
    section[data-testid="stSidebar"] .stMarkdown p,
    section[data-testid="stSidebar"] .stRadio label,
    section[data-testid="stSidebar"] .stRadio label span,
    section[data-testid="stSidebar"] .stRadio label p,
    section[data-testid="stSidebar"] .stRadio div[role="radiogroup"] label,
    section[data-testid="stSidebar"] .stRadio div[role="radiogroup"] label div,
    section[data-testid="stSidebar"] .stRadio div[role="radiogroup"] label div p,
    section[data-testid="stSidebar"] [data-testid="stWidgetLabel"],
    section[data-testid="stSidebar"] [data-testid="stWidgetLabel"] p,
    section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"],
    section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p,
    section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] span,
    section[data-testid="stSidebar"] .stCaption,
    section[data-testid="stSidebar"] .stCaption p {{
        color: {t["text"]} !important;
    }}
    section[data-testid="stSidebar"] hr {{
        border-color: {t["divider"]} !important;
    }}

    /* ── Expanders ──────────────────────────────── */
    .stExpander {{ border: 1px solid {t["border"]}; background-color: {t["surface"]}; }}
    details summary span {{ color: {t["text"]} !important; }}
    .stExpander .stMarkdown, .stExpander .stMarkdown p {{ color: {t["text"]}; }}

    /* ── Tabs ───────────────────────────────────── */
    .stTabs [data-baseweb="tab"] {{ color: {t["muted"]}; }}
    .stTabs [aria-selected="true"] {{
        color: {t["accent"]} !important;
        border-bottom-color: {t["accent"]} !important;
    }}

    /* ── Metrics ────────────────────────────────── */
    div[data-testid="stMetricValue"] {{ font-size: 1.5rem; color: {t["text"]}; }}
    div[data-testid="stMetricLabel"] {{ color: {t["metric_label"]}; }}

    /* ── Dividers ───────────────────────────────── */
    hr {{ border-color: {t["divider"]} !important; }}

    /* ── Typography ─────────────────────────────── */
    .stMarkdown, .stMarkdown p, .stMarkdown li {{ color: {t["text"]}; }}
    .stMarkdown h1, .stMarkdown h2, .stMarkdown h3, .stMarkdown h4 {{ color: {t["text"]}; }}
    .stCaption, .stCaption p {{ color: {t["muted"]} !important; }}
    label {{ color: {t["text"]} !important; }}

    /* ── Inputs ─────────────────────────────────── */
    .stTextInput input {{
        background-color: {t["input_bg"]} !important;
        color: {t["text"]} !important;
        border-color: {t["input_border"]} !important;
    }}
    .stSelectbox [data-baseweb="select"] > div {{
        background-color: {t["input_bg"]};
        color: {t["text"]};
    }}
    </style>
    """, unsafe_allow_html=True)


def render_lifecycle_header(state: str):
    t = get_active_theme()
    states = ["INSTANTIATED", "ACQUIRING", "AUDITED", "DECIDED", "SEALED"]
    chips = ""
    for i, s in enumerate(states):
        s_idx = states.index(state) if state in states else 0
        if i < s_idx:
            color, txt = t["chip_done_bg"], t["chip_done_text"]
        elif i == s_idx:
            color, txt = t["accent"], t["chip_active_text"]
        else:
            color, txt = t["chip_future_bg"], t["chip_future_text"]
        weight = "bold" if i == s_idx else "normal"
        chips += (
            f'<span style="background:{color};color:{txt};padding:4px 12px;margin:2px;'
            f'border-radius:12px;font-size:0.7rem;font-weight:{weight};'
            f'display:inline-block;">{s}</span>'
        )
    st.markdown(
        f'<div style="text-align:center;margin-bottom:15px;">{chips}</div>',
        unsafe_allow_html=True,
    )
