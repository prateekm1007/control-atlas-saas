"""
dashboard/usage_widget.py

Toscanini B.6 — Credits / Quota Sidebar Widget

Design:
  - Polls /usage via BACKEND_URL env var (set in docker-compose face service)
  - 30-second Streamlit cache — never blocks the UI
  - Graceful degradation: brain offline → warning, not crash
  - Tier-aware colour: free=amber, pro=blue, enterprise=green
  - No PDB coordinates stored — metadata display only
"""
from __future__ import annotations

import os
import time
from typing import Optional

import requests
import streamlit as st

# Match dashboard.py: uses BACKEND_URL, not BRAIN_URL
_BACKEND    = os.environ.get("BACKEND_URL", "http://brain:8000")
_API_KEY    = os.environ.get("TOSCANINI_API_KEY", "")
_CACHE_TTL  = int(os.environ.get("USAGE_CACHE_TTL", "30"))
_TIMEOUT    = 4   # seconds — never block UI on brain slowness

_TIER_CFG = {
    "free":       {"colour": "#F59E0B", "icon": "🔶", "label": "Free"},
    "pro":        {"colour": "#3B82F6", "icon": "💎", "label": "Pro"},
    "enterprise": {"colour": "#10B981", "icon": "🏛",  "label": "Enterprise"},
}
_TIER_FALLBACK = {"colour": "#6B7280", "icon": "❓", "label": "Unknown"}


@st.cache_data(ttl=_CACHE_TTL, show_spinner=False)
def _fetch_usage(api_key: str) -> Optional[dict]:
    """
    Fetch /usage for api_key. Cached _CACHE_TTL seconds.
    Returns None on any failure — caller renders graceful fallback.
    """
    if not api_key:
        return None
    try:
        resp = requests.get(
            f"{_BACKEND}/usage",
            headers={"X-API-Key": api_key},
            timeout=_TIMEOUT,
        )
        return resp.json() if resp.status_code == 200 else None
    except requests.exceptions.RequestException:
        return None


def render_usage_sidebar(api_key: str = "") -> None:
    """
    Render the Credits Remaining panel in the Streamlit sidebar.
    Call inside `with st.sidebar:` in dashboard.py.

    Parameters
    ----------
    api_key : str   Active API key. Falls back to TOSCANINI_API_KEY env var.
    """
    key = api_key or _API_KEY
    st.sidebar.markdown("---")
    st.sidebar.markdown("### ⚡ Account Status")

    data = _fetch_usage(key)

    if data is None:
        st.sidebar.warning(
            "⚠️ **Usage unavailable**\n\n"
            "Brain not responding. Your key is still active.\n"
            "Display recovers automatically.",
        )
        return

    tier_key = data.get("tier", "free").lower()
    cfg      = _TIER_CFG.get(tier_key, _TIER_FALLBACK)

    # Tier badge
    st.sidebar.markdown(
        f"<span style='color:{cfg['colour']};font-weight:bold;'>"
        f"{cfg['icon']} {cfg['label']} Tier</span>",
        unsafe_allow_html=True,
    )

    # Credit bar
    used  = int(data.get("credits_used",  0))
    total = int(data.get("credits_total", 1))
    rem   = max(0, total - used)
    fill  = rem / total if total else 0.0
    pct   = fill * 100

    st.sidebar.metric(
        label="Credits Remaining",
        value=f"{rem:,}",
        delta=f"{pct:.0f}% of {total:,}",
    )
    st.sidebar.progress(fill, text=f"{used:,} used / {total:,} total")

    # Quota detail
    with st.sidebar.expander("📊 Quota Details", expanded=False):
        q = data.get("quotas", {})
        for label, val in [
            ("Audits / day",      q.get("audits_per_day", "—")),
            ("GPU runs",          q.get("gpu_runs",        "—")),
            ("Batch max",         q.get("batch_max",       "—")),
            ("Max residues",      q.get("max_residues",    "—")),
            ("Max file (MB)",     q.get("max_file_mb",     "—")),
        ]:
            c1, c2 = st.sidebar.columns([3, 1])
            c1.caption(label)
            c2.caption(str(val))

    st.sidebar.caption(
        f"⏱ Refreshes every {_CACHE_TTL}s · "
        f"{time.strftime('%H:%M:%S')}"
    )
