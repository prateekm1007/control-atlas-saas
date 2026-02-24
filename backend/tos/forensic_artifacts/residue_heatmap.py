"""
Residue-Level Confidence Heatmap (v1)
1D strip + line trace for PDF Page 4. Palette locked to structure_viz.py.
"""
import io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

BAND_COLORS = {
    "Very High": "#0066CC", "High": "#4DA6FF",
    "Medium": "#FFB84D", "Low": "#CC0000",
}
CMAP = ListedColormap([
    BAND_COLORS["Low"], BAND_COLORS["Medium"],
    BAND_COLORS["High"], BAND_COLORS["Very High"],
])
NORM = BoundaryNorm([0, 50, 70, 90, 100], CMAP.N)


def generate_residue_heatmap_png(plddt_scores, title="Per-Residue Confidence (pLDDT)",
                                  figwidth=7.0, figheight=1.8, dpi=200):
    """Return PNG bytes of the residue-level confidence heatmap."""
    n = len(plddt_scores)
    if n == 0:
        return _empty_png(figwidth, dpi)

    scores = np.clip(np.array(plddt_scores, dtype=np.float64), 0, 100)
    res_ids = np.arange(1, n + 1)

    fig, (ax_t, ax_s) = plt.subplots(
        2, 1, figsize=(figwidth, figheight),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05}, dpi=dpi)
    fig.patch.set_facecolor("white")

    # Strip
    ax_s.imshow(scores.reshape(1, -1), aspect="auto", cmap=CMAP, norm=NORM,
                extent=[0.5, n + 0.5, 0, 1], interpolation="nearest")
    ax_s.set_xlim(0.5, n + 0.5)
    ax_s.set_yticks([])
    ax_s.set_xlabel("Residue Index", fontsize=7, color="#333")
    step = 5 if n <= 50 else 20 if n <= 200 else 50 if n <= 500 else 100
    ticks = np.arange(step, n + 1, step)
    ticks = np.insert(ticks, 0, 1)
    ax_s.set_xticks(ticks)
    ax_s.tick_params(axis="x", labelsize=6)
    for sp in ["top", "right", "left"]:
        ax_s.spines[sp].set_visible(False)
    ax_s.spines["bottom"].set_color("#E0E0E0")

    # Trace
    ax_t.fill_between(res_ids, scores, alpha=0.06, color="#1A1A1A")
    ax_t.plot(res_ids, scores, linewidth=0.6, color="#1A1A1A", alpha=0.85)
    for th in [50, 70, 90]:
        ax_t.axhline(y=th, linewidth=0.4, linestyle="--", color="#E0E0E0", zorder=0)
        ax_t.text(n + 0.8, th, str(th), fontsize=5, color="#999", va="center")
    ax_t.set_xlim(0.5, n + 0.5)
    ax_t.set_ylim(0, 105)
    ax_t.set_ylabel("pLDDT", fontsize=7, color="#333")
    ax_t.set_title(title, fontsize=9, color="#333", pad=6)
    ax_t.tick_params(axis="both", labelsize=6)
    ax_t.set_xticks([])
    ax_t.set_facecolor("white")
    ax_t.spines["bottom"].set_visible(False)
    for sp in ["top", "right"]:
        ax_t.spines[sp].set_visible(False)
    ax_t.spines["left"].set_color("#E0E0E0")

    legend = [Patch(facecolor=BAND_COLORS[k], label=k) for k in BAND_COLORS]
    ax_t.legend(handles=legend, loc="upper right", fontsize=5.5,
                frameon=True, framealpha=0.9, edgecolor="#E0E0E0", ncol=4)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def _empty_png(figwidth, dpi):
    fig, ax = plt.subplots(figsize=(figwidth, 0.6), dpi=dpi)
    ax.text(0.5, 0.5, "Per-residue confidence data unavailable",
            ha="center", va="center", fontsize=9, color="#999", transform=ax.transAxes)
    ax.axis("off")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return buf.read()
