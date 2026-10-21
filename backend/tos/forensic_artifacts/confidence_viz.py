"""
Confidence Visualization Module v3
Institutional-grade pLDDT distribution bar for PDF embedding.
"""
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def generate_confidence_bar(bin_counts, coverage_pct: float, det_passed: int, det_total: int) -> bytes:
    """
    pLDDT bin histogram for PDF embedding, computed from MEASURED
    per-residue pLDDT observations (previous version fabricated bin
    proportions from aggregate coverage via fixed multipliers).

    Args:
        bin_counts: measured per-residue counts per band
            {"very_high": N, "high": N, "low": N} from the structure's
            own pLDDT/B-factor observations, or None when unmeasurable
            (the bar then renders an accurately labelled atom-coverage
            indicator instead of synthetic proportions).
        coverage_pct: atom coverage percentage (OD-04 basis).
        det_passed / det_total: deterministic compliance duo.
    """
    fig, ax = plt.subplots(figsize=(10, 1.6))
    fig.patch.set_facecolor('white')

    bins = [
        ("Very High (>=90)",  "#0050A0", 90, 100),
        ("High (70-90)",      "#4DA6FF", 70,  90),
        ("Low (<70)",         "#CC2200",  0,  70),
    ]

    if bin_counts and sum(bin_counts.values()) > 0:
        measured_total = float(sum(bin_counts.values()))
        fracs = [bin_counts.get(key, 0) / measured_total * 100.0
                 for key in ("very_high", "high", "low")]
    else:
        # No measured distributions available: do NOT generate synthetic
        # bin proportions — swap to the accurate coverage indicator.
        return _atom_coverage_indicator(coverage_pct)

    left = 0.0
    for (label, color, lo, hi), frac in zip(bins, fracs):
        if frac < 0.3:
            left += frac
            continue
        ax.barh(0, frac, left=left, height=0.55,
                color=color, edgecolor='white', linewidth=1.5)
        if frac >= 3:
            ax.text(left + frac / 2, 0,
                    "%.0f%%" % frac,
                    ha='center', va='center',
                    fontsize=10, fontweight='bold', color='white',
                    clip_on=True)
        left += frac

    ax.axvline(x=70, color='#D32F2F', linestyle='--',
               linewidth=2.0, zorder=5)
    ax.text(70, 0.38, '70% Threshold',
            fontsize=8, color='#D32F2F', va='center',
            ha='center', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor='#D32F2F', alpha=0.9))

    cov_color = '#2E7D32' if coverage_pct >= 70 else '#D32F2F'
    cov_label = 'SUFFICIENT' if coverage_pct >= 70 else 'INSUFFICIENT'
    ax.text(101, 0, '%.1f%%\n%s' % (coverage_pct, cov_label),
            ha='left', va='center', fontsize=9,
            fontweight='bold', color=cov_color,
            linespacing=1.4)

    legend_patches = [
        mpatches.Patch(color=c, label=lbl)
        for lbl, c, _, _ in bins
    ]
    ax.legend(handles=legend_patches,
              loc='upper center',
              bbox_to_anchor=(0.42, -0.25),
              ncol=3, fontsize=8, frameon=False,
              handlelength=1.5, handletextpad=0.5,
              columnspacing=1.5)

    ax.set_xlim(0, 130)
    ax.set_ylim(-0.5, 0.65)
    ax.axis('off')

    plt.title('pLDDT Confidence Distribution (measured, per residue)',
              fontsize=11, fontweight='bold', loc='left', pad=8)
    plt.tight_layout(pad=0.4)

    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=160,
                bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def _atom_coverage_indicator(coverage_pct: float) -> bytes:
    """Accurately labelled atom-coverage bar (no synthetic pLDDT bins).

    Used when per-residue pLDDT observations are unavailable: the segment
    shows the measured atom coverage share, labelled as such.
    """
    fig, ax = plt.subplots(figsize=(10, 1.6))
    fig.patch.set_facecolor('white')

    measured = min(max(coverage_pct, 0.0), 100.0)
    ax.barh(0, measured, left=0.0, height=0.55,
            color='#0050A0' if measured >= 70 else '#CC2200',
            edgecolor='white', linewidth=1.5)
    ax.barh(0, 100.0 - measured, left=measured, height=0.55,
            color='#E0E0E0', edgecolor='white', linewidth=1.5)

    ax.axvline(x=70, color='#D32F2F', linestyle='--',
               linewidth=2.0, zorder=5)
    ax.text(70, 0.38, '70% Threshold',
            fontsize=8, color='#D32F2F', va='center',
            ha='center', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor='#D32F2F', alpha=0.9))

    cov_color = '#2E7D32' if measured >= 70 else '#D32F2F'
    cov_label = 'SUFFICIENT' if measured >= 70 else 'INSUFFICIENT'
    ax.text(101, 0, '%.1f%%\n%s' % (measured, cov_label),
            ha='left', va='center', fontsize=9,
            fontweight='bold', color=cov_color,
            linespacing=1.4)

    ax.set_xlim(0, 130)
    ax.set_ylim(-0.5, 0.65)
    ax.axis('off')

    plt.title('Atom coverage (share of atoms in pLDDT >= 70 residues) - '
              'per-residue pLDDT distribution unavailable',
              fontsize=9, fontweight='bold', loc='left', pad=8)
    plt.tight_layout(pad=0.4)

    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=160,
                bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def generate_deterministic_gauge(det_passed: int, det_total: int) -> bytes:
    """Donut gauge showing deterministic law compliance."""
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    failed = det_total - det_passed
    if failed > 0:
        colors = ['#2E7D32', '#D32F2F']
        sizes = [det_passed, failed]
        labels = ['%d PASS' % det_passed, '%d FAIL' % failed]
    else:
        colors = ['#2E7D32']
        sizes = [det_passed]
        labels = ['%d PASS' % det_passed]

    wedges, texts = ax.pie(sizes, colors=colors, startangle=90,
                           wedgeprops=dict(width=0.4, edgecolor='white'))
    ax.text(0, 0, '%d/%d' % (det_passed, det_total),
            ha='center', va='center', fontsize=16, fontweight='bold')
    ax.set_title('Deterministic\nCompliance', fontsize=10, fontweight='bold')
    patches = [mpatches.Patch(color=c, label=l)
               for c, l in zip(colors, labels)]
    ax.legend(handles=patches, loc='lower center', fontsize=8,
              bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    buf.seek(0)
    return buf.read()
