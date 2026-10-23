"""
Reliability-coverage visualization for PDF embedding.

This bar displays the single measured LAW-105 residue-coverage percentage.
It deliberately does not invent pLDDT bin proportions.
"""
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def coverage_segments(coverage_pct: float):
    """Return the measured covered and below-coverage percentages."""
    covered = min(max(float(coverage_pct), 0.0), 100.0)
    return covered, 100.0 - covered


def generate_confidence_bar(coverage_pct: float, det_passed: int, det_total: int) -> bytes:
    """
    Full-width reliability-coverage bar with the 70% minimum marker.
    Returns PNG bytes for PDF embedding.
    """
    fig, ax = plt.subplots(figsize=(10, 1.6))
    fig.patch.set_facecolor('white')

    covered, below = coverage_segments(coverage_pct)
    colors = ['#2E7D32', '#D32F2F'] if covered >= 70.0 else ['#D32F2F', '#9E9E9E']

    ax.barh(0, covered, height=0.55, color=colors[0],
            edgecolor='white', linewidth=1.5)
    ax.barh(0, below, left=covered, height=0.55, color=colors[1],
            edgecolor='white', linewidth=1.5)

    for start, width in ((0.0, covered), (covered, below)):
        if width >= 3:
            ax.text(start + width / 2, 0, "%.0f%%" % width,
                    ha='center', va='center', fontsize=10,
                    fontweight='bold', color='white', clip_on=True)

    ax.axvline(x=70, color='#D32F2F', linestyle='--',
               linewidth=2.0, zorder=5)
    ax.text(70, 0.38, '70% Minimum',
            fontsize=8, color='#D32F2F', va='center',
            ha='center', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      edgecolor='#D32F2F', alpha=0.9))

    cov_color = '#2E7D32' if covered >= 70 else '#D32F2F'
    cov_label = 'SUFFICIENT' if covered >= 70 else 'INSUFFICIENT'
    ax.text(101, 0, '%.1f%%\n%s' % (covered, cov_label),
            ha='left', va='center', fontsize=9,
            fontweight='bold', color=cov_color,
            linespacing=1.4)

    legend_patches = [
        mpatches.Patch(color=colors[0], label='Residues covered'),
        mpatches.Patch(color=colors[1], label='Residues below coverage'),
    ]
    ax.legend(handles=legend_patches,
              loc='upper center',
              bbox_to_anchor=(0.42, -0.25),
              ncol=2, fontsize=8, frameon=False,
              handlelength=1.5, handletextpad=0.5,
              columnspacing=1.5)

    ax.set_xlim(0, 130)
    ax.set_ylim(-0.5, 0.65)
    ax.axis('off')

    plt.title('Measured Residue Coverage (LAW-105)',
              fontsize=11, fontweight='bold', loc='left', pad=8)
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
