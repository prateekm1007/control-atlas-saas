"""
Confidence Visualization Module
Generates matplotlib-based confidence bar charts for PDF embedding.
"""
import io
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server use
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def generate_confidence_bar(coverage_pct: float, det_passed: int, det_total: int) -> bytes:
    """
    Generate a horizontal confidence bar showing coverage level.
    Returns PNG bytes for PDF embedding.
    """
    fig, ax = plt.subplots(figsize=(6, 0.8))
    
    # Color based on coverage threshold
    if coverage_pct >= 70:
        bar_color = '#2E7D32'  # Green - sufficient coverage
        label = 'SUFFICIENT'
    elif coverage_pct >= 50:
        bar_color = '#F57C00'  # Orange - marginal
        label = 'MARGINAL'
    else:
        bar_color = '#757575'  # Grey - insufficient
        label = 'INSUFFICIENT'
    
    # Background bar (100%)
    ax.barh(0, 100, height=0.6, color='#E0E0E0', edgecolor='none')
    
    # Coverage bar
    ax.barh(0, coverage_pct, height=0.6, color=bar_color, edgecolor='none')
    
    # 70% threshold line
    ax.axvline(x=70, color='#D32F2F', linestyle='--', linewidth=2, label='70% Threshold')
    
    # Coverage percentage text
    ax.text(coverage_pct / 2, 0, f'{coverage_pct:.1f}%', 
            ha='center', va='center', fontsize=12, fontweight='bold', color='white')
    
    # Label on right
    ax.text(102, 0, label, ha='left', va='center', fontsize=10, fontweight='bold', color=bar_color)
    
    # Styling
    ax.set_xlim(0, 130)
    ax.set_ylim(-0.5, 0.5)
    ax.axis('off')
    ax.set_title('Coverage Distribution vs. 70% Threshold', fontsize=11, fontweight='bold', loc='left', pad=10)
    
    plt.tight_layout()
    
    # Save to bytes
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def generate_deterministic_gauge(det_passed: int, det_total: int) -> bytes:
    """
    Generate a simple gauge showing deterministic law compliance.
    Returns PNG bytes for PDF embedding.
    """
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    
    passed_pct = (det_passed / det_total * 100) if det_total > 0 else 0
    failed = det_total - det_passed
    
    # Pie chart as gauge
    if failed > 0:
        colors = ['#2E7D32', '#D32F2F']  # Green for pass, red for fail
        sizes = [det_passed, failed]
        labels = [f'{det_passed} PASS', f'{failed} FAIL']
    else:
        colors = ['#2E7D32']
        sizes = [det_passed]
        labels = [f'{det_passed} PASS']
    
    wedges, texts = ax.pie(sizes, colors=colors, startangle=90,
                           wedgeprops=dict(width=0.4, edgecolor='white'))
    
    # Center text
    ax.text(0, 0, f'{det_passed}/{det_total}', ha='center', va='center',
            fontsize=16, fontweight='bold')
    
    ax.set_title('Deterministic\nCompliance', fontsize=10, fontweight='bold')
    
    # Legend
    patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, labels)]
    ax.legend(handles=patches, loc='lower center', fontsize=8, 
              bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
    
    plt.tight_layout()
    
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    buf.seek(0)
    return buf.read()
