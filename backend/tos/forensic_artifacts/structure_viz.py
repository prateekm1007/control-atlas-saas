"""
3D Structure Visualization Module
Generates matplotlib-based CA trace scatter plots colored by confidence.
"""
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def generate_structure_render(coords: list, confidences: list, title: str = "Structure",
                              metric_label="pLDDT") -> bytes:
    """
    Generate 3D scatter plot of CA atoms, optionally colored by pLDDT.

    coords: list of (x, y, z) tuples
    confidences: list of measured pLDDT values; ``None``/non-finite values
                 are rendered neutral and are never replaced with a score.
    metric_label: "pLDDT" enables the pLDDT colour bands. Any other value
                  renders the trace without a confidence colour scale.
    Returns: PNG bytes
    """
    if not coords or len(coords) < 3:
        # Not enough atoms for visualization
        return generate_placeholder()
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Extract coordinates
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    
    finite = [c for c in confidences
              if c is not None and np.isfinite(float(c))]
    colors = []
    if metric_label != "pLDDT" or not finite:
        colors = ['#808080'] * len(coords)
    else:
        for i in range(len(coords)):
            conf = confidences[i] if i < len(confidences) else None
            if conf is None or not np.isfinite(float(conf)):
                colors.append('#808080')  # missing observation, never fabricated
            elif conf >= 90:
                colors.append('#0066CC')  # Dark blue - very high confidence
            elif conf >= 70:
                colors.append('#4DA6FF')  # Light blue - high confidence
            elif conf >= 50:
                colors.append('#FFB84D')  # Orange - medium confidence
            else:
                colors.append('#CC0000')  # Red - low confidence
    
    # Plot atoms as spheres
    ax.scatter(xs, ys, zs, c=colors, s=50, alpha=0.8, edgecolors='white', linewidth=0.5)
    
    # Connect with CA trace line
    ax.plot(xs, ys, zs, color='grey', alpha=0.3, linewidth=1)
    
    # Styling
    ax.set_xlabel('X (Å)', fontsize=10)
    ax.set_ylabel('Y (Å)', fontsize=10)
    ax.set_zlabel('Z (Å)', fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # Legend only when the observed metric is pLDDT
    if metric_label == "pLDDT" and finite:
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#0066CC', label='Very High (≥90)'),
            Patch(facecolor='#4DA6FF', label='High (70-90)'),
            Patch(facecolor='#FFB84D', label='Medium (50-70)'),
            Patch(facecolor='#CC0000', label='Low (<50)'),
            Patch(facecolor='#808080', label='Unavailable'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    # Rotate for better view
    ax.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    
    # Save to bytes
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def generate_placeholder() -> bytes:
    """Generate placeholder image when structure data insufficient."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(0.5, 0.5, 'Structure visualization unavailable\n(insufficient coordinate data)',
            ha='center', va='center', fontsize=14, color='grey',
            bbox=dict(boxstyle='round', facecolor='#F0F0F0', alpha=0.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    buf.seek(0)
    return buf.read()
