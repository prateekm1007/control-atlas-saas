import streamlit as st
import json
from pathlib import Path
import tempfile
import os
import sys
import numpy as np
from typing import Dict, Any

# Import your tools
sys.path.insert(0, 'tools')
from native_audit import SovereignJudge
from kinetic_rescue import run_rescue

st.set_page_config(
    page_title="Sovereign Sieve",
    page_icon="🛡️",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🛡️ Sovereign Sieve")
st.markdown("**Falsify Before You Fabricate** - Physics-based validation for AI protein designs")

# Sidebar
st.sidebar.header("Upload Structure")
uploaded_file = st.sidebar.file_uploader(
    "Choose CIF/PDB file", 
    type=['cif', 'pdb'],
    help="Upload your Chai-1, RFdiffusion, or AlphaFold output"
)

target_chain = st.sidebar.text_input(
    "Target Chain ID", 
    value="A", 
    help="Chain ID for target protein (usually A)"
)

binder_chain = st.sidebar.text_input(
    "Binder Chain ID", 
    value="B", 
    help="Chain ID for peptide binder (usually B)"
)

run_minimization = st.sidebar.checkbox(
    "Run Energy Minimization", 
    value=True,
    help="Apply Amber14/OBC2 relaxation to resolve static clashes"
)

if st.sidebar.button("Falsify Structure", type="primary"):
    if uploaded_file is not None:
        # Save uploaded file
        temp_dir = Path(tempfile.mkdtemp())
        temp_cif = temp_dir / "input.cif"
        temp_cif.write_bytes(uploaded_file.read())
        
        with st.spinner("Running geometric audit..."):
            # Step 1: Geometric audit
            judge = SovereignJudge()
            audit_result = judge.calculate_contact_density(
                str(temp_cif),
                target_chain=target_chain,
                binder_chain=binder_chain
            )
        
        # Step 2: Energy minimization
        if run_minimization:
            with st.spinner("Running energy minimization..."):
                minimized_path = temp_dir / "minimized.pdb"
                rescue_result = run_rescue(str(temp_cif), str(minimized_path))
                
                # Re-audit minimized structure
                post_audit = judge.calculate_contact_density(
                    str(minimized_path),
                    target_chain=target_chain,
                    binder_chain=binder_chain
                )
                
                audit_result["pre_minimization"] = audit_result
                audit_result = post_audit
                audit_result["minimization"] = {
                    "pre_rho": rescue_result.get("pre", {}).get("rho", audit_result["rho"]),
                    "post_rho": audit_result["rho"],
                    "rho_loss": ((rescue_result.get("pre", {}).get("rho", 0) - audit_result["rho"]) / rescue_result.get("pre", {}).get("rho", 1)) * 100
                }
        else:
            audit_result["minimization"] = {"skipped": True}
        
        # Display results
        st.subheader("Falsification Report")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("Status", audit_result["status"])
            st.metric("Contact Density (ρ)", audit_result["rho"])
            st.metric("Min Distance", f"{audit_result['min_distance_A']:.2f} Å")
            st.metric("Clashes (< 2.5Å)", audit_result["clashes"])
        
        with col2:
            if "minimization" in audit_result and "skipped" not in audit_result["minimization"]:
                st.metric("ρ Loss (%)", f"{audit_result['minimization']['rho_loss']:.1f}%")
                if audit_result['minimization']['rho_loss'] < 5:
                    st.success("✅ Industrial Grade (stable)")
                else:
                    st.warning("⚠️  Metastable (repositioning)")
            else:
                st.info("No minimization performed")
        
        # Detailed report
        st.subheader("Detailed Metrics")
        st.json(audit_result)
        
        # Download
        report_json = json.dumps(audit_result, indent=2)
        st.download_button(
            label="Download Report (JSON)",
            data=report_json,
            file_name="sovereign_sieve_report.json",
            mime="application/json"
        )
        
        if run_minimization and Path(temp_dir / "minimized.pdb").exists():
            with open(temp_dir / "minimized.pdb", "rb") as f:
                st.download_button(
                    label="Download Minimized PDB",
                    data=f.read(),
                    file_name="minimized.pdb",
                    mime="application/octet-stream"
                )
        
        # Recommendation
        st.subheader("Recommendation")
        if audit_result["status"] == "GEOMETRIC_PASS" and audit_result.get("minimization", {}).get("rho_loss", 100) < 5:
            st.success("✅ **SOVEREIGN PASS** - Ready for wet-lab synthesis")
            st.markdown("""
            Your structure passes geometric audit and maintains stability during energy minimization.
            This is an industrial-grade lead ready for experimental validation.
            """)
        else:
            st.error("❌ **VETO** - Requires refinement")
            st.markdown("""
            Your structure has geometric clashes or significant contact loss during minimization.
            Recommended actions:
            - Adjust linker length/spacing
            - Modify warhead positioning
            - Try alternative helical prop (EAAAK vs AEAAK)
            - Reduce peptide length
            """)
        
    else:
        st.warning("Upload a CIF/PDB file to begin falsification")

# Sidebar info
st.sidebar.markdown("""
### How It Works
1. **Geometric Audit**: Calculates contact density and clashes between chains
2. **Energy Minimization**: Relaxes structure with Amber14/OBC2 (100 steps)
3. **Stability Check**: Measures contact retention (< 5% loss = Industrial Grade)

### Status Guide
- **SOVEREIGN PASS**: 0 clashes, <5% contact loss
- **VETO**: Clashes persist or >20% contact loss
- **Metastable**: Passes geometry but repositions during minimization

### About
Built on Control Atlas v5.0. Validates AI-generated protein structures for wet-lab readiness.
""")

if __name__ == "__main__":
    st.run()
