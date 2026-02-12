import gradio as gr
import json
from pathlib import Path
import tempfile
from typing import Dict, Any
import os

# Import your tools (adjust paths)
import sys
sys.path.insert(0, 'tools')
from native_audit import SovereignJudge
from kinetic_rescue import run_rescue

def falsify_structure(
    uploaded_file: Path,
    target_chain: str = "A",
    binder_chain: str = "B",
    run_minimization: bool = True
) -> Dict[str, Any]:
    """
    Falsify a structure: Geometric audit + optional energy minimization.
    """
    if not uploaded_file:
        return {"error": "No file uploaded"}
    
    # Save uploaded file
    temp_dir = Path(tempfile.mkdtemp())
    temp_cif = temp_dir / "input.cif"
    temp_cif.write_bytes(uploaded_file.read_bytes())
    
    try:
        # Step 1: Geometric audit
        judge = SovereignJudge()
        audit_result = judge.calculate_contact_density(
            str(temp_cif),
            target_chain=target_chain,
            binder_chain=binder_chain
        )
        
        # Step 2: Energy minimization (if requested)
        if run_minimization:
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
        
        # Step 3: Generate report
        report = {
            "status": audit_result["status"],
            "metrics": {
                "rho": audit_result["rho"],
                "min_distance": audit_result["min_distance"],
                "clashes": audit_result["clashes"],
                "target_atoms": audit_result["target_atoms"],
                "binder_atoms": audit_result["binder_atoms"]
            },
            "recommendation": "PASS" if audit_result["status"] == "PASS" else "VETO - Refine design"
        }
        
        # Save results
        report_path = temp_dir / "report.json"
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)
        
        return {
            "report": json.dumps(report, indent=2),
            "download": str(report_path),
            "status": report["status"]
        }
        
    except Exception as e:
        return {"error": str(e)}

# Gradio interface
iface = gr.Interface(
    fn=falsify_structure,
    inputs=[
        gr.File(label="Upload CIF/PDB file", file_types=[".cif", ".pdb"]),
        gr.Textbox(label="Target Chain ID", value="A", placeholder="A"),
        gr.Textbox(label="Binder Chain ID", value="B", placeholder="B"),
        gr.Checkbox(label="Run Energy Minimization", value=True)
    ],
    outputs=[
        gr.Textbox(label="Falsification Report (JSON)", lines=20),
        gr.File(label="Download Report + Minimized PDB")
    ],
    title="Sovereign Sieve: Physics-Based Structure Validation",
    description="Upload AI-generated protein structures. Get geometric audit and energy minimization. VETO or PASS for wet-lab readiness.",
    theme=gr.themes.Soft()
)

if __name__ == "__main__":
    iface.launch(share=True)
