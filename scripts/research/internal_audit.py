import os, sys, json, requests
import gemmi
from fpdf import FPDF

def audit():
    print("--- 🧠 BRAIN INTERNAL DNA AUDIT ---")
    
    # H-1 (accepted ruling, ported Round 17): the external narrative
    # provider and its SDK were removed from the tree. No provider
    # credential is probed or required.

    # 2. Physics Engine Dependency
    try:
        doc = gemmi.read_pdb_string("ATOM      1  N   MET A   1       0.000   0.000   0.000")
        print(f"✅ Gemmi Engine: Operational (v{gemmi.__version__})")
    except:
        print("❌ Gemmi Engine: Failed to parse string.")

    # 2. PDF Engine Dependency
    try:
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("Helvetica", size=12)
        pdf.cell(200, 10, txt="Audit", ln=True)
        out = pdf.output(dest='S')
        print(f"✅ FPDF2 Engine: Operational ({len(out)} bytes generated)")
    except Exception as e:
        print(f"❌ FPDF2 Engine: Failed. {e}")

    # 4. NKG Filesystem
    nkg_path = "/app/nkg/piu_moat.jsonl"
    if os.access(os.path.dirname(nkg_path), os.W_OK):
        print(f"✅ NKG Filesystem: Writable at {os.path.dirname(nkg_path)}")
    else:
        print(f"❌ NKG Filesystem: NO WRITE ACCESS to {os.path.dirname(nkg_path)}")

if __name__ == "__main__":
    audit()
