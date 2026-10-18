import os, sys, json, requests
import google.generativeai as genai
import gemmi
from fpdf import FPDF

def audit():
    print("--- 🧠 BRAIN INTERNAL DNA AUDIT ---")
    
    # 1. API Supply Chain
    key = os.getenv("GEMINI_API_KEY")
    if not key or key == "dummy":
        print("❌ GEMINI_API_KEY: Missing or Placeholder.")
    else:
        try:
            genai.configure(api_key=key)
            models = [m.name for m in genai.list_models()]
            print(f"✅ Gemini API: Authorized. {len(models)} models available.")
        except Exception as e:
            print(f"❌ Gemini API: Authorization Failed. {e}")

    # 2. Physics Engine Dependency
    try:
        doc = gemmi.read_pdb_string("ATOM      1  N   MET A   1       0.000   0.000   0.000")
        print(f"✅ Gemmi Engine: Operational (v{gemmi.__version__})")
    except:
        print("❌ Gemmi Engine: Failed to parse string.")

    # 3. PDF Engine Dependency
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
