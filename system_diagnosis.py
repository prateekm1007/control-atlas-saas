import os, sys, requests, json
import google.generativeai as genai

def diagnose():
    print("--- 🛡️ TOSCANINI FORENSIC DIAGNOSIS ---")
    
    # 1. API Key Check — key is read from the environment only, never hardcoded.
    key = os.environ.get("GEMINI_API_KEY")
    if not key:
        print("❌ Gemini API: GEMINI_API_KEY is not set. Export it (without printing it) and retry.")
        sys.exit(1)
    try:
        genai.configure(api_key=key)
        models = [m.name for m in genai.list_models()]
        print(f"✅ Gemini API: Authorized. Found {len(models)} models.")
    except Exception as e:
        print(f"❌ Gemini API: Failed. {str(e)}")

    # 2. AFDB Connectivity
    try:
        res = requests.get("https://alphafold.ebi.ac.uk/api/prediction/P01308", timeout=10)
        if res.status_code == 200:
            print(f"✅ AFDB API: Reachable. Found {len(res.json())} entries for P01308.")
    except:
        print("❌ AFDB API: Unreachable.")

    # 3. Dependency Check
    deps = ['gemmi', 'fpdf', 'numpy', 'fastapi']
    for d in deps:
        try:
            __import__(d)
            print(f"✅ Dependency: {d} found.")
        except:
            print(f"❌ Dependency: {d} MISSING.")

if __name__ == "__main__":
    diagnose()
