import os, sys, requests, json

def diagnose():
    print("--- 🛡️ TOSCANINI FORENSIC DIAGNOSIS ---")

    # H-1 (Round 14): the external narrative provider and its SDK were
    # removed from the tree (owner instruction; auditor direction 3,
    # P0 stack). This diagnostic no longer probes that provider.

    # 1. AFDB Connectivity
    try:
        res = requests.get("https://alphafold.ebi.ac.uk/api/prediction/P01308", timeout=10)
        if res.status_code == 200:
            print(f"✅ AFDB API: Reachable. Found {len(res.json())} entries for P01308.")
    except:
        print("❌ AFDB API: Unreachable.")

    # 2. Dependency Check
    deps = ['gemmi', 'fpdf', 'numpy', 'fastapi']
    for d in deps:
        try:
            __import__(d)
            print(f"✅ Dependency: {d} found.")
        except:
            print(f"❌ Dependency: {d} MISSING.")

if __name__ == "__main__":
    diagnose()
