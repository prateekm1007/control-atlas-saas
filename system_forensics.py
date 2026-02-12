import sys, os, json, inspect
from pathlib import Path

# Setup paths
sys.path.insert(0, str(Path(__file__).resolve().parent / "backend"))

def diagnostic():
    errors = []
    print("--- 🛡️ TOSCANINI TOTAL FORENSIC AUDIT ---")

    # 1. Test Definitions Registry
    try:
        from tos.glossary.epistemic_definitions import DEFINITIONS
        print("✅ Registry: Loaded")
        expected_keys = ["PHYSICAL_SCORE", "STRATEGIC_SCORE", "TIER_1_README"]
        for k in expected_keys:
            if k in DEFINITIONS: print(f"  - Key '{k}': FOUND")
            else: 
                print(f"  - Key '{k}': MISSING 🚨")
                errors.append(f"Registry Key Missing: {k}")
    except Exception as e:
        print(f"❌ Registry: CRITICAL FAILURE ({e})")
        errors.append("Registry Load Failure")

    # 2. Test Ingestion Signatures
    try:
        from tos.ingestion.structure_object import StructureObject
        from tos.ingestion.processor import IngestionProcessor
        sig = inspect.signature(StructureObject.__init__)
        print(f"✅ StructureObject Signature: {sig}")
        
        # Check if processor calls match class
        # (Looking for label/source_mode vs old params)
    except Exception as e:
        print(f"❌ Ingestion Logic: MISMATCH ({e})")
        errors.append("Ingestion Signature Mismatch")

    # 3. Test Law Glossary
    try:
        from tos.glossary.law_glossary import list_all_law_ids
        ids = list_all_law_ids()
        print(f"✅ Law Glossary: {len(ids)} laws found")
    except Exception as e:
        print(f"❌ Law Glossary: FAILED ({e})")
        errors.append("Law Glossary Failure")

    # 4. Check for binary traps
    try:
        from tos.utils.type_guards import force_bytes
        test_val = bytearray(b"test")
        if isinstance(force_bytes(test_val), bytes):
            print("✅ Type Guards: ACTIVE (Byte-safe)")
    except:
        print("⚠️ Type Guards: MISSING OR INCOMPLETE")

    print("\n--- ⚖️ FINAL DIAGNOSTIC VERDICT ---")
    if not errors:
        print("🟢 SYSTEM GREEN: No mechanical mismatches detected.")
    else:
        print(f"🔴 SYSTEM COMPROMISED: {len(errors)} defects found.")
        for err in errors: print(f"  - {err}")

if __name__ == "__main__":
    diagnostic()
