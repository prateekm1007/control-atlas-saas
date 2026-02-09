import json
from pathlib import Path

def verify():
    mdi_path = Path("entries/_094_mdi/mdi_v2_5.json")
    if not mdi_path.exists():
        print("❌ MDI ledger missing at entries/_094_mdi/")
        return

    data = json.loads(mdi_path.read_text())
    laws = len(data.get("locked_doctrines", []))

    verified = len(list(Path("ledger/motifs/verified").glob("*")))
    rejected = len(list(Path("ledger/motifs/rejected").glob("*")))

    print("=== TRACEABILITY AUDIT ===")
    print(f"MDI Laws Found      : {laws}")
    print(f"Verified Motifs     : {verified}")
    print(f"Rejected Motifs     : {rejected}")

    if laws < 50:
        print("⚠️ WARNING: Law count unusually low")
    else:
        print("✅ Law corpus sufficient for claims")

if __name__ == "__main__":
    verify()
