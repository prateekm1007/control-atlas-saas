import json
from pathlib import Path

MDI_PATH = Path("entries/094_mdi/mdi_v2_9.json")

def load_mdi():
    if not MDI_PATH.exists():
        return []

    data = json.loads(MDI_PATH.read_text())
    return data.get("locked_doctrines", [])

if __name__ == "__main__":
    doctrines = load_mdi()
    print(f"✅ Loaded {len(doctrines)} locked doctrines")
    for d in doctrines:
        print(f" - {d['doctrine_id']}: {d['title']}")
