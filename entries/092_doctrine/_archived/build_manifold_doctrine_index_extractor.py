import json
from pathlib import Path

SNAPSHOT_DIR = Path("entries/090_negative_knowledge/snapshots")
OUT_PATH = Path("entries/092_doctrine/manifold_doctrine_index_v2_3.json")

def build_mdi():
    snaps = sorted(SNAPSHOT_DIR.glob("nkg_v2_3_*_annotated.json"))
    if not snaps:
        raise RuntimeError("No annotated v2.3 snapshot found")

    snap = json.loads(snaps[-1].read_text())
    doctrines = []

    for e in snap.get("nkg", []):
        if (
            e.get("requires_geometry") is True
            and e.get("generator") is not None
            and e.get("summary") is not None
        ):
            doctrines.append({
                "entry_id": e.get("entry_id"),
                "target": e.get("target"),
                "generator": e.get("generator"),
                "failure_class": e.get("failure_class"),
                "summary": e.get("summary"),
                "notes": e.get("notes"),
                "confidence": e.get("confidence", "UNSPECIFIED"),
                "consultation_hint": "INTERPRETIVE_REASONING_REQUIRED"
            })

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    out = {
        "meta": {
            "source_snapshot": snaps[-1].name,
            "version": "v2.3-doctrine",
            "authority": "READ_ONLY_DERIVED",
            "doctrine_count": len(doctrines)
        },
        "doctrines": doctrines
    }

    OUT_PATH.write_text(json.dumps(out, indent=2))

    print("✅ Manifold Doctrine Index built")
    print(f"   Snapshot: {snaps[-1].name}")
    print(f"   Doctrines: {len(doctrines)}")
    print(f"   Output: {OUT_PATH}")

if __name__ == "__main__":
    build_mdi()
