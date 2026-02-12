import json
from pathlib import Path

SNAPSHOT_DIR = Path("entries/090_negative_knowledge/snapshots")
OUT_PATH = Path("entries/091_enforcement/manifold_law_index_v2_3.json")

def build_mli():
    snaps = sorted(SNAPSHOT_DIR.glob("nkg_v2_3_*_annotated.json"))
    if not snaps:
        raise RuntimeError("No annotated v2.3 snapshot found")

    snap = json.loads(snaps[-1].read_text())
    laws = []

    for e in snap.get("nkg", []):
        if (
            e.get("tier") in ("TIER_1", "TIER_2")
            and e.get("requires_geometry") is True
        ):
            laws.append({
                "entry_id": e.get("entry_id"),
                "target": e.get("target"),
                "generator": e.get("generator"),
                "failure_class": e.get("failure_class"),
                "enforcement_stage": e.get("enforcement_stage"),
                "summary": e.get("summary"),
                "notes": e.get("notes"),
            })

    out = {
        "meta": {
            "source_snapshot": snaps[-1].name,
            "version": "v2.3-manifold-law",
            "authority": "READ_ONLY_DERIVED",
            "law_count": len(laws),
        },
        "laws": laws,
    }

    OUT_PATH.write_text(json.dumps(out, indent=2))

    print("✅ Manifold Law Index built")
    print(f"   Snapshot: {snaps[-1].name}")
    print(f"   Laws: {len(laws)}")
    print(f"   Output: {OUT_PATH}")

if __name__ == "__main__":
    build_mli()
