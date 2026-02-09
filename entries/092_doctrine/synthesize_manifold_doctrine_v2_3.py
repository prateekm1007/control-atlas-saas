import json
from collections import defaultdict
from pathlib import Path

SNAPSHOT_DIR = Path("entries/090_negative_knowledge/snapshots")
OUT_PATH = Path("entries/092_doctrine/manifold_doctrine_index_v2_3.json")

def synthesize():
    snaps = sorted(SNAPSHOT_DIR.glob("nkg_v2_3_*_annotated.json"))
    if not snaps:
        raise RuntimeError("No annotated v2.3 snapshot found")

    snap = json.loads(snaps[-1].read_text())
    buckets = defaultdict(list)

    for e in snap.get("nkg", []):
        ctx = e.get("context") or {}
        if not ctx and e.get("forensics"):
            ctx = {"forensics_keys": sorted(e["forensics"].keys())}

        key = (
            e.get("failure_class"),
            e.get("failure_tag"),
            json.dumps(ctx, sort_keys=True),
        )

        buckets[key].append(e)

    doctrines = []
    for (failure_class, failure_tag, ctx_sig), entries in buckets.items():
        if len(entries) < 2:
            continue  # doctrine requires repetition

        avg_conf = sum(e.get("confidence", 0.5) for e in entries) / len(entries)
        max_sev = max(e.get("severity_rank", 0) for e in entries)

        doctrines.append({
            "doctrine_id": f"MDI-{len(doctrines)+1:03d}",
            "failure_class": failure_class,
            "failure_tag": failure_tag,
            "context_signature": ctx_sig,
            "supporting_events": len(entries),
            "avg_confidence": round(avg_conf, 3),
            "max_severity_rank": max_sev,
            "authority": "DERIVED_MANIFOLD_SYNTHESIS",
            "consultation_hint": "INTERPRETIVE_REASONING_REQUIRED"
        })

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    out = {
        "meta": {
            "source_snapshot": snaps[-1].name,
            "version": "v2.3-doctrine-synthesized",
            "doctrine_count": len(doctrines)
        },
        "doctrines": doctrines
    }

    OUT_PATH.write_text(json.dumps(out, indent=2))

    print("✅ Manifold Doctrine Index synthesized")
    print(f"   Snapshot: {snaps[-1].name}")
    print(f"   Doctrines: {len(doctrines)}")
    print(f"   Output: {OUT_PATH}")

if __name__ == "__main__":
    synthesize()
