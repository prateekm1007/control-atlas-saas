#!/usr/bin/env python3
"""U9 (Round-11): MODRES census over the v23 40-target experimental cohort.

Question: do any of the 40 experimental coordinate files in the recorded
v23 cohort carry MODRES (modified residue) annotation records?

Method:
  - cohort + per-file digest binding taken from the Round-9/10 item-level
    table (n7_od16_40target_table.json), which records the pinned
    dataset_v23.py blob SHA-1 it was derived from (verified below);
  - each local coordinate file's SHA-256 is verified against the table's
    c5_coordinate_sha256_fresh BEFORE counting (provenance binding);
  - counts PDB `MODRES ` records and, for context, HETATM occurrences of
    the engine's MODIFIED_RESIDUE_PARENTS names (MSE, SEP, TPO, ...).

Output: JSON results + this stdout log. Documentation only — no cohort,
code, or threshold change (U9).
"""
import hashlib
import json
import sys
from pathlib import Path

TABLE = Path("/home/z/my-project/work/audit-evidence/n7_od16_40target_table.json")
COORDS = Path("/home/z/my-project/work/audit-evidence/coordinates")
OUT_JSON = Path("/home/z/my-project/work/r11-outputs/evidence/modres-census/"
                "modres_census_results.json")

# Mirror of backend/tos/ingestion/processor.py MODIFIED_RESIDUE_PARENTS keys
# (science branch blob 3a535837... @ ebc1fd79) — for the HETATM context count.
MODIFIED_RESIDUE_NAMES = [
    "MSE", "SEC", "PYL", "HYP", "CME", "CSO", "CSD", "OCS", "CSX", "CAS",
    "CSM", "MLY", "M3L", "KCX", "SEP", "TPO", "PTR", "PCA", "NLE",
]
PINNED_DATASET_BLOB = "3a5358374b2037b9475bd2d494754944485da5ac"


def main() -> int:
    table = json.loads(TABLE.read_text())
    blob_ok = table["dataset_v23_blob_sha1"] == PINNED_DATASET_BLOB
    rows = table["rows"]
    print(f"table dataset_v23 blob: {table['dataset_v23_blob_sha1']} "
          f"({'MATCHES pinned' if blob_ok else 'MISMATCH'})")
    if not blob_ok:
        print("FAIL: table not bound to the pinned dataset blob")
        return 1
    print(f"cohort rows: {len(rows)}")

    results = []
    total_modres = 0
    total_het_mod = 0
    digest_failures = 0
    for r in rows:
        pid = r["pdb_id"]
        expected = r["c5_coordinate_sha256_fresh"]
        path = COORDS / f"{pid}.pdb"
        if not path.exists():
            results.append({"pdb_id": pid, "error": "coordinate file missing"})
            digest_failures += 1
            continue
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        bound = digest == expected
        if not bound:
            digest_failures += 1
        text = data.decode("utf-8", errors="replace")
        modres_lines = [l for l in text.splitlines() if l.startswith("MODRES ")]
        het_counts = {}
        for l in text.splitlines():
            if l.startswith("HETATM"):
                name = l[17:20].strip()
                if name in MODIFIED_RESIDUE_NAMES:
                    het_counts[name] = het_counts.get(name, 0) + 1
        total_modres += len(modres_lines)
        total_het_mod += sum(het_counts.values())
        results.append({
            "pdb_id": pid,
            "sha256": digest,
            "sha256_matches_table": bound,
            "modres_record_count": len(modres_lines),
            "modres_records": modres_lines[:10],
            "hetatm_modified_residue_atoms": het_counts,
        })
        flag = "" if bound else "  <-- DIGEST MISMATCH"
        print(f"{pid}: modres={len(modres_lines)} "
              f"hetatm_modified_atoms={sum(het_counts.values())} "
              f"digest_bound={bound}{flag}")

    print()
    print(f"TOTAL: {len(rows)} targets, {total_modres} MODRES records, "
          f"{total_het_mod} HETATM atoms of modified-residue names, "
          f"{digest_failures} digest failures")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps({
        "census": "MODRES census of the v23 40-target experimental cohort",
        "purpose": "U9 Round-11 documentation-only census",
        "table_source": str(TABLE),
        "dataset_v23_blob_sha1": table["dataset_v23_blob_sha1"],
        "coordinate_digest_field": "c5_coordinate_sha256_fresh",
        "totals": {
            "targets": len(rows),
            "modres_records": total_modres,
            "hetatm_modified_residue_atoms": total_het_mod,
            "digest_failures": digest_failures,
        },
        "results": results,
    }, indent=2))
    print(f"written: {OUT_JSON}")
    return 0 if digest_failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
