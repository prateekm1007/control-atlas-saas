#!/usr/bin/env python3
"""MODRES census v2 — v23 40-target experimental cohort (Round-12 fix).

Question: do any of the 40 experimental coordinate files in the recorded
v23 cohort carry MODRES (modified residue) annotation records?

v2 portability fixes (Round-12 direction 4, Round-13 delivery):
  - the input artifact (the Round-9/10 item-level table) is read from a
    REPOSITORY-RELATIVE path, committed alongside this script, and its
    SHA-256 is recorded in the results;
  - coordinate bytes are fetched from the canonical public source
    (https://files.rcsb.org/download/<pdb_id>.pdb) and every file's
    SHA-256 is verified against the table's c5_coordinate_sha256_fresh
    BEFORE counting, so the census is bound to the exact recorded bytes
    without depending on any coder's local filesystem;
  - no absolute local paths appear anywhere in the outputs.

Method (unchanged from v1): counts PDB `MODRES ` records and, for
context, HETATM occurrences of the engine's MODIFIED_RESIDUE_PARENTS
names (MSE, SEP, TPO, ...).

Documentation only — no cohort, code, or threshold change.
"""
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLE = HERE / "n7_od16_40target_table.json"
OUT_JSON = HERE / "modres_census_results.json"
RCSB = "https://files.rcsb.org/download/{pdb_id}.pdb"
TABLE_REPO_RELATIVE = "docs/audit-evidence/2026-10-20/modres-census-v2/n7_od16_40target_table.json"

# Mirror of backend/tos/ingestion/processor.py MODIFIED_RESIDUE_PARENTS keys
# (dataset pinned at blob 3a535837... — the table records the same blob and
# is verified against it below) — for the HETATM context count.
MODIFIED_RESIDUE_NAMES = [
    "MSE", "SEC", "PYL", "HYP", "CME", "CSO", "CSD", "OCS", "CSX", "CAS",
    "CSM", "MLY", "M3L", "KCX", "SEP", "TPO", "PTR", "PCA", "NLE",
]
PINNED_DATASET_BLOB = "3a5358374b2037b9475bd2d494754944485da5ac"


def fetch(pdb_id: str) -> bytes:
    with urllib.request.urlopen(RCSB.format(pdb_id=pdb_id), timeout=60) as r:
        return r.read()


def main() -> int:
    table_bytes = TABLE.read_bytes()
    table_sha256 = hashlib.sha256(table_bytes).hexdigest()
    table = json.loads(table_bytes)
    blob_ok = table["dataset_v23_blob_sha1"] == PINNED_DATASET_BLOB
    print(f"input artifact: {TABLE_REPO_RELATIVE}")
    print(f"input artifact sha256: {table_sha256}")
    print(f"table dataset_v23 blob: {table['dataset_v23_blob_sha1']} "
          f"({'MATCHES pinned' if blob_ok else 'MISMATCH'})")
    if not blob_ok:
        print("FAIL: table not bound to the pinned dataset blob")
        return 1
    rows = table["rows"]
    print(f"cohort rows: {len(rows)}")

    results = []
    total_modres = 0
    total_het_mod = 0
    digest_failures = 0
    for r in rows:
        pid = r["pdb_id"]
        expected = r["c5_coordinate_sha256_fresh"]
        try:
            data = fetch(pid)
        except Exception as e:
            results.append({"pdb_id": pid, "error": f"fetch failed: {e}"})
            digest_failures += 1
            continue
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
            "source": f"https://files.rcsb.org/download/{pid}.pdb",
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

    OUT_JSON.write_text(json.dumps({
        "census": "MODRES census of the v23 40-target experimental cohort (v2)",
        "purpose": "documentation-only census; Round-12 portability fix",
        "input_artifact": {
            "path": TABLE_REPO_RELATIVE,
            "sha256": table_sha256,
            "rows": len(rows),
            "dataset_v23_blob_sha1": table["dataset_v23_blob_sha1"],
        },
        "coordinate_digest_field": "c5_coordinate_sha256_fresh",
        "coordinate_source": "https://files.rcsb.org/download/<pdb_id>.pdb "
                             "(canonical; per-target sha256 verified against "
                             "the table before counting)",
        "totals": {
            "targets": len(rows),
            "modres_records": total_modres,
            "hetatm_modified_residue_atoms": total_het_mod,
            "digest_failures": digest_failures,
        },
        "results": results,
    }, indent=2))
    print(f"written: {OUT_JSON.name} (repo-relative: "
          f"docs/audit-evidence/2026-10-20/modres-census-v2/)")
    return 0 if digest_failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
