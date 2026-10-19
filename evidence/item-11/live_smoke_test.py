"""Live smoke test of the metadata-first AFDB acquisition (P0 item 11).

Runs the REAL code path (tos.generation.afdb.fetch_prediction) against
the live API for a small set of accessions, including the two known
fail-closed cohort cases. No secrets, read-only public API.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from tos.generation.afdb import fetch_prediction, AFDBUnavailable, uniprot_entry_type

CASES = [
    "AF-Q5VSL9-F1",   # expected: acquired (v6, cif)
    "AF-P0DTC2-F1",   # expected: fail closed, migrated id recorded
    "AF-P62988-F1",   # expected: fail closed, UniProt Inactive recorded
    "AF-P69905-F1",   # expected: acquired (hemoglobin alpha, classic)
]

for case in CASES:
    print(f"=== {case} ===")
    try:
        pred = fetch_prediction(case)
        print(f"  ACQUIRED  entity={pred.model_entity_id} v{pred.latest_version}")
        print(f"  url={pred.source_url}")
        print(f"  fmt={pred.fmt}  bytes={len(pred.content)}  sha256={pred.sha256[:16]}...")
        print(f"  uniprot={pred.uniprot_accession}  globalMetric={pred.global_metric_value}")
        assert "model_v4" not in pred.source_url, "must not use dead v4 URLs"
        assert pred.content[:100].lstrip().startswith((b"data_", b"#")), "cif content"
    except AFDBUnavailable as exc:
        print(f"  FAIL-CLOSED  {exc}")
        if exc.record.get("migrated_identifier"):
            print(f"  migrated_identifier={exc.record['migrated_identifier']}")
        acc = case.replace("AF-", "").split("-F")[0]
        et = uniprot_entry_type(acc)
        print(f"  uniprot_entry_type={et}")
print("done")
