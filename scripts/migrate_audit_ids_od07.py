#!/usr/bin/env python3
"""OD-07 migration: re-key stored audit records to the canonical identity.

Round 14, direction 2 (auditor ruling, reaffirmed): ``audit_id`` is the
SHA-256 of one canonical serialization binding the coordinate hash, the
canon hash, and the matrix hash. Records created under the earlier
schemes are re-keyed by THIS EXPLICIT OPERATOR SCRIPT — never by a read
path (the legacy lookup is read-only; no writes on GET).

What it migrates (files under ``$TOSCANINI_DATA_DIR/audits/*.json``):

1. WP-5-era records stored under the bare ``sha256(coord_hash)``
   filename: re-keyed to the canonical combined ID derived from the
   record's own ``provenance.hash`` + ``governance_fingerprint``; the
   record's ``governance.audit_id`` is updated to match.
2. Pre-WP-5 records stored under 8-char filenames that carry enough
   metadata (provenance hash + fingerprint): same treatment, and the
   8-char name is recorded as ``governance.legacy_audit_id`` when absent.

Records lacking the metadata needed to derive the canonical ID are
REPORTED and left untouched (never guessed). Key attribution
(``created_by_key_hash``) is always preserved.

By default the script is a DRY RUN: it prints the plan and writes
nothing. ``--apply`` performs the rewrites; the old files are kept by
default so old references keep resolving during transition, and can be
pruned afterward with ``--remove-legacy`` (only together with
``--apply``).

Usage:
    python scripts/migrate_audit_ids_od07.py            # dry run
    python scripts/migrate_audit_ids_od07.py --apply    # migrate, keep old files
    python scripts/migrate_audit_ids_od07.py --apply --remove-legacy
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))

from tos.governance.audit_identity import (  # noqa: E402
    audit_id_from_record, is_full_audit_id,
)
from tos.storage.audit_store import (  # noqa: E402
    _get_audit_dir, compute_legacy_audit_id,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true",
                    help="perform the rewrites (default: dry run)")
    ap.add_argument("--remove-legacy", action="store_true",
                    help="with --apply: delete the old files after re-keying")
    args = ap.parse_args()
    if args.remove_legacy and not args.apply:
        ap.error("--remove-legacy requires --apply")

    audit_dir = _get_audit_dir()
    migrated, skipped, unchanged = [], [], []
    for path in sorted(audit_dir.glob("*.json")):
        key = path.stem
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as e:
            skipped.append((key, f"unreadable: {e}"))
            continue
        canonical = audit_id_from_record(record)
        if canonical is None:
            skipped.append((key, "insufficient metadata (no provenance "
                                 "hash + fingerprint); left untouched"))
            continue
        if key == canonical and \
                record.get("governance", {}).get("audit_id") == canonical:
            unchanged.append(key)
            continue
        gov = record.setdefault("governance", {})
        gov.setdefault("legacy_audit_id", None)
        if not gov["legacy_audit_id"]:
            prov_hash = record.get("provenance", {}).get("hash")
            if prov_hash:
                gov["legacy_audit_id"] = compute_legacy_audit_id(prov_hash)
        migrated.append((key, canonical))
        if args.apply:
            gov["audit_id"] = canonical
            new_path = path.parent / f"{canonical}.json"
            new_path.write_text(json.dumps(
                record, indent=2, sort_keys=True, default=str))
            if args.remove_legacy:
                path.unlink()

    print(f"audit dir: {audit_dir}")
    print(f"already canonical (untouched): {len(unchanged)}")
    print(f"to migrate / migrated: {len(migrated)}")
    for old, new in migrated:
        print(f"  {old} -> {new}")
    print(f"skipped (reported, untouched): {len(skipped)}")
    for key, why in skipped:
        print(f"  {key}: {why}")
    print(f"mode: {'APPLIED' if args.apply else 'DRY RUN (nothing written)'}")
    if skipped and args.apply:
        print("note: skipped records were NOT modified; migrate them "
              "manually after supplying their metadata.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
