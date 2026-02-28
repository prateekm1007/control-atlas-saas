#!/usr/bin/env python3
"""
scripts/provision_pilot_key.py
Operator tool — provision one pilot API key manually.

Usage:
    TOSCANINI_DATA_DIR=/app/data python3 scripts/provision_pilot_key.py \
        --tier pro --note "researcher@university.edu"
"""
import argparse, json, os, sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier", choices=["free","pro","enterprise"], default="pro")
    parser.add_argument("--note", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    os.environ.setdefault("TOSCANINI_DATA_DIR", "/app/data")
    from tos.saas.key_store import create_key, hash_key, GPU_TIER_CREDITS, TIER_LIMITS

    if args.dry_run:
        print(f"[DRY RUN] Would create {args.tier} key for: {args.note}")
        return

    key      = create_key(args.tier)
    key_hash = hash_key(key)
    limits   = TIER_LIMITS[args.tier]
    gpu      = GPU_TIER_CREDITS[args.tier]

    record = {
        "provisioned_at": datetime.now(timezone.utc).isoformat(),
        "tier":           args.tier,
        "note":           args.note,
        "key_hash":       key_hash,
        "limits":         {"daily_audits": limits["daily"], "gpu_runs": gpu},
    }
    roster = Path(__file__).parent / "pilot_roster.jsonl"
    with open(roster, "a") as f:
        f.write(json.dumps(record) + "\n")

    print()
    print("=" * 60)
    print("  PILOT KEY — COPY NOW, NOT SHOWN AGAIN")
    print("=" * 60)
    print(f"  Key  : {key}")
    print(f"  Tier : {args.tier}")
    print(f"  Note : {args.note}")
    print(f"  GPU  : {gpu} runs")
    print("=" * 60)

if __name__ == "__main__":
    main()
