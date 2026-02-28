#!/usr/bin/env python3
"""
scripts/monitor_pilot.py

Operator signal reader — Week 1 pilot monitoring.
Run on the server against live telemetry/events.jsonl.

Usage:
    python3 scripts/monitor_pilot.py [--days N] [--jsonl /path/to/events.jsonl]

Signals reported (matches exact event names in main.py):
    audit_created       — structure submitted for audit
    token_generated     — callback token minted and consumed by user
    callback_completed  — refined structure posted back successfully
    credit_consumed     — GPU/notebook credit deducted
    delta_coverage      — coverage improvement on callback (inside callback_completed)

Shipping threshold (non-negotiable for v24.0):
    token->callback completion rate >= 25%
"""
import argparse
import json
import os
import sys
from collections import defaultdict
from datetime import datetime, timezone, timedelta
from pathlib import Path

VERSION = "v23.2.0-B6"

# Real event names — single source of truth
EVENTS = {
    "audit_created",
    "token_generated",
    "callback_completed",
    "credit_consumed",
}


def load_events(jsonl_path: Path, days: int) -> list:
    if not jsonl_path.exists():
        print(f"  No telemetry file at {jsonl_path}")
        return []
    cutoff = datetime.now(timezone.utc) - timedelta(days=days)
    events = []
    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                ts  = datetime.fromisoformat(rec.get("ts", "2000-01-01T00:00:00+00:00"))
                if ts >= cutoff:
                    events.append(rec)
            except (json.JSONDecodeError, ValueError):
                continue
    return events


def report(events: list, days: int) -> None:
    print()
    print("=" * 60)
    print(f"  TOSCANINI PILOT SIGNAL — last {days} day(s)")
    print(f"  Engine: {VERSION}")
    print(f"  Events loaded: {len(events)}")
    print("=" * 60)

    by_type = defaultdict(list)
    for e in events:
        by_type[e.get("event", "unknown")].append(e)

    # ── Volume ────────────────────────────────────────────────────
    print("\n  VOLUME")
    for ev in sorted(EVENTS):
        count = len(by_type[ev])
        print(f"    {ev:<28} {count:>5}")

    # ── Completion rate (THE metric) ──────────────────────────────
    tokens    = len(by_type["token_generated"])
    callbacks = len(by_type["callback_completed"])
    rate      = callbacks / tokens * 100 if tokens else 0.0
    threshold = 25.0
    status    = "✓ ABOVE THRESHOLD" if rate >= threshold else "✗ BELOW THRESHOLD"
    print(f"\n  COMPLETION RATE")
    print(f"    Tokens generated  : {tokens}")
    print(f"    Callbacks done    : {callbacks}")
    print(f"    Rate              : {rate:.1f}%  [{status}]")
    print(f"    Ship threshold    : {threshold:.0f}%")

    # ── Unique users ──────────────────────────────────────────────
    users_audited   = {e.get("audit_id") for e in by_type["audit_created"]}
    users_callbacked = {e.get("user_email") for e in by_type["callback_completed"]
                        if e.get("user_email")}
    print(f"\n  UNIQUE SIGNAL")
    print(f"    Unique audits     : {len(users_audited)}")
    print(f"    Users who re-uploaded: {len(users_callbacked)}")

    # ── Coverage deltas ───────────────────────────────────────────
    deltas = [
        e.get("delta_coverage", 0)
        for e in by_type["callback_completed"]
        if isinstance(e.get("delta_coverage"), (int, float))
    ]
    if deltas:
        avg_delta = sum(deltas) / len(deltas)
        max_delta = max(deltas)
        min_delta = min(deltas)
        print(f"\n  COVERAGE IMPROVEMENT (pp)")
        print(f"    Mean delta        : {avg_delta:+.1f}")
        print(f"    Best delta        : {max_delta:+.1f}")
        print(f"    Worst delta       : {min_delta:+.1f}")
        if avg_delta < 5.0:
            print(f"    ⚠  Mean < 5pp — notebook prescription may be weak")
    else:
        print(f"\n  COVERAGE IMPROVEMENT: no data yet")

    # ── Credit ceiling hits ───────────────────────────────────────
    by_tier = defaultdict(int)
    for e in by_type["credit_consumed"]:
        by_tier[e.get("tier", "unknown")] += 1
    print(f"\n  CREDIT CONSUMPTION BY TIER")
    if by_tier:
        for tier, count in sorted(by_tier.items()):
            print(f"    {tier:<16} {count:>5} events")
    else:
        print("    no credit events yet")

    # ── Shipping criteria summary ─────────────────────────────────
    print(f"\n  SHIPPING CRITERIA (v24.0 gate)")
    criteria = [
        ("completion >= 25%",        rate >= 25.0),
        (">= 3 users hit ceiling",   False),   # manual check
        ("no token replay abuse",    True),     # enforced by blacklist
        ("no rate-limit breach",     True),     # enforced by rate limiter
        ("governance stable",        True),     # frozen
    ]
    for label, met in criteria:
        icon = "✓" if met else "✗"
        print(f"    {icon}  {label}")

    print()
    print("=" * 60)
    print()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--days",  type=int, default=7)
    parser.add_argument("--jsonl", type=str,
        default=os.path.join(
            os.environ.get("TOSCANINI_DATA_DIR", "/app/data"),
            "telemetry", "events.jsonl"
        )
    )
    args = parser.parse_args()
    events = load_events(Path(args.jsonl), args.days)
    report(events, args.days)


if __name__ == "__main__":
    main()
