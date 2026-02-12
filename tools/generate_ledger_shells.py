from pathlib import Path

base = Path("ledger/motifs")
motifs_file = base / "WEEK1_MOTIFS.txt"

if not motifs_file.exists():
    raise SystemExit(f"❌ Missing {motifs_file}")

motifs = motifs_file.read_text().splitlines()

for m in motifs:
    m = m.strip()
    if not m or m.startswith("#"):
        continue

    target = base / "pending" / m
    target.mkdir(parents=True, exist_ok=True)

    layers = {
        "INTENT.md": "Design strategy and sequence hypothesis.",
        "OUTCOME.md": "Raw results from Chai-1 or explicit failure.",
        "JUDGMENT.md": "Clearance, energy, MD stability.",
        "DOCTRINE.md": "Linked laws (LAW-XXX) and NKG references."
    }

    for name, desc in layers.items():
        f = target / name
        if not f.exists():
            f.write_text(
                f"# {name.replace('.md','')} — {m}\n\n"
                f"{desc}\n\n"
                "- **Status:** PENDING\n"
            )

print(f"✅ Generated governed shells for {len([m for m in motifs if m and not m.startswith('#')])} motifs.")
