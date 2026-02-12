import sys
from pathlib import Path

status = Path("GENERATOR_STATUS.md").read_text()

if "UNBOUND" in status:
    print("⛔ Atom minting disabled — generator backend is UNBOUND")
    sys.exit(1)

print("✅ Generator backend bound")
