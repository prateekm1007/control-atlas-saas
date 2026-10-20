#!/usr/bin/env python3
"""Runner for the acceptance tests AT-01..AT-34 (EC-20 dogfood entry point).

The tests live in test_acceptance.py as pytest-compatible node IDs
(scripts/ci/tests/test_acceptance.py::test_atNN_*), so `pytest` can execute and
collect them (that is what the EC-20 mechanism executes). This wrapper runs the
same functions directly when pytest is unavailable and prints explicit counts.

Exit 0 = all pass; 1 = any failure.
"""
import importlib.util
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "test_acceptance.py")


def main() -> int:
    if (r := subprocess.run([sys.executable, "-m", "pytest", SRC, "-q", "--no-header"],
                            capture_output=True, text=True)).returncode == 0:
        print(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else "pytest: ok")
        print("ALL ACCEPTANCE TESTS PASS (AT-01..AT-34, pytest runner)")
        return 0
    if "No module named pytest" not in r.stderr:
        print(r.stdout)
        print(r.stderr)
        print(f"FAILURES via pytest runner (exit {r.returncode})")
        return 1
    # pytest unavailable: execute the test functions directly
    spec = importlib.util.spec_from_file_location("test_acceptance", SRC)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    fns = [getattr(mod, n) for n in sorted(dir(mod)) if re.fullmatch(r"test_at\d+_\w+", n)]
    for fn in fns:
        try:
            fn()
        except AssertionError as e:
            mod.FAIL.append(f"{fn.__name__} raised: {e}")
        except Exception as e:  # noqa: BLE001
            mod.FAIL.append(f"{fn.__name__} ERROR: {type(e).__name__}: {e}")
    print(f"PASS {len(mod.PASS)}: " + "; ".join(mod.PASS))
    if mod.FAIL:
        print(f"FAIL {len(mod.FAIL)}: " + "; ".join(mod.FAIL))
        return 1
    print(f"ALL ACCEPTANCE TESTS PASS ({len(fns)} groups, direct runner)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
