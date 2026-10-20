import sys
from entries._094_mdi.check_candidate import violates_mdi, all_violations

def parse_args(argv):
    args = {"target": argv[1], "generator": argv[2]}
    if "--scaffold" in argv:
        idx = argv.index("--scaffold")
        args["scaffold"] = argv[idx + 1]
    if "--audit" in argv:
        args["audit"] = True
    return args

def plan_or_abort(candidate):
    if candidate.get("audit"):
        violations = all_violations(candidate)
        if violations:
            print("⛔ AUDIT REPORT — MULTIPLE VIOLATIONS")
            for law in violations:
                print(f" - {law['doctrine_id']} ({law.get('title')})")
            sys.exit(1)
        else:
            print("✅ AUDIT PASS — no violations")
            return True

    hit = violates_mdi(candidate)
    if hit:
        print("⛔ PLANNER ABORT")
        print(f"   Violated Doctrine: {hit['doctrine_id']}")
        print(f"   Verdict: {hit['verdict']}")
        sys.exit(1)

    print("✅ PLANNER PASS")
    print("   No locked doctrine violated")
    return True

if __name__ == "__main__":
    candidate = parse_args(sys.argv)
    plan_or_abort(candidate)
