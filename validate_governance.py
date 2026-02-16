import sys, os, subprocess, json

TARGETS = [
    ("4HHB", "Hemoglobin (Crystal)"),
    ("1CRN", "Crambin (Crystal)"),
    ("AF-P01308-F1", "Insulin (AlphaFold)"),
    ("AF-P68871-F1", "Hemoglobin (AlphaFold)"),
]

print("\n" + "=" * 80)
print("  TOSCANINI v22.5.2 GOVERNANCE VALIDATION REPORT")
print("  Calibration Mode: CONSERVATIVE (12 Deterministic Laws)")
print("  PIL-CAL-02: LAW-100 advisory for experimental structures")
print("=" * 80)

results_table = []
for sid, name in TARGETS:
    # Use the /ingest API directly — this is the canonical scoring path
    cmd = (
        f'curl -s -X POST "http://localhost:8000/ingest" '
        f'-F "mode=Benchmark" '
        f'-F "candidate_id={sid}" '
        f'-F "t3_category=NONE"'
    )
    try:
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=180)
        if proc.returncode != 0:
            results_table.append({"name": name, "verdict": "CRASH", "error": proc.stderr.strip()[-200:]})
            continue
        d = json.loads(proc.stdout)
        v = d.get("verdict", {})
        c = d.get("characterization", {})
        laws = d.get("tier1", {}).get("laws", [])

        # Count deterministic laws honestly — exclude advisory_experimental
        det_laws = [l for l in laws if l["method"] == "deterministic"]
        det_passed = sum(1 for l in det_laws if l["status"] == "PASS")
        det_total = len(det_laws)

        # Count advisory_experimental separately
        adv_exp = [l for l in laws if l["method"] == "advisory_experimental"]
        adv_exp_info = ""
        if adv_exp:
            adv_exp_info = f" +{len(adv_exp)} advisory"

        results_table.append({
            "name": name,
            "source": "EXP" if c.get("source_type") == "experimental" else "AF",
            "resolution": c.get("resolution", "N/A"),
            "coverage": v.get("coverage_pct", "N/A"),
            "clash": next((l["observed"] for l in laws if l["law_id"] == "LAW-130"), "N/A"),
            "det_pass": det_passed,
            "det_total": det_total,
            "adv_note": adv_exp_info,
            "verdict": v.get("binary", "ERROR"),
            "law100_obs": next((l["observed"] for l in laws if l["law_id"] == "LAW-100"), "N/A"),
            "law100_method": next((l["method"] for l in laws if l["law_id"] == "LAW-100"), "?"),
        })
    except Exception as e:
        results_table.append({"name": name, "verdict": "CRASH", "error": str(e)[:200]})

print("\n" + "=" * 145)
print(f"{'Target':<24} | {'Src':<4} | {'Res':<5} | {'Cov':<7} | {'Clash':<8} | {'Det Pass':<12} | {'LAW-100':<20} | {'Verdict'}")
print("-" * 145)
for r in results_table:
    if r["verdict"] == "CRASH":
        print(f"{r['name']:<24} | {'?':<4} | {'?':<5} | {'?':<7} | {'?':<8} | {'?':<12} | {'CRASH':<20} | {r.get('error','')[:40]}")
        continue
    det_str = f"{r['det_pass']}/{r['det_total']}{r.get('adv_note','')}"
    law100_str = f"{r['law100_obs']} ({r['law100_method'][:8]})"
    res_str = str(r.get('resolution', 'N/A'))
    print(
        f"{r['name']:<24} | {r['source']:<4} | {res_str:<5} | "
        f"{str(r['coverage']):<7} | {str(r['clash']):<8} | "
        f"{det_str:<12} | {law100_str:<20} | {r['verdict']}"
    )
print("=" * 145)

# Summary
crystals = [r for r in results_table if r.get("source") == "EXP" and r["verdict"] != "CRASH"]
af = [r for r in results_table if r.get("source") == "AF" and r["verdict"] != "CRASH"]
crystal_pass = sum(1 for r in crystals if r["verdict"] == "PASS")
af_pass = sum(1 for r in af if r["verdict"] == "PASS")
false_veto = sum(1 for r in crystals if r["verdict"] != "PASS")

print(f"\nCrystal structures: {crystal_pass}/{len(crystals)} PASS (false veto rate: {false_veto}/{len(crystals)})")
print(f"AlphaFold models:   {af_pass}/{len(af)} PASS")
if false_veto == 0 and crystals:
    print(f"\n✅ No false veto events against {len(crystals)} experimental crystal structures.")
print()
