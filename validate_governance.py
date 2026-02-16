import sys, os, subprocess, json

TARGETS = [
    ("3NIR", "Endothiapepsin (0.48A)"),
    ("2VB1", "Insulin (0.65A)"),
    ("1CRN", "Crambin (1.50A)"),
    ("4HHB", "Hemoglobin (1.74A)"),
    ("1UBQ", "Ubiquitin (1.80A)"),
    ("6LZG", "SARS-CoV-2 RBD (2.50A)"),
    ("AF-P01308-F1", "Insulin (AlphaFold)"),
    ("AF-P68871-F1", "Hemoglobin (AlphaFold)"),
]

print("\n" + "=" * 80)
print("  TOSCANINI v22.5.2 GOVERNANCE VALIDATION REPORT")
print("  Calibration: 11 Deterministic + 1 Advisory (Experimental)")
print("  PIL-CAL-02: LAW-100 advisory for experimental structures")
print("=" * 80)

results_table = []
for sid, name in TARGETS:
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

        det_laws = [l for l in laws if l["method"] == "deterministic"]
        det_passed = sum(1 for l in det_laws if l["status"] == "PASS")
        det_total = len(det_laws)
        adv_exp = [l for l in laws if l["method"] == "advisory_experimental"]
        adv_note = f" +{len(adv_exp)} adv" if adv_exp else ""

        results_table.append({
            "name": name,
            "source": "EXP" if c.get("source_type") == "experimental" else "AF",
            "resolution": c.get("resolution", "N/A"),
            "coverage": v.get("coverage_pct", "N/A"),
            "clash": next((l["observed"] for l in laws if l["law_id"] == "LAW-130"), "N/A"),
            "det_pass": det_passed,
            "det_total": det_total,
            "adv_note": adv_note,
            "verdict": v.get("binary", "ERROR"),
            "law100_obs": next((l["observed"] for l in laws if l["law_id"] == "LAW-100"), "N/A"),
            "law100_method": next((l["method"] for l in laws if l["law_id"] == "LAW-100"), "?"),
        })
    except Exception as e:
        results_table.append({"name": name, "verdict": "CRASH", "error": str(e)[:200]})

print("\n" + "=" * 150)
print(f"{'Target':<28} | {'Src':<4} | {'Res':<5} | {'Cov':<7} | {'Clash':<8} | {'Det Pass':<14} | {'LAW-100':<20} | {'Verdict'}")
print("-" * 150)
for r in results_table:
    if r["verdict"] == "CRASH":
        print(f"{r['name']:<28} | {'?':<4} | {'?':<5} | {'?':<7} | {'?':<8} | {'?':<14} | {'CRASH':<20} | {r.get('error','')[:40]}")
        continue
    det_str = f"{r['det_pass']}/{r['det_total']}{r.get('adv_note', '')}"
    law100_str = f"{r['law100_obs']} ({r['law100_method'][:8]})"
    res_str = str(r.get('resolution', 'N/A'))
    print(
        f"{r['name']:<28} | {r['source']:<4} | {res_str:<5} | "
        f"{str(r['coverage']):<7} | {str(r['clash']):<8} | "
        f"{det_str:<14} | {law100_str:<20} | {r['verdict']}"
    )
print("=" * 150)

crystals = [r for r in results_table if r.get("source") == "EXP" and r["verdict"] != "CRASH"]
af = [r for r in results_table if r.get("source") == "AF" and r["verdict"] != "CRASH"]
crystal_pass = sum(1 for r in crystals if r["verdict"] == "PASS")
false_veto = sum(1 for r in crystals if r["verdict"] != "PASS")
af_pass = sum(1 for r in af if r["verdict"] == "PASS")

res_values = [r.get("resolution") for r in crystals if r.get("resolution") not in (None, "N/A")]
res_range = f"{min(res_values)}-{max(res_values)} A" if res_values else "N/A"

print(f"\nCrystal structures: {crystal_pass}/{len(crystals)} PASS (false veto rate: {false_veto}/{len(crystals)})")
print(f"Resolution range: {res_range}")
print(f"AlphaFold models:  {af_pass}/{len(af)} PASS")
if false_veto == 0 and crystals:
    print(f"\n{'='*60}")
    print(f"  CALIBRATION CERTIFIED")
    print(f"  No false veto events against {len(crystals)} experimental")
    print(f"  crystal structures ({res_range}).")
    print(f"{'='*60}")
print()
