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
print("=" * 80)

results_table = []
for sid, name in TARGETS:
    inner_cmd = (
        "import sys, json; "
        "sys.path.insert(0, '/app'); "
        "from tos.generation.dispatcher import GenerationDispatcher; "
        "from tos.ingestion.processor import IngestionProcessor; "
        "from tos.engine.tier1_measurements import Tier1Measurements; "
        "from tos.governance.station_sop import LAW_METHOD_CLASSIFICATIONS; "
        "from tos.utils.type_guards import force_bytes; "
        f"res = GenerationDispatcher.acquire('{sid}', None); "
        "struct = IngestionProcessor.run(force_bytes(res[0]), 't.pdb', 'Audit', 'Discovery'); "
        "audit, cov, ff = Tier1Measurements.run_full_audit(struct); "
        "det_pass = sum(1 for lid, r in audit.items() "
        "if r['status']=='PASS' and LAW_METHOD_CLASSIFICATIONS.get(lid)=='deterministic'); "
        "src = 'EXP' if struct.confidence.is_experimental else 'AF'; "
        f"print(json.dumps({{'name': '{name}', 'pass': det_pass, "
        f"'verdict': 'PASS' if det_pass==12 else 'VETO', "
        f"'coverage': audit.get('LAW-105', {{}}).get('observed', 'N/A'), "
        f"'clash': audit.get('LAW-130', {{}}).get('observed', 'N/A'), "
        f"'source': src}}))"
    )
    cmd = f'docker exec brain python3 -c "{inner_cmd}"'
    try:
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        if proc.returncode != 0:
            results_table.append({"name": name, "verdict": "CRASH", "error": proc.stderr.strip()[-200:]})
            continue
        found = False
        for line in reversed(proc.stdout.strip().split("\n")):
            line = line.strip()
            if line.startswith("{"):
                results_table.append(json.loads(line))
                found = True
                break
        if not found:
            results_table.append({"name": name, "verdict": "CRASH", "error": "No JSON in stdout"})
    except Exception as e:
        results_table.append({"name": name, "verdict": "CRASH", "error": str(e)[:200]})

print("\n" + "=" * 130)
print(f"{'Target':<24} | {'Source':<6} | {'Coverage':<12} | {'Clash':<12} | {'Pass':<6} | {'Verdict'}")
print("-" * 130)
for r in results_table:
    v = r["verdict"]
    if v == "CRASH":
        v = f"CRASH: {r.get('error', '')[:60]}"
    print(
        f"{r['name']:<24} | {str(r.get('source', '?')):<6} | "
        f"{str(r.get('coverage', 'N/A')):<12} | {str(r.get('clash', 'N/A')):<12} | "
        f"{str(r.get('pass', 'N/A')):<6} | {v}"
    )
print("=" * 130 + "\n")
