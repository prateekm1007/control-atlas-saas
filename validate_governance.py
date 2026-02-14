import sys, os, subprocess, json

TARGETS = [
    ("AF-P01308-F1", "Insulin (Standard)"),
    ("AF-P06213-F1", "Insulin (AFDB)"),
    ("AF-O43526-F1", "Disordered (NRIP1)"),
    ("AF-P68871-F1", "Hemoglobin (Beta)")
]

print("\n" + "═"*80)
print("  TOSCANINI v21.48.0 GOVERNANCE VALIDATION REPORT")
print("  Calibration Mode: CONSERVATIVE (12 Deterministic Laws)")
print("═"*80)

results_table = []
for sid, name in TARGETS:
    inner_cmd = (
        "import sys, json; "
        "sys.path.insert(0, '/app'); "
        "from tos.generation.dispatcher import GenerationDispatcher; "
        "from tos.ingestion.processor import IngestionProcessor; "
        "from tos.engine.tier1_measurements import Tier1Measurements; "
        "from tos.utils.type_guards import force_bytes; "
        f"res = GenerationDispatcher.acquire('{sid}', None); "
        "struct = IngestionProcessor.run(force_bytes(res[0]), 't.pdb', 'Audit', 'Discovery'); "
        "audit = Tier1Measurements.run_full_audit(struct); "
        "det_pass = sum(1 for lid, r in audit.items() if r[0]=='PASS' and r[3]=='deterministic'); "
        f"print(json.dumps({{'name': '{name}', 'pass': det_pass, 'verdict': 'PASS' if det_pass==12 else 'VETO', "
        f"'coverage': audit['LAW-105'][1], 'clash': audit['LAW-130'][1]}}))"
    )
    cmd = f"docker exec brain python3 -c \"{inner_cmd}\""
    try:
        out = subprocess.check_output(cmd, shell=True).decode()
        results_table.append(json.loads(out))
    except: results_table.append({'name': name, 'verdict': 'CRASH'})

print("\n" + "═"*115)
print(f"{'Target':<20} | {'Coverage (LAW-105)':<35} | {'Clash (LAW-130)':<20} | {'Pass':<6} | {'Verdict'}")
print("-" * 115)
for r in results_table:
    print(f"{r['name']:<20} | {r.get('coverage','N/A'):<35} | {r.get('clash','N/A'):<20} | {r.get('pass','N/A'):<6} | {r['verdict']}")
print("═"*115 + "\n")
