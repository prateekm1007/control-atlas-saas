import sys, os, json

def calculate_3x3_matrix(results):
    states = ["PASS", "VETO", "INDET"]
    matrix = {s: {s2: 0 for s2 in states} for s in states}
    
    for r in results:
        matrix[r['actual']][r['exp']] += 1
    
    print("\n3x3 CONCORDANCE MATRIX (Actual vs Expected):")
    print(f"{'':<8} | {'EXP PASS':<10} | {'EXP VETO':<10} | {'EXP INDET':<10}")
    for actual in states:
        print(f"ACT {actual:<4} | {matrix[actual]['PASS']:<10} | {matrix[actual]['VETO']:<10} | {matrix[actual]['INDET']:<10}")

if __name__ == "__main__":
    # Pilot targets (50 high-res PDB study simulated here)
    mock_results = [
        {"id": "1A2B", "actual": "PASS", "exp": "PASS"},
        {"id": "BAD1", "actual": "VETO", "exp": "VETO"},
        {"id": "AF_D", "actual": "INDET", "exp": "INDET"},
    ]
    calculate_3x3_matrix(mock_results)
