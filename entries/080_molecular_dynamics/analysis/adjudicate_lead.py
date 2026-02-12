import os
import pandas as pd
import numpy as np

def adjudicate_energy(log_file):
    if not os.path.exists(log_file):
        return "WAITING", 0.0
        
    df = pd.read_csv(log_file)
    # We ignore the first 10% of the run as equilibration noise
    pe = df['Potential Energy (kJ/mole)'].values
    offset = len(pe) // 10
    stable_pe = pe[offset:]
    
    # Simple linear drift check (slope should oscillate around 0)
    # Positive slope > 1.0 indicates steady energy leakage (Fail)
    slope = np.polyfit(range(len(stable_pe)), stable_pe, 1)[0]

    if slope > 1.0:
        return "FAIL_ENERGY_LEAK", slope
    return "PASS", slope

if __name__ == "__main__":
    log_file = "../outputs/md_log.csv"
    status, metric = adjudicate_energy(log_file)
    print(f"⚖️ ADJUDICATION RESULT: {status}")
    print(f"📊 Energy slope: {metric:.5f}")
