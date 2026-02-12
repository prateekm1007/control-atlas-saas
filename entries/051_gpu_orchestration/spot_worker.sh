#!/bin/bash
# Worker script for Spot Instances (Azure/AWS)
# Usage: ./spot_worker.sh <batch_file.json>

BATCH_FILE=$1
RESULTS_DIR="rfaa_results"

mkdir -p $RESULTS_DIR

echo "[*] Processing Batch: $BATCH_FILE"

# Parse JSON loop (simplified for bash/python hybrid env)
# In prod, this would be a python script importing rfaa_wrapper
# Here we demonstrate the invocation

python3 -c "
import json
import sys
import os
# Add library path
sys.path.append('../../library/rfaa')
from rfaa_wrapper import RFAAEngine
from pathlib import Path

with open('$BATCH_FILE') as f:
    ids = json.load(f)

engine = RFAAEngine(Path('$RESULTS_DIR'), use_docker=True)

for hid in ids:
    target, compound = hid.split('_', 1)
    # Mock lookup of FASTA/SDF (In prod: fetch from Atlas/Library)
    fasta = f'{target}.fasta'
    sdf = f'{compound}.sdf'
    
    # Run RFAA
    print(f' -> Validating {hid}...')
    try:
        # Checkpoint check
        if (Path('$RESULTS_DIR') / f'{target}_complex/result.json').exists():
            print('    Already done.')
            continue
            
        # Creating dummy inputs for demo if missing
        if not os.path.exists(fasta): Path(fasta).touch()
        if not os.path.exists(sdf).touch()
        
        # Run
        res = engine.validate_complex(hid, fasta, sdf)
        
        # Save result for KG ingestion later
        with open(Path('$RESULTS_DIR') / f'{target}_complex/result.json', 'w') as f:
            json.dump(res, f)
            
    except Exception as e:
        print(f'    Failed: {e}')
"

echo "[+] Batch Complete. Sync results to S3/Blob now."
