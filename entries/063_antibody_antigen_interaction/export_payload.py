#!/usr/bin/env python3
import json

MANIFEST = "interaction_jobs_063.json"

def main():
    with open(MANIFEST) as f:
        jobs = json.load(f)

    print("\n" + "="*20 + " COPY BELOW THIS LINE " + "="*20)
    print("JOBS_PAYLOAD = [")
    for job in jobs:
        entry = {
            "id": job["job_id"],
            "antibody_seq": job["antibody"]["sequence"],
            "antigen_seq": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQ"
        }
        print(f"    {json.dumps(entry)},")
    print("]")
    print("="*20 + " COPY ABOVE THIS LINE " + "="*20 + "\n")

if __name__ == "__main__":
    main()
