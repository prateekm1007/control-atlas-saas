import requests

class GenerationDispatcher:
    AFDB_VERSIONS = ["v4", "v6", "v5"] # v4 is often more stable for canonicals like Insulin
    HEADERS = {"User-Agent": "Mozilla/5.0", "Accept": "*/*"}

    @staticmethod
    def acquire(candidate_id: str, sequence: str = None):
        # --- PATH 1: AlphaFold DB (Evidence) ---
        if candidate_id.startswith("AFDB"):
            parts = candidate_id.split(":")
            uid, isoform = parts[1], parts[2][1:]
            
            for ver in GenerationDispatcher.AFDB_VERSIONS:
                for ext in ["pdb", "cif"]:
                    url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F{isoform}-model_{ver}.{ext}"
                    try:
                        res = requests.get(url, headers=GenerationDispatcher.HEADERS, timeout=8)
                        if res.status_code == 200: return res.content, f"AlphaFold DB {ver}", ext
                    except: continue
            return None, "Evidence not found in AFDB (Check Accession ID)", "pdb"

        # --- PATH 2: ESM Atlas (Hosted Inference) ---
        elif candidate_id.startswith("ESMATLAS"):
            if not sequence: return None, "Please provide a sequence for the Meta API", "pdb"
            try:
                res = requests.post("https://api.esmatlas.com/foldSequence/v1/pdb/", 
                                   data=sequence.strip(), headers={"Content-Type": "text/plain"}, timeout=30)
                if res.status_code == 200: return res.content, "ESM Atlas (Meta Cloud)", "pdb"
                return None, f"Meta API Error: {res.status_code}", "pdb"
            except: return None, "Meta Cloud unreachable (Timeout)", "pdb"

        return None, "Invalid Protocol", "pdb"
