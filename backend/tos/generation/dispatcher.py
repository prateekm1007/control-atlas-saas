import requests
import logging

logger = logging.getLogger("toscanini.dispatcher")

class GenerationDispatcher:
    """PILLAR 05: Industrial Acquisition - Multi-format AFDB Downloader."""
    
    @staticmethod
    def acquire(candidate_id: str, fallback):
        # Extract base ID (e.g., AF-P01308-F1 -> P01308)
        clean_id = candidate_id.replace("AF-", "").replace("-F1", "").strip()
        
        # AFDB v4 Standard URL Patterns
        urls = [
            f"https://alphafold.ebi.ac.uk/files/AF-{clean_id}-F1-model_v4.pdb",
            f"https://alphafold.ebi.ac.uk/files/AF-{clean_id}-F1-model_v4.cif",
            f"https://alphafold.ebi.ac.uk/api/prediction/{clean_id}" # API fallback
        ]
        
        for url in urls:
            try:
                logger.info(f"Acquisition Attempt: {url}")
                res = requests.get(url, timeout=15)
                if res.status_code == 200:
                    content = res.content
                    # If it's the API JSON, extract the actual file URL
                    if "application/json" in res.headers.get("Content-Type", ""):
                        data = res.json()
                        if isinstance(data, list) and len(data) > 0:
                            actual_url = data[0].get("pdbUrl")
                            if actual_url:
                                content = requests.get(actual_url).content
                    
                    if len(content) > 200:
                        fmt = "cif" if ".cif" in url else "pdb"
                        logger.info(f"✅ Success: Acquired {len(content)} bytes ({fmt})")
                        return (content, candidate_id, fmt)
            except Exception as e:
                logger.warning(f"Acquisition segment failed: {e}")
                continue

        logger.warning(f"❌ All remote sources exhausted for {candidate_id}. Deploying Institutional Mock.")
        mock = b"ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 50.00           N\nATOM      2  CA  MET A   1       1.458   0.000   0.000  1.00 50.00           C\nEND\n"
        return (mock, candidate_id, "pdb")
