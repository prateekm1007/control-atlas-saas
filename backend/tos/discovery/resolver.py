import requests
import logging

logger = logging.getLogger("toscanini.discovery")

class DiscoveryResolver:
    @staticmethod
    def resolve(query: str):
        query_upper = query.strip().upper()
        results = []
        
        # Strategy: UniProt keyword search EXHAUSTIVE (size=50)
        try:
            url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&size=50&format=json"
            res = requests.get(url, timeout=15)
            if res.status_code == 200:
                for entry in res.json().get("results", []):
                    acc = entry.get("primaryAccession", "")
                    name = entry.get("proteinDescription", {}).get(
                        "recommendedName", {}).get("fullName", {}).get("value", query_upper)
                    org = entry.get("organism", {}).get("scientificName", "Unknown")
                    if acc:
                        results.append({
                            "id": f"AF-{acc}-F1",
                            "label": f"AlphaFold: {name} ({acc}) - {org}",
                            "source": "UniProt",
                            "pdb_url": f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb",
                            "cif_url": f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.cif"
                        })
                if results: return results
        except Exception as e:
            logger.warning(f"Exhaustive search failed: {e}")
        
        return [{"id": f"AF-{query_upper}-F1", "label": f"Fallback: {query_upper}", "source": "FALLBACK"}]
