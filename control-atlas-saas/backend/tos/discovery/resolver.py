import requests, re

class DiscoveryResolver:
    UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search?query={}&format=json&size=15"

    @staticmethod
    def resolve(query: str):
        query = query.strip()
        candidates = []

        # 1. ACCESSION MATCH (Direct to AFDB)
        acc_match = re.search(r"([A-Z0-9]{6,10})", query.upper())
        if acc_match:
            uid = acc_match.group(1)
            candidates.append({
                "id": f"AFDB:{uid}:F1", "label": f"📦 AlphaFold DB Primary ({uid})",
                "source": "AFDB", "type": "Evidence", "hint": "Canonical structural evidence"
            })

        # 2. SEMANTIC SEARCH (UniProt Matches)
        try:
            resp = requests.get(DiscoveryResolver.UNIPROT_SEARCH_URL.format(query), timeout=5)
            if resp.status_code == 200:
                for entry in resp.json().get('results', []):
                    uid = entry['primaryAccession']
                    name = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
                    if any(c['id'].split(":")[1] == uid for c in candidates): continue
                    candidates.append({
                        "id": f"AFDB:{uid}:F1", "label": f"📦 {name} ({uid})",
                        "source": "AFDB", "type": "Evidence", "hint": "High-fidelity evidence"
                    })
        except: pass

        # 3. HOSTED INFERENCE (ESM Atlas - Working Meta API)
        candidates.append({
            "id": "ESMATLAS:OnDemand", "label": "🔬 ESM Atlas Fast Fold (Hosted Meta API)",
            "source": "ESMAtlas", "type": "Inference", "hint": "Instant generation via Meta cloud"
        })

        return candidates
