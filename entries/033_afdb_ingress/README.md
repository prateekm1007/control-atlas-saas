# Entry 033 — AlphaFold DB Ingress

Proteome-scale ingestion engine.

1. **Fetches** structures from AlphaFold DB.
2. **Filters** by global pLDDT (default > 70).
3. **Detects** pockets using Entry 027 physics.
4. **Indexes** VALIDATED pockets to `library/atlas_index`.

## Usage
```bash
python ingest_proteome.py --ids my_targets.txt
