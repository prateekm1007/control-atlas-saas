# Entry 031 — DUD-E Benchmark

Validates Control Atlas against the DUD-E benchmark dataset.

## Usage

1. Download data: python download_dude.py
2. Run benchmark: python run_benchmark.py
3. View results: cat results/benchmark_results.json

## Metrics

- TPR: True Positive Rate (actives correctly accepted)
- TNR: True Negative Rate (decoys correctly rejected)
- Enrichment Factor: How much better than random

## Targets

EGFR, BRAF, CDK2, SRC, VEGFR2 (Round 1)
