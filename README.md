# Toscanini OS — control-atlas-saas

Structural-governance API for AI-generated protein structures. Submit a
PDB/mmCIF structure and the service measures it against **15 pinned
physical laws**, returning a verdict: `PASS`, `VETO`, or `INDETERMINATE`.
Every verdict is reproducible and hash-anchored to the exact physics that
produced it.

- **Runtime source of truth:** `backend/tos/governance/station_sop.py`
  (`LAW_CANON`) — protected by canon hash `f24ec581b3e04732`.
- **Law documentation:** `governance/law_registry.yaml` (WP-8 registry;
  glossary text derives from it).
- **Version:** see `STATION_METADATA["version"]` in `station_sop.py`.

## What the service does

1. **Ingest** — you upload a structure or name a public source (RCSB PDB
   ID like `1CRN`, or AlphaFold DB accession). If no source can provide
   the structure, the API returns **503** and nothing is scored or
   recorded (no fabricated structures are ever measured).
2. **Measure** — deterministic geometric and chemical invariants (bond
   integrity, coverage, backbone gaps, bond angles, Ramachandran,
   clashscore, omega planarity, chirality, rotamers, chain integrity,
   residue identity, disulfide geometry) plus three heuristic scores
   (voxel occupancy LAW-155, hydrophobic burial LAW-182, packing quality
   LAW-200).
3. **Adjudicate** — deterministic failures veto; heuristic failures are
   advisory; coverage below 70% yields `INDETERMINATE` rather than a
   pass/fail.
4. **Record** — the audit result is stored, an NKG ledger line is
   appended (validated input only), and a v21 PDF dossier is generated.

## What a verdict means — and what it does not mean

A `PASS` means the structure did not violate any deterministic pinned
law at the pinned thresholds and cleared the coverage gate. It does
**not** mean the protein will fold, express, bind, or be manufacturable.
`VETO` means at least one deterministic law failed — the `verdict`
object names the primary failing law with observed vs threshold values.
`INDETERMINATE` means the measurement confidence was insufficient
(usually low pLDDT coverage), not that the structure is bad. Heuristic
law failures (LAW-155/182/200) are advisory and never by themselves
produce a `VETO`.

## Endpoints and authentication

All requests (except the three public reads below) require an API key
header: `X-API-Key: <key>`.

| Method | Path | Auth | Purpose |
|---|---|---|---|
| GET | `/health` | public | Token system + storage + canon-hash check |
| GET | `/laws` | public | The pinned law canon (reference data) |
| GET | `/definitions` | public | Score definitions (reference data) |
| POST | `/ingest` | key | Audit a structure (upload or candidate ID) |
| POST | `/v1/batch` | key | Audit a ZIP of structures |
| GET | `/audit/{audit_id}` | key | Fetch a stored audit |
| GET | `/nkg` | key | Recent NKG outcomes (whitelisted fields only) |
| POST | `/search` | key | Resolve a candidate ID to a source |
| GET | `/comparison/{o}/{r}` | key | Stored comparison metadata |
| GET | `/comparisons/by-original/{id}` | key | Refinements of an audit |
| POST | `/compare` | key | Delta analysis between two audits |
| GET | `/refinement/pdb/{audit_id}` | key | Original PDB for an audit |
| POST | `/refinement/token` | key | Single-use B1 callback token |
| POST | `/refinement/callback` | token | B1 refined-structure upload (JWT-gated) |
| POST | `/refinement/submit` | key | B2 managed GPU refinement job |
| GET | `/refinement/status/{job_id}` | key | Job state |
| GET | `/refinement/jobs/{audit_id}` | key | Jobs for an audit |
| GET | `/history/{user_email}` | key | Key-scoped refinement history |
| GET | `/usage` | key | Quota / credit status for the key |
| POST | `/refinement/internal-callback` | HMAC | GPU worker only — signed with `INTERNAL_CALLBACK_SECRET`; must also be denied at the proxy (see `docs/DEPLOY.md`) |

Authentication modes (fail-closed, WP-2):

- Keys configured (`TOSCANINI_API_KEYS` / `TOSCANINI_API_KEY`, or
  DB-backed keys): protected routes require a valid key; quota and tier
  enforcement applies.
- No keys and no escape hatch: protected routes return **503** — the
  service never silently runs open.
- `TOSCANINI_ALLOW_OPEN_MODE=1` (dev only, never production) explicitly
  enables unauthenticated open mode when no keys are configured.

`/history` is **key-scoped**: a valid key sees only comparisons created
with that key. Audit IDs are the full SHA-256 of the canonical
coordinate hash — identical structures share one ID; job IDs are full
UUIDs.

## Environment variables

Scanned from `os.environ` / `os.getenv` across `backend/`, `dashboard/`,
`gpu_worker/` (including WP-2/WP-3 additions):

| Variable | Consumer | Purpose |
|---|---|---|
| `TOSCANINI_API_KEYS` | brain | Comma-separated static API keys (or use DB-backed keys) |
| `TOSCANINI_API_KEY` | brain | Single-key variant |
| `TOSCANINI_ALLOW_OPEN_MODE` | brain | `1` = explicit dev open mode; never set in production |
| `ALLOWED_ORIGINS` | brain | Comma-separated CORS origins; empty = no cross-origin browser access |
| `INTERNAL_CALLBACK_SECRET` | brain + gpu_worker | Shared HMAC secret for `/refinement/internal-callback` (WP-3) |
| `TOSCANINI_JWT_SECRET` | brain | B1 callback-token signing secret — **required, no default** (WP-13): missing or <32 chars fails closed (503 on token routes) |
| `CALLBACK_RATE_LIMIT` | brain | Max B1 callbacks per key per hour (default 10) |
| `TOSCANINI_DATA_DIR` | brain | Writable data dir (audits, comparisons, NKG, telemetry) |
| `TOSCANINI_NKG_PATH` | brain | Optional NKG ledger path override |
| `REDIS_URL` | brain + gpu_worker | Job store / Celery broker |
| `BRAIN_URL` | gpu_worker | Brain base URL for the signed callback |
| `BACKEND_URL`, `USAGE_CACHE_TTL` | dashboard | Dashboard API target + widget cache |
| `GITHUB_ORG`, `GITHUB_REPO` | dashboard | BYOC notebook template links |

`docs/DEPLOY.md` has the full deployment checklist including the
network-layer deny rule for the internal callback.

## Run locally

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r backend/requirements-lock.txt -r backend/requirements-dev.txt
export TOSCANINI_DATA_DIR=$(mktemp -d)
export TOSCANINI_API_KEYS=dev-key-1
cd backend && uvicorn main:app --port 8000
# then: curl -H "X-API-Key: dev-key-1" -F mode=experimental -F candidate_id=1CRN \
#            -F file=@your.pdb http://localhost:8000/ingest
```

## Run tests

```bash
cd backend && python -m pytest tests -q
```

The suite is green with zero failures and zero errors; every skipped
test carries a reason string. CI (`.github/workflows/toscanini-gate.yml`)
runs the suite plus a pinned gitleaks secret scan and verifies the canon
(`f24ec581b3e04732`) and matrix
(`ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4`)
hashes on every run.

## Deploy

`.github/workflows/deploy.yml` deploys to EC2 on every push to `main`
(a merge IS a production deploy). Before deploying a WP-2/WP-3 build,
set the environment variables above on the instance — otherwise the
fail-closed design serves 503 on protected routes by design. Follow
`docs/DEPLOY.md` for the required proxy rules.

## Research materials

The research corpus (entry notes, run results, archives, and analysis
scripts) lives in the `control-atlas` research repository, not in this
product repo. Research scripts retained here for reproducibility are
under `scripts/research/` and are not imported by the product.
