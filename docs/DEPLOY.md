# Deployment Guide — Toscanini SaaS (control-atlas-saas)

This document covers the environment variables and network-layer controls
required to run the service safely. It is the reference for WP-2 / WP-3
hardening requirements.

## Required environment variables

| Variable | Set on | Purpose | Effect if missing |
|---|---|---|---|
| `TOSCANINI_API_KEYS` | brain (EC2) | Comma-separated static API keys. DB-backed keys (`backend/tos/saas/keys.db`) are the alternative. | Protected routes return **503** (fail-closed) |
| `TOSCANINI_ALLOW_OPEN_MODE` | brain (dev only) | `1` enables unauthenticated OPEN_MODE when no keys are configured. **Never set in production.** | Open mode stays off (safe) |
| `ALLOWED_ORIGINS` | brain | Comma-separated browser origins allowed by CORS, e.g. `https://app.example.com`. | No cross-origin browser access |
| `INTERNAL_CALLBACK_SECRET` | brain **and** gpu_worker | Shared HMAC secret for `/refinement/internal-callback`. Must be identical on both sides. | Brain: callback returns 503. Worker: refuses to send (fail-closed) |
| `REDIS_URL` | brain + gpu_worker | Job store / Celery broker. | Defaults to `redis://redis:6379/0` |
| `TOSCANINI_DATA_DIR` | brain | Writable data directory (audit store, comparisons, NKG, telemetry). | Defaults to `/app/data` in Docker |

Set variables before merging WP-2/WP-3 — a deploy without
`TOSCANINI_API_KEYS` serves 503 on every protected route by design.

## Network-layer restriction on the internal callback (WP-3)

`POST /refinement/internal-callback` is authenticated with an HMAC-SHA256
signature (`X-Timestamp` / `X-Signature` over `timestamp + "." + raw_body`,
5-minute freshness window). It must ALSO be unreachable from the public
internet — defense in depth. Deny it at the load balancer or reverse proxy:

### nginx (reverse proxy in front of the brain)

```nginx
location = /refinement/internal-callback {
    # Only the GPU worker may reach this endpoint.
    allow <gpu_worker_private_ip>;
    deny all;
    proxy_pass http://brain:8000;
}
```

### AWS ALB / security groups

- Listener rule: `IF path == /refinement/internal-callback THEN return 403`
  (fixed-response), or
- Security group: restrict port 8000 on the brain instance to the gpu_worker
  security group + the ALB, and deny the path at the ALB as above.

### docker-compose (single host)

The brain's 8000 port must not be published on `0.0.0.0`. Publish it only on
the internal compose network and, if a host port is unavoidable, bind to
`127.0.0.1` and front it with the nginx rules above.

## Verifying the restriction

```bash
# from OUTSIDE the cluster — must be refused (403/404/timeout, never 200)
curl -i -X POST https://<public-host>/refinement/internal-callback \
     -H 'Content-Type: application/json' -d '{}'
```

A valid signed request is only expected from the gpu_worker container; any
200 from a public interface means the proxy rules above are not applied.
