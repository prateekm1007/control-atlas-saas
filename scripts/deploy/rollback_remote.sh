#!/usr/bin/env bash
# WP-14: remote rollback script — runs ON the EC2 host via `ssh bash -s`.
# Restores the service to the pre-deploy SHA recorded by deploy_remote.sh,
# rebuilds, and restarts. Used ONLY by the automated rollback step after a
# failed bounded health verification.
set -euo pipefail

APP_DIR="${HOME}/deployments/toscanini-demo/control-atlas-saas"
STATE_DIR="${HOME}/deployments/toscanini-demo"

PRE_DEPLOY_SHA="$(cat "${STATE_DIR}/.pre_deploy_sha" 2>/dev/null || true)"
if [ -z "${PRE_DEPLOY_SHA}" ]; then
  echo "FATAL: no pre-deploy rollback point recorded — cannot roll back safely."
  exit 1
fi

cd "${APP_DIR}"
echo "rolling back to pre-deploy SHA: ${PRE_DEPLOY_SHA}"
git fetch origin
git checkout --force "${PRE_DEPLOY_SHA}"

ACTUAL="$(git rev-parse HEAD)"
if [ "${ACTUAL}" != "${PRE_DEPLOY_SHA}" ]; then
  echo "FATAL: rollback checkout mismatch (got ${ACTUAL})"
  exit 1
fi

echo "── Rebuilding at rollback point"
docker compose down
docker compose build --no-cache
docker compose up -d
echo "rollback complete: service restored to ${ACTUAL}"
