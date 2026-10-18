#!/usr/bin/env bash
# WP-14: remote deployment script — runs ON the EC2 host via `ssh bash -s`.
# Deploys the EXACT SHA the CI gate verified. Never resets to origin/main:
# a later push must never change what this deployment serves.
#
# Usage (from the runner): ssh ... "bash -s" -- "${DEPLOY_SHA}" < deploy_remote.sh
set -euo pipefail

DEPLOY_SHA="${1:?FATAL: deploy SHA argument required}"
APP_DIR="${HOME}/deployments/toscanini-demo/control-atlas-saas"
STATE_DIR="${HOME}/deployments/toscanini-demo"

cd "${APP_DIR}"

echo "── Enforcing clean Git state"
git fetch origin
git clean -fd

PRE_DEPLOY_SHA="$(git rev-parse HEAD)"
echo "${PRE_DEPLOY_SHA}" > "${STATE_DIR}/.pre_deploy_sha"
echo "pre-deploy SHA (rollback point): ${PRE_DEPLOY_SHA}"
echo "deploying SHA                  : ${DEPLOY_SHA}"

# Pin the working tree to the exact verified SHA (NOT origin/main).
git checkout --force "${DEPLOY_SHA}"
ACTUAL="$(git rev-parse HEAD)"
if [ "${ACTUAL}" != "${DEPLOY_SHA}" ]; then
  echo "FATAL: remote checkout SHA mismatch (got ${ACTUAL})"
  exit 1
fi

echo "── Rebuilding Ledger"
docker compose down
docker compose build --no-cache
docker compose up -d

echo "── Deployment fingerprint"
TREE="$(git rev-parse HEAD^{tree})"
BRAIN_IMAGE_ID="$(docker image inspect --format '{{.Id}}' "$(docker compose images -q brain | head -1)" 2>/dev/null || echo unavailable)"
echo "deployed SHA   : ${ACTUAL}"
echo "deployed tree  : ${TREE}"
echo "brain image ID : ${BRAIN_IMAGE_ID}"
echo "Deploy step finished for ${ACTUAL}"
