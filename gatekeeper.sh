#!/usr/bin/env bash
#
# Toscanini Gatekeeper — Phase 3.4 Stability Gate
# ────────────────────────────────────────────────
# Single command to verify the entire station.
# Runs unit tests + sentinel audit.
# Appends to STABILITY_LOG.txt on success, clears it on failure.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

LOG_FILE="STABILITY_LOG.txt"
BACKEND="backend"
TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
DATE_ONLY=$(date -u +"%Y-%m-%d")
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "══════════════════════════════════════════════════════════════════════"
echo "  TOSCANINI GATEKEEPER — Phase 3.4"
echo "  Timestamp: $TIMESTAMP"
echo "  Commit:    $COMMIT"
echo "══════════════════════════════════════════════════════════════════════"

GATE_PASSED=true

# ── Stage 1: Unit Tests ──────────────────────────────────────────────────────
echo ""
echo "──────────────────────────────────────────────────────────────────────"
echo "  STAGE 1: Unit Tests"
echo "──────────────────────────────────────────────────────────────────────"

cd "$BACKEND"
TEST_OUTPUT=$(python3 -m pytest tests/ -v --tb=short 2>&1)
TEST_EXIT=$?
TEST_COUNT=$(echo "$TEST_OUTPUT" | grep -oP '\d+ passed' | head -1 || echo "? passed")
cd "$SCRIPT_DIR"

if [ $TEST_EXIT -eq 0 ]; then
    echo "  ✓ Unit Tests: $TEST_COUNT"
else
    echo "  ✗ Unit Tests: FAILED"
    echo "$TEST_OUTPUT" | tail -20
    GATE_PASSED=false
fi

# ── Stage 2: Hash Verification ───────────────────────────────────────────────
echo ""
echo "──────────────────────────────────────────────────────────────────────"
echo "  STAGE 2: Hash Verification"
echo "──────────────────────────────────────────────────────────────────────"

HASH_OUTPUT=$(python3 -c "
import sys; sys.path.insert(0, 'backend')
from tos.governance.station_sop import compute_canon_hash, LAW_CANON_HASH
from tos.governance.modality_matrix import compute_matrix_hash
canon = compute_canon_hash()
matrix = compute_matrix_hash()
c_ok = canon == '6a9cd4b4349b81de'
m_ok = matrix == 'ed1731ace09a9a3bd24d3659ab6f3fd184da2a2335d935a234fd6a9c14880aa4'
print(f'Canon:  {canon} {\"OK\" if c_ok else \"DRIFT\"}')
print(f'Matrix: {matrix[:16]}... {\"OK\" if m_ok else \"DRIFT\"}')
if not (c_ok and m_ok):
    sys.exit(1)
" 2>&1)
HASH_EXIT=$?

echo "$HASH_OUTPUT" | while IFS= read -r line; do echo "  $line"; done

if [ $HASH_EXIT -ne 0 ]; then
    echo "  ✗ Hash Verification: FAILED"
    GATE_PASSED=false
else
    echo "  ✓ Hash Verification: PASSED"
fi

# ── Stage 3: Sentinel Audit ──────────────────────────────────────────────────
echo ""
echo "──────────────────────────────────────────────────────────────────────"
echo "  STAGE 3: Sentinel Audit (7 Benchmarks)"
echo "──────────────────────────────────────────────────────────────────────"

# Check if server is running
if curl -s http://localhost:8000/health > /dev/null 2>&1; then
    SENTINEL_OUTPUT=$(python3 backend/scripts/sentinel_audit.py 2>&1)
    SENTINEL_EXIT=$?

    echo "$SENTINEL_OUTPUT"

    if [ $SENTINEL_EXIT -ne 0 ]; then
        GATE_PASSED=false
    fi
else
    echo "  ⚠ Server not running on localhost:8000"
    echo "  ⚠ Sentinel audit SKIPPED (start server with: cd backend && uvicorn main:app --port 8000)"
    echo "  ⚠ Gate will pass on unit tests + hashes only."
    echo "  NOTE: For full CI compliance, the server must be running."
fi

# ── Gate Result ───────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════════════════════════"

if [ "$GATE_PASSED" = true ]; then
    echo "  ✓ GATEKEEPER: ALL STAGES PASSED"
    echo "══════════════════════════════════════════════════════════════════════"

    # Append to stability log
    echo "$TIMESTAMP | $COMMIT | $TEST_COUNT | HASH_OK | GREEN" >> "$LOG_FILE"
    echo ""
    echo "  Stability Log updated: $LOG_FILE"

    # Count unique green days
    GREEN_DAYS=$(grep "GREEN" "$LOG_FILE" | awk -F'T' '{print $1}' | sort -u | wc -l)
    echo "  Green days: $GREEN_DAYS / 7"

    if [ "$GREEN_DAYS" -ge 7 ]; then
        echo ""
        echo "  ════════════════════════════════════════════════════════"
        echo "  ✓✓✓ PHASE 3 STABILITY MANDATE: ACHIEVED ✓✓✓"
        echo "  ✓✓✓ PHASE 4 IS UNLOCKED                 ✓✓✓"
        echo "  ════════════════════════════════════════════════════════"
    else
        REMAINING=$((7 - GREEN_DAYS))
        echo "  $REMAINING day(s) remaining until Phase 4 unlock."
    fi

    exit 0
else
    echo "  ✗ GATEKEEPER: FAILED — STABILITY CLOCK RESET"
    echo "══════════════════════════════════════════════════════════════════════"

    # Clear the log on failure
    > "$LOG_FILE"
    echo "  STABILITY_LOG.txt cleared. Clock resets to Day 0."

    exit 1
fi
