#!/usr/bin/env bash
# U3 (Round-11): assert that the evidence verifier's TIPS array contains
# EXACTLY the six expected keys, each mapped to a bundle name (its file
# name in 06-git-bundles/) and that file's pinned tip SHA-1.
#
# Usage: test_tips_array.sh <verifier.sh> [extraction-dir]
#   <verifier.sh>      path to the verifier.sh under test
#   [extraction-dir]   optional directory containing 06-git-bundles/
#                      (enables bundle-file existence checks)
#
# Exit 0 = all assertions pass. Requires bash only.
#
# Anti-display-artifact note: the key "main.bundle" is bracketed as
# [main.bundle] in the verifier source. Some terminal/rendering layers
# swallow the two bytes "[m" (ESC[m = ANSI reset), which makes the line
# LOOK like "ain.bundle]=..." on screen. This test therefore never relies
# on visual inspection: it parses the bytes, validates a strict per-entry
# grammar, and evaluates the array in a subshell. A genuinely corrupted
# entry (a line starting "ain.bundle]=" with no opening bracket) fails
# the grammar check and the key-set check.
set -uo pipefail

VERIFIER="${1:?usage: test_tips_array.sh <verifier.sh> [extraction-dir]}"
EXTRACT_DIR="${2:-}"

EXPECTED_KEYS="chore-hygiene-2026-10-19.bundle
integration-p0-hardening-2026-10-19.bundle
main.bundle
policy-engineering-constitution.bundle
science-law-validity-wp15.bundle
science-p0-correctness-2026-10-19.bundle"

# Pinned tips (fresh-fetch verified 2026-10-19T21:48:59Z, this session).
declare -A EXPECTED_TIPS=(
  [main.bundle]="2195a048e87dcec4aa63e8dc7981b3bb6d237ca0"
  [science-p0-correctness-2026-10-19.bundle]="ebc1fd79f20b9a6b6a262a70ce6942ba069d19e9"
  [science-law-validity-wp15.bundle]="296f10c4f2c87df2003c1a9d0150edc0b54a13d9"
  [chore-hygiene-2026-10-19.bundle]="abc450ea45183e9cebf3da8a6e20878d2d18a1bf"
  [policy-engineering-constitution.bundle]="4993cb3b6954ffb8ecad281817b2c59f1f8fc041"
  [integration-p0-hardening-2026-10-19.bundle]="0eae6d575cd91376588fe6b36d1b9212352df7fc"
)

FAIL=0
note() { printf '%s\n' "$*"; }
bad()  { printf 'FAIL: %s\n' "$*"; FAIL=1; }

[ -f "$VERIFIER" ] || { bad "verifier not found: $VERIFIER"; exit 1; }

# --- 1. extract the TIPS block --------------------------------------------
BLOCK="$(sed -n '/^declare -A TIPS=(/,/^)$/p' "$VERIFIER")"
[ -n "$BLOCK" ] || { bad "no 'declare -A TIPS=( ... )' block found"; exit 1; }

# --- 2. grammar check: every entry line must be  [key]="40-hex"  ----------
ENTRY_RE='^[[:space:]]*\[[A-Za-z0-9._-]+\.bundle\]="[0-9a-f]{40}"$'
ln=0
badlines=0
while IFS= read -r line; do
  ln=$((ln+1))
  case "$line" in
    'declare -A TIPS=(' | ')' | '') continue ;;
  esac
  if ! printf '%s' "$line" | grep -qE "$ENTRY_RE"; then
    bad "entry line $ln violates grammar: $line"
    badlines=$((badlines+1))
  fi
done <<EOF
$BLOCK
EOF
if [ "$badlines" -eq 0 ]; then
  note "grammar: all entry lines match  [name.bundle]=\"40-hex-tip\"  exactly"
  note "grammar: no line begins with a bracket-less fragment (corruption) "
fi

# --- 3. evaluate the array in a clean subshell ----------------------------
TIPS_KEYS="$(eval "$BLOCK" >/dev/null 2>&1; eval "
$(printf '%s\n' "$BLOCK" | sed 's/declare -A TIPS=(/declare -A TIPS=(/' )
for k in \"\${!TIPS[@]}\"; do echo \"\$k\"; done" 2>/dev/null | LC_ALL=C sort)"
[ -n "$TIPS_KEYS" ] || { bad "could not evaluate TIPS array"; exit 1; }

N_KEYS="$(printf '%s\n' "$TIPS_KEYS" | wc -l | tr -d ' ')"
if [ "$N_KEYS" -eq 6 ]; then
  note "key count: exactly 6"
else
  bad "key count is $N_KEYS, expected 6"
fi

EXP_KEYS_SORTED="$(printf '%s\n' "$EXPECTED_KEYS" | LC_ALL=C sort)"
if [ "$TIPS_KEYS" = "$EXP_KEYS_SORTED" ]; then
  note "key set: equals the six expected bundle names (no extra, none missing)"
else
  bad "key set differs from expected:"
  printf '%s\n' "--- expected ---" "$EXP_KEYS_SORTED" "--- actual ---" "$TIPS_KEYS"
fi
if printf '%s\n' "$TIPS_KEYS" | grep -qx 'ain.bundle'; then
  bad "mangled key 'ain.bundle' present (missing opening bracket)"
fi

# --- 4. per-key: value is 40-hex AND equals the pinned tip -----------------
for k in $TIPS_KEYS; do
  v="$(eval "$(printf '%s\n' "$BLOCK")" >/dev/null 2>&1; eval "
$(printf '%s\n' "$BLOCK")
printf '%s' \"\${TIPS[$k]}\"" 2>/dev/null)"
  if ! printf '%s' "$v" | grep -qE '^[0-9a-f]{40}$'; then
    bad "value of $k is not a 40-hex SHA-1: $v"
    continue
  fi
  if [ "$v" = "${EXPECTED_TIPS[$k]}" ]; then
    note "tip pin: $k -> $v  OK"
  else
    bad "tip pin: $k -> $v  != pinned ${EXPECTED_TIPS[$k]}"
  fi
done

# --- 5. bundle files exist (when extraction dir given) ----------------------
if [ -n "$EXTRACT_DIR" ]; then
  for k in $TIPS_KEYS; do
    if [ -f "$EXTRACT_DIR/06-git-bundles/$k" ]; then
      note "bundle file present: 06-git-bundles/$k"
    else
      bad "bundle file missing: $EXTRACT_DIR/06-git-bundles/$k"
    fi
  done
else
  note "(bundle-file existence skipped: no extraction dir given)"
fi

note ""
if [ "$FAIL" -eq 0 ]; then
  note "RESULT: TIPS ARRAY TEST PASSED — exactly six keys, each a bundle"
  note "name mapped to its pinned tip; no mangled entries."
else
  note "RESULT: TIPS ARRAY TEST FAILED (see above)"
fi
exit "$FAIL"
