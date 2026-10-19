#!/usr/bin/env bash
# U5 audit (Round-11): confirm in the current delivery tree that
#   (1) retired D-xx labels appear ONLY in INTERNAL_CROSSWALK.md
#   (2) every OD-xx ID appears with exactly one meaning
#       (targeted collision checks for the known historical collisions:
#        OD-16 must never mean metadata integrity; OD-17 must never mean
#        merge/deploy; the auditor's ad-hoc D-4/D-11 labels must be gone
#        from user-facing text)
# Scope: the CURRENT user-facing delivery tree (outputs staging). Sealed
# historical archives (r1/r2 tarballs) are out of scope by design and are
# documented as such in the audit output.
set -uo pipefail
TREE="${1:?usage: od_register_audit.sh <outputs-tree>}"
FAIL=0
note() { printf '%s\n' "$*"; }
bad()  { printf 'FAIL: %s\n' "$*"; FAIL=1; }

note "== U5 register audit of tree: $TREE =="
note "scope: current user-facing delivery; sealed r1/r2 archives excluded (historical)"

# --- 1. D-xx labels outside the internal crosswalk -------------------------
note ""
note "-- [1] retired D-xx labels outside INTERNAL_CROSSWALK.md --"
DXX="$(grep -rnE '\bD-[0-9]{1,2}\b' "$TREE" --include='*.md' \
        | grep -v 'INTERNAL_CROSSWALK.md' || true)"
if [ -z "$DXX" ]; then
  note "OK: no D-xx label appears in any user-facing document"
else
  bad "D-xx labels found outside the crosswalk:"
  printf '%s\n' "$DXX"
fi

# --- 2. OD occurrence inventory ---------------------------------------------
note ""
note "-- [2] OD-xx occurrence inventory (files mentioning each ID) --"
for i in $(seq 1 21); do
  id="$(printf 'OD-%02d' "$i")"
  files="$(grep -rl "$id" "$TREE" --include='*.md' 2>/dev/null \
           | sed "s|^$TREE/||" | LC_ALL=C sort | tr '\n' ' ')"
  note "$id: ${files:-<none>}"
done

# --- 3. targeted collision checks -------------------------------------------
note ""
note "-- [3] targeted collision checks --"
# 3a. OD-16 must not be used for metadata integrity (old coder-side collision)
#     Scope note: INTERNAL_CROSSWALK.md is the sanctioned internal annex that
#     DOCUMENTS the collision history (U5: retired labels live only there);
#     the check therefore applies to user-facing files, mirroring check [1].
hits="$(grep -rn 'OD-16' "$TREE" --include='*.md' | grep -v 'INTERNAL_CROSSWALK.md' \
        | grep -iE 'metadata|integrity' || true)"
if [ -z "$hits" ]; then
  note "OK: OD-16 never co-occurs with metadata/integrity in user-facing files"
  note "    (the collision is documented only inside INTERNAL_CROSSWALK.md, as U5 permits)"
else bad "OD-16 used near metadata/integrity:"; printf '%s\n' "$hits"; fi
# 3b. OD-17 must not be used for merge/deploy (same user-facing scope)
hits="$(grep -rn 'OD-17' "$TREE" --include='*.md' | grep -v 'INTERNAL_CROSSWALK.md' \
        | grep -iE 'merge|deploy' || true)"
if [ -z "$hits" ]; then note "OK: OD-17 never co-occurs with merge/deploy"
else bad "OD-17 used near merge/deploy:"; printf '%s\n' "$hits"; fi
# 3c. each OD-xx that appears is within OD-01..OD-21
odd="$(grep -rhoE '\bOD-[0-9]+\b' "$TREE" --include='*.md' | LC_ALL=C sort -u \
       | awk -F- '$2 < 1 || $2 > 21' || true)"
if [ -z "$odd" ]; then note "OK: every OD ID used is within OD-01..OD-21"
else bad "OD IDs outside the consolidated range:"; printf '%s\n' "$odd"; fi
# 3d. register table row count sanity: 21 decision rows
rows="$(grep -cE '^\| OD-[0-9]{2} \|' "$TREE/OWNER_DECISIONS.md" || true)"
if [ "$rows" -eq 21 ]; then note "OK: register table has exactly 21 OD rows"
else bad "register table has $rows OD rows, expected 21"; fi
# 3e. all statuses OPEN in the table
nonopen="$(grep -E '^\| OD-[0-9]{2} \|' "$TREE/OWNER_DECISIONS.md" \
           | grep -vc '\*\*OPEN\*\*' || true)"
if [ "$nonopen" -eq 0 ]; then note "OK: every register row shows status OPEN"
else bad "$nonopen register rows are not OPEN"; fi

note ""
if [ "$FAIL" -eq 0 ]; then
  note "RESULT: U5 REGISTER AUDIT PASSED"
else
  note "RESULT: U5 REGISTER AUDIT FAILED (see above)"
fi
exit "$FAIL"
