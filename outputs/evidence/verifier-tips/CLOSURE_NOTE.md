# TIPS-array question — closure (U3, Round-11)

**Question.** During Round-10 inspection, the verifier's TIPS array
displayed a line reading `ain.bundle]="2195a048…"`, suggesting the key
`[main.bundle]` might have lost its opening bracket bytes (`[m`) in the
shipped file. The Round-10 conclusion ("terminal display artifact") was
explicitly NOT to be relied upon. Round-11 U3 required byte-level proof.

**Method (evidence in this directory).**

1. `verifier.sh` was extracted from BOTH sealed archives
   (`tar -xzOf … ./verifier.sh`), so the bytes under test are the shipped
   bytes, not staging-tree copies:
   - r1 (SHA-256 `42d8e881…` archive): `r1-verifier.sh`,
     SHA-256 `af7c899d925e40eca56fdbcb27e69785e7737685a7bf4d662e78b58ae37904d6`
   - r2 (SHA-256 `5ac52b93…` archive): `r2-verifier.sh`,
     SHA-256 `b05972de3f49ebfb95e5f85d653f206af591f73979b2114c608c40445e95880d`
2. `cat -A` dumps of both files were saved
   (`r1_from_sealed_archive_cat-A.txt`, `r2_from_sealed_archive_cat-A.txt`).
   NOTE: when a human or an agent VIEWS these dumps through a rendering
   layer, the two bytes `[m` can still be swallowed on screen (ESC `[m` is
   an ANSI reset sequence), which is precisely how the artifact arises.
   Visual inspection is therefore explicitly NOT the proof.
3. `od -c` byte dumps of the first TIPS entry line were saved
   (`r1_od-c_TIPS_first_entry.txt`, `r2_od-c_TIPS_first_entry.txt`). Octal
   dumps space out every byte, so no escape-swallowing is possible. Both
   files show the literal bytes:
   `[ m a i n . b u n d l e ] = " 2 1 9 5 a 0 4 8 …`
   i.e. `[main.bundle]="2195a048…"` — the opening bracket and `m` ARE
   present in the shipped bytes of both archives.
4. A functional test was added (`test_tips_array.sh`, this directory) that
   never relies on display: it grammar-checks every entry line against
   `^[[:space:]]*[name.bundle]="[0-9a-f]{40}"$` (a bracket-less fragment
   such as `ain.bundle]=…` fails this grammar), evaluates the array in a
   subshell, and asserts EXACTLY the six expected keys, each mapped to a
   bundle name whose value equals that bundle's pinned tip.

**Results.**

- r1 sealed verifier: `test_run_r1_sealed.log` — **PASSED** (6/6 keys,
  6/6 tip pins, grammar clean, bundle files present).
- r2 sealed verifier: `test_run_r2_sealed.log` — **PASSED** (6/6 keys,
  6/6 tip pins, grammar clean, bundle files present).
- Full offline verifier re-run on an independent extraction of the sealed
  r2 archive: `verifier_full_run_r2_independent_extraction.log` —
  **ALL CHECKS PASSED** (manifest, bundles, pinned blobs, counting
  conventions).
- The r2 staging copy (`work/verify-r9/r2-test/verifier.sh`) hashes
  identically to the sealed r2 verifier (`b05972de…`), so staging and
  shipped bytes agree.

**Verdict.** The `ain.bundle]=` sighting is a display-layer artifact of
the `[m` byte pair, not file corruption. No key is wrong. Per U3's rule
("If any key is wrong, fix it in r2 and rebuild r2") **no fix and no
rebuild is required**: r2 remains byte-identical at SHA-256
`5ac52b9332c2fe2c2510b443c6aeb26ce785eccf61ad6d44cebce6151a3257cd`, and
r1 remains untouched at
`42d8e8819bc7024bfad1d9a511c7707046c3c2a0b3ff38a07fc49d4861297528`.

**Standing caveat for any future reader.** Any tool that renders
`[main.bundle]` through an ANSI-aware layer may show `ain.bundle]`. The
disproof is always: `od -c` the line, or run `test_tips_array.sh`.
