#!/usr/bin/env python3
"""Session-3 (Round-11, U1) constitution count re-derivation.

Re-derives the byte/char/word/line counts of BOTH constitutions from the
git-object dumps made this session (work/session3-read/), under the single
pinned counting convention (auditor T1 / Round-10):

  bytes  = len(data) octets
  chars  = len(utf-8 decode)
  words  = len(text.split())   (whitespace split; wc -w equivalent)
  lines  = LF count / split('\n') elements / non-blank lines
  blob   = git blob SHA-1 (header-prefixed)

Locale pinned by the caller (LC_ALL=C.UTF-8, LANG=C.UTF-8).
Prints one line per document and a MATCH/MISMATCH verdict vs the pinned
expected tuples.
"""
import hashlib
import os
import sys

EXPECTED = {
    "constitution.domain.md": (
        622, 615, 84, 11, 12, 9,
        "211ad2c5d4b40059d9f3f03de12e0701f9235bc7",
    ),
    "ENGINEERING_CONSTITUTION.md": (
        20572, 20564, 2784, 178, 179, 129,
        "b6bb35e84956d84a8c9d9e65734c212baf3ed075",
    ),
}

READ_DIR = "/home/z/my-project/work/session3-read"


def counts(path: str):
    data = open(path, "rb").read()
    text = data.decode("utf-8")
    return (
        len(data),                                   # bytes
        len(text),                                   # UTF-8 chars
        len(text.split()),                           # words (whitespace split)
        text.count("\n"),                            # LF lines
        len(text.split("\n")),                       # split('\n') elements
        len([l for l in text.split("\n") if l.strip()]),  # non-blank lines
        hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest(),
    )


def main() -> int:
    ok = True
    for name, exp in EXPECTED.items():
        got = counts(os.path.join(READ_DIR, name))
        if got == exp:
            print(f"MATCH {name}: bytes={got[0]} chars={got[1]} words={got[2]} "
                  f"lines={got[3]}/{got[4]}/{got[5]} blob={got[6]}")
        else:
            print(f"MISMATCH {name}: got {got} != expected {exp}")
            ok = False
    print("convention: bytes | UTF-8 chars | whitespace-split words | "
          "LF / split-elements / non-blank lines")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
