#!/usr/bin/env python3
"""P0 item 4 evidence: 1PPE insertion-code ordering, from the real bytes.

Source of truth: RCSB download of 1PPE (toyocamycin-free proteinase E /
porcine pancreatic elastase complex).
  URL:    https://files.rcsb.org/download/1PPE.pdb
  sha256: fe51b4280055199c644214049f59a752e6124d858eaea29c1a2f721fb03ed0bb
  bytes:  203958

What this script proves, from the deposited coordinates alone:

  1. 1PPE chain E lists insertion-coded residues BEFORE their base number
     in file order: 183, 184A, 184, 185 / 187, 188A, 188, 189 /
     220, 221A, 221, 222.
  2. In FILE ORDER, every consecutive residue pair around those
     insertions is joined by a genuine peptide bond (C-N ~ 1.3 A).
  3. Under the engine's previous (res_seq, icode) SORT, the same
     coordinates produce C-N distances of 4-7 A -- fabricated backbone
     gaps. This is the artifact behind the v24 report's "1PPE genuine
     7.12 A coordinate break" conclusion (now superseded).

Usage:
  python3 1ppe_evidence.py /path/to/1PPE.pdb
"""
import hashlib
import sys

EXPECTED_SHA256 = "fe51b4280055199c644214049f59a752e6124d858eaea29c1a2f721fb03ed0bb"
EXPECTED_BYTES = 203958

INSERTION_WINDOWS = [(180, 190), (216, 224)]  # chain E residue ranges


def parse_residues(lines, chain):
    """File-order residue list for one chain: [(seq, icode, resname, atoms)]."""
    residues = []
    seen = set()
    for line in lines:
        if not line.startswith("ATOM  ") and not line.startswith("HETATM"):
            continue
        if line[21:22] != chain:
            continue
        try:
            seq = int(line[22:26])
        except ValueError:
            continue
        icode = line[26:27].strip()
        key = (seq, icode)
        name = line[12:16].strip()
        pos = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if key not in seen:
            seen.add(key)
            residues.append([seq, icode, line[17:20].strip(), {}])
        for r in residues:
            if (r[0], r[1]) == key:
                r[3][name] = pos
                break
    return [tuple(r) for r in residues]


def dist(a, b):
    return sum((ai - bi) ** 2 for ai, bi in zip(a, b)) ** 0.5


def main(path):
    data = open(path, "rb").read()
    digest = hashlib.sha256(data).hexdigest()
    print(f"file: {path}")
    print(f"sha256: {digest}")
    print(f"bytes: {len(data)}")
    print(f"digest matches expected: {digest == EXPECTED_SHA256}")
    print(f"byte count matches expected: {len(data) == EXPECTED_BYTES}")
    if digest != EXPECTED_SHA256:
        print("REFUSING TO CONTINUE: bytes differ from the audited evidence file")
        return 1

    lines = data.decode("ascii", errors="ignore").splitlines()
    residues = parse_residues(lines, "E")

    print("\n== chain E insertion-code regions, FILE ORDER ==")
    for seq, icode, resname, _ in residues:
        for lo, hi in INSERTION_WINDOWS:
            if lo <= seq <= hi:
                print(f"  {seq}{icode:<2s} {resname}")
                break

    print("\n== bonded C-N distances, FILE ORDER (window pairs only) ==")
    for r1, r2 in zip(residues, residues[1:]):
        in_window = any(lo <= r1[0] <= hi or lo <= r2[0] <= hi
                        for lo, hi in INSERTION_WINDOWS)
        if in_window and r1[3].get("C") and r2[3].get("N"):
            d = dist(r1[3]["C"], r2[3]["N"])
            mark = "PEPTIDE BOND" if d <= 2.0 else "NOT BONDED"
            print(f"  {r1[0]}{r1[1]:<2s} -> {r2[0]}{r2[1]:<2s}  C-N = {d:6.3f} A   {mark}")

    sorted_res = sorted(residues, key=lambda r: (r[0], r[1]))
    print("\n== same coordinates under the REMOVED (res_seq, icode) sort ==")
    for r1, r2 in zip(sorted_res, sorted_res[1:]):
        in_window = any(lo <= r1[0] <= hi or lo <= r2[0] <= hi
                        for lo, hi in INSERTION_WINDOWS)
        if in_window and r1[3].get("C") and r2[3].get("N"):
            d = dist(r1[3]["C"], r2[3]["N"])
            mark = "peptide bond" if d <= 2.0 else "FABRICATED GAP"
            print(f"  {r1[0]}{r1[1]:<2s} -> {r2[0]}{r2[1]:<2s}  C-N = {d:6.3f} A   {mark}")

    print("\n== CA-CA across the insertion pairs (file order) ==")
    for r1, r2 in zip(residues, residues[1:]):
        if r1[1] or r2[1]:  # pair touching an insertion-coded residue
            if r1[3].get("CA") and r2[3].get("CA"):
                d = dist(r1[3]["CA"], r2[3]["CA"])
                print(f"  {r1[0]}{r1[1]:<2s} -> {r2[0]}{r2[1]:<2s}  CA-CA = {d:6.3f} A")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
