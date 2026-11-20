#!/usr/bin/env python3
"""
95_extract_frames.py -- step 1 of the from-scratch BACKBONE generator.

Goal
----
Turn our blind corpus of PDB structures into a training set of backbone frames
for an SE(3) flow-matching model (FrameDiff / FrameFlow style). Each residue is
a rigid frame built from its N, CA, C atoms; the generator later learns to flow
random frames into protein-like frames. This script only does the data side:
pull clean, contiguous, single-chain backbone segments and write their N/CA/C
coordinates. The frames (rotation + translation) are built from those three
atoms at train time -- we store coordinates here so the file is pure text and
runs on the bare login node (no numpy).

Design decisions (unconditional monomers first)
-----------------------------------------------
* One segment per structure: the chain with the LONGEST clean segment.
* "Clean / contiguous" is tested GEOMETRICALLY, not by residue numbering:
  consecutive residues must each have N, CA, C and consecutive CA-CA distance
  < --break-dist (4.5 A). This catches missing-residue gaps that residue
  numbers lie about.
* Length window [--min-len, --max-len]; segments longer than max-len are
  cropped to the first max-len residues (the loader may random-crop later).
* SELF-VERIFYING BLIND: re-runs 80's MDM2/MDMX/p53 holdout guard over every
  file at startup and ABORTS on any hit, exactly like 81. A memorised fold is
  as much a check-that-cannot-fail as a memorised sequence.
* SANITY GATE (a check that can fail): reports mean N-CA, CA-C, CA-CA(i,i+1)
  bond lengths over everything emitted and flags if they stray from physical
  values (1.46 / 1.52 / 3.80 A). Wrong column parsing shows up here loudly.

Output: one TSV, one row per emitted residue:
  name  idx  Nx Ny Nz  CAx CAy CAz  Cx Cy Cz
`name` = PDB id (file stem), unique per structure (one segment each).

Usage (login node; stdlib only, no container needed):
  cd ~/rfd/scripts
  python3 model/95_extract_frames.py ~/rfd/corpus_v1/raw \
      --keep-list ~/rfd/corpus_v1/keep_unique.txt \
      --out ~/rfd/corpus_v1/backbones.tsv
"""

import argparse
import glob
import importlib.util
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "fetch_corpus", os.path.join(_HERE, "80_fetch_corpus.py"))
fc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fc)

# ideal backbone bond lengths (angstrom), for the sanity gate
IDEAL = {"N-CA": 1.458, "CA-C": 1.525, "CA-CA": 3.803}
TOL = 0.15  # allowed deviation of the MEAN before we flag a parsing problem


def _dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                     + (a[2] - b[2]) ** 2)


def backbone_by_chain(text):
    """First model only. -> {chain: [(key, {'N':xyz,'CA':xyz,'C':xyz}), ...]}
    in file order, residues keyed by (resSeq+iCode), only standard residues
    that carry all three backbone atoms."""
    acc = {}      # chain -> dict key -> atom dict  (dedupe residues/altlocs)
    order = {}    # chain -> list of keys in first-seen order
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"):
            break
        if not ln.startswith("ATOM"):
            continue
        atom = ln[12:16].strip()
        if atom not in ("N", "CA", "C"):
            continue
        if ln[17:20].strip() not in fc.THREE2ONE:
            continue
        alt = ln[16]
        if alt not in (" ", "A"):          # keep one altloc only
            continue
        ch = ln[21]
        key = ln[22:27]                    # resSeq + iCode
        try:
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        res = acc.setdefault(ch, {}).setdefault(key, {})
        if atom not in res:                # first occurrence wins
            res[atom] = xyz
            if len(res) == 1:
                order.setdefault(ch, []).append(key)
    out = {}
    for ch, keys in order.items():
        seq = []
        for k in keys:
            r = acc[ch][k]
            if "N" in r and "CA" in r and "C" in r:
                seq.append((k, r))
        if seq:
            out[ch] = seq
    return out


def longest_segment(residues, break_dist):
    """Longest run of residues whose consecutive CA-CA distance < break_dist."""
    if not residues:
        return []
    best = cur = [residues[0]]
    for prev, nxt in zip(residues, residues[1:]):
        if _dist(prev[1]["CA"], nxt[1]["CA"]) < break_dist:
            cur.append(nxt)
        else:
            if len(cur) > len(best):
                best = cur
            cur = [nxt]
    if len(cur) > len(best):
        best = cur
    return best


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", help="corpus raw/ directory of .pdb files")
    ap.add_argument("--out", required=True, help="output TSV of backbone coords")
    ap.add_argument("--keep-list", default=None,
                    help="restrict to these filenames (90 --emit-keep output)")
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=128)
    ap.add_argument("--break-dist", type=float, default=4.5,
                    help="CA-CA distance (A) above which a chain is cut")
    ap.add_argument("--mdm2-cut", type=float, default=0.30)
    ap.add_argument("--p53-cut", type=float, default=0.60)
    args = ap.parse_args(argv)

    if os.path.exists(args.out):
        raise SystemExit("%s exists; pass a fresh --out (no silent append)."
                         % args.out)

    files = sorted(glob.glob(os.path.join(args.raw, "*.pdb")))
    if not files:
        raise SystemExit("no .pdb files in %s" % args.raw)
    if args.keep_list:
        keep = set(l.strip() for l in open(args.keep_list) if l.strip())
        files = [f for f in files if os.path.basename(f) in keep]
        print("keep-list: %d of the corpus files" % len(files))

    # --- blind guard over the real data (abort on any MDM2/MDMX/p53 hit) ---
    print("leak-check: scanning %d files for held-out targets ..." % len(files))
    for p in files:
        try:
            cseq = fc.chains_from_pdb(open(p).read())
        except Exception:
            continue
        hold, why = fc.is_holdout(cseq, args.mdm2_cut, args.p53_cut)
        if hold:
            raise SystemExit("LEAK: %s -> %s. Corpus not blind; aborting."
                             % (os.path.basename(p), why))
    print("leak-check: 0 of %d held-out -> corpus clean" % len(files))

    n_emit = n_skip_short = n_no_bb = 0
    len_hist = []
    s_nca = s_cac = s_caca = 0.0
    c_nca = c_cac = c_caca = 0
    with open(args.out, "w") as out:
        out.write("name\tidx\tNx\tNy\tNz\tCAx\tCAy\tCAz\tCx\tCy\tCz\n")
        for p in files:
            name = os.path.splitext(os.path.basename(p))[0]
            try:
                chains = backbone_by_chain(open(p).read())
            except Exception:
                continue
            if not chains:
                n_no_bb += 1
                continue
            # best segment across all chains of this structure
            seg = max((longest_segment(r, args.break_dist)
                       for r in chains.values()), key=len, default=[])
            if len(seg) < args.min_len:
                n_skip_short += 1
                continue
            seg = seg[:args.max_len]
            len_hist.append(len(seg))
            prev_ca = None
            for idx, (_, r) in enumerate(seg):
                N, CA, C = r["N"], r["CA"], r["C"]
                s_nca += _dist(N, CA); c_nca += 1
                s_cac += _dist(CA, C); c_cac += 1
                if prev_ca is not None:
                    s_caca += _dist(prev_ca, CA); c_caca += 1
                prev_ca = CA
                out.write("%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t"
                          "%.3f\t%.3f\t%.3f\n"
                          % (name, idx, N[0], N[1], N[2],
                             CA[0], CA[1], CA[2], C[0], C[1], C[2]))
            n_emit += 1

    len_hist.sort()
    print("\n=== frame extraction: %s ===" % args.raw)
    print("structures emitted        : %d" % n_emit)
    print("skipped (segment < min)   : %d" % n_skip_short)
    print("skipped (no backbone)     : %d" % n_no_bb)
    if len_hist:
        med = len_hist[len(len_hist) // 2]
        print("segment length min/med/max: %d / %d / %d"
              % (len_hist[0], med, len_hist[-1]))
        print("total residues (frames)   : %d" % sum(len_hist))

    print("\n--- bond-length sanity gate (mean over all frames) ---")
    ok = True
    for lab, s, c in (("N-CA", s_nca, c_nca), ("CA-C", s_cac, c_cac),
                      ("CA-CA", s_caca, c_caca)):
        m = s / c if c else 0.0
        bad = abs(m - IDEAL[lab]) > TOL
        ok = ok and not bad
        print("  %-5s mean %.3f A  (ideal %.3f)%s"
              % (lab, m, IDEAL[lab], "   <-- OFF" if bad else ""))
    print("  gate: %s" % ("PASS -- geometry looks like real backbone"
                          if ok else "FAIL -- check column parsing / inputs"))
    print("\nwrote %s" % args.out)
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
