#!/usr/bin/env python3
"""
thread_monomer.py -- thread ProteinMPNN sequences onto generated MONOMER
backbones, for the backbone-generator designability pipeline (step 4).

Monomer version of thread_multi.py: no binder/target chains, no gate1 filter.
For each generated backbone PDB (one chain, any letter) and each MPNN-designed
sequence, write a PDB with the backbone coordinates and the designed residue
identities, ready for AF2 to fold. Threading rewrites only the residue name in
columns 18-20; fixed-column edits have silently corrupted data on this project
before, so every lookup and every line length is asserted.

  python model/thread_monomer.py --backbones ~/rfd/corpus_v1/bb_samples \
      --seqs ~/rfd/corpus_v1/mpnn_bb/seqs --out ~/rfd/corpus_v1/bb_threaded
"""

import argparse
import glob
import os
import sys

AA1TO3 = dict(zip(
    "ACDEFGHIKLMNPQRSTVWY",
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
))
assert len(AA1TO3) == 20, "amino acid table is incomplete"


def read_fasta(path):
    """MPNN writes the native sequence first, then one '>T=' header per sampled
    design. Return only the sampled designs."""
    lines = open(path).read().split("\n")
    out = []
    for i, head in enumerate(lines):
        if head.startswith(">T=") and i + 1 < len(lines):
            s = lines[i + 1].strip()
            if s:
                out.append(s)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbones", required=True)
    ap.add_argument("--seqs", required=True, help="MPNN seqs/ directory of .fa")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    bbs = sorted(glob.glob(os.path.join(args.backbones, "*.pdb")))
    if not bbs:
        raise SystemExit("no backbones in %s" % args.backbones)

    problems = []
    written = 0
    lens_seen = set()
    for bb_path in bbs:
        tag = os.path.splitext(os.path.basename(bb_path))[0]
        fa = os.path.join(args.seqs, "%s.fa" % tag)
        if not os.path.exists(fa):
            hits = glob.glob(os.path.join(args.seqs, "*%s*.fa" % tag))
            if len(hits) == 1:
                fa = hits[0]
            else:
                problems.append("%s: no unique .fa" % tag); continue
        bb = [l for l in open(bb_path) if l[:4] == "ATOM" and len(l) > 26]
        if not bb:
            problems.append("%s: no ATOM lines" % tag); continue
        resnums = sorted({int(l[22:26]) for l in bb})
        seqs = read_fasta(fa)
        if not seqs:
            problems.append("%s: zero sampled sequences" % tag); continue
        ok = 0
        for k, s in enumerate(seqs, 1):
            if len(s) != len(resnums):
                problems.append("%s: seq %d is %d aa, backbone has %d res"
                                % (tag, k, len(s), len(resnums)))
                continue
            mapping = dict(zip(resnums, s))
            lines, bad = [], False
            for l in bb:
                three = AA1TO3.get(mapping[int(l[22:26])])
                if three is None:
                    problems.append("%s seq %d unknown residue" % (tag, k))
                    bad = True; break
                nl = l[:17] + three + l[20:]
                if len(nl) != len(l):
                    problems.append("%s seq %d line length changed" % (tag, k))
                    bad = True; break
                lines.append(nl); lens_seen.add(len(nl.rstrip("\n")))
            if bad:
                continue
            with open(os.path.join(args.out, "%s_s%d.pdb" % (tag, k)), "w") as fh:
                fh.write("".join(lines)); fh.write("TER\nEND\n")
            ok += 1; written += 1
        print("  %-24s %3d res, %2d/%2d threaded" % (tag, len(resnums), ok, len(seqs)))

    print("=" * 60)
    print("threaded %d structures -> %s" % (written, args.out))
    if len(lens_seen) > 1:
        problems.append("MIXED ATOM line lengths %s - columns shifted"
                        % sorted(lens_seen))
    if problems:
        print("%d THING(S) TO LOOK AT:" % len(problems))
        for p in problems[:20]:
            print("  - %s" % p)
        return 1 if written == 0 else 0
    print("no complaints; all ATOM lines one length")
    return 0


if __name__ == "__main__":
    sys.exit(main())
