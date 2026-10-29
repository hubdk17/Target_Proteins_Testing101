#!/usr/bin/env python3
"""
GATE 2, part 1: thread each ProteinMPNN sequence back onto its own RFdiffusion
backbone, for many targets at once.

Method and thresholds: /mnt/project-files/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md
This is the multi-target version of the threading script in that document. Same
asserts, same integrity checks; it just keeps targets in separate directories so
names cannot collide at twelve targets.

Threading rewrites the residue name in columns 18-20 of each backbone ATOM line.
Fixed-column edits have silently corrupted data twice on this project, so every
lookup and every line length is asserted. Compile it before running so a mangled
paste fails loudly rather than quietly:

    cd ~/rfd && python3 -m py_compile thread_multi.py && echo INTACT || echo MANGLED
    cd ~/rfd && python3 thread_multi.py            # all targets
    cd ~/rfd && python3 thread_multi.py pd1        # one target

Reads   gate1/<name>.pass       the backbones that passed gate 1
        mpnn_out_<name>/seqs/   or mpnn_out/<name>/seqs/
Writes  sc_in/<name>/<backbone>_s<k>.pdb    one file per sequence

Only backbones listed in gate1/<name>.pass are threaded. Spending AF2 time on a
design that already failed the compactness gate is the mistake this ordering
exists to prevent.
"""

import glob
import os
import sys

AA1TO3 = dict(zip(
    "ACDEFGHIKLMNPQRSTVWY",
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
))
assert len(AA1TO3) == 20, "amino acid table is incomplete"

ROOT = os.path.expanduser("~/rfd")
BINDER_CHAIN = "B"


def seq_dir(name):
    """MPNN output location, whichever layout was used."""
    for cand in ("mpnn_out_%s/seqs" % name,
                 "mpnn_out/%s/seqs" % name,
                 "mpnn_out_%s" % name):
        full = os.path.join(ROOT, cand)
        if os.path.isdir(full) and glob.glob(os.path.join(full, "*.fa")):
            return full
    return None


def read_fasta(path):
    """MPNN writes one sequence per '>T=' header. Returns list of sequences."""
    lines = open(path).read().split("\n")
    out = []
    for i, head in enumerate(lines):
        if head.startswith(">T=") and i + 1 < len(lines):
            s = lines[i + 1].strip()
            if s:
                out.append(s)
    return out


def thread_one(name, problems):
    pass_file = os.path.join(ROOT, "gate1", "%s.pass" % name)
    if not os.path.exists(pass_file):
        problems.append("%s: no gate1/%s.pass - run 32_gate1_compactness.sh first"
                        % (name, name))
        return 0
    backbones = [l.strip() for l in open(pass_file) if l.strip()]
    if not backbones:
        problems.append("%s: nothing passed gate 1, nothing to thread" % name)
        return 0

    sdir = seq_dir(name)
    if sdir is None:
        problems.append("%s: no MPNN .fa files found (looked for mpnn_out_%s/seqs)"
                        % (name, name))
        return 0

    outdir = os.path.join(ROOT, "sc_in", name)
    os.makedirs(outdir, exist_ok=True)
    written = 0
    print("-" * 68)
    print("%s: %d backbone(s) passed gate 1, sequences from %s"
          % (name, len(backbones), os.path.relpath(sdir, ROOT)))

    for bb_path in backbones:
        tag = os.path.basename(bb_path)
        if tag.endswith(".pdb"):
            tag = tag[:-4]

        fa = os.path.join(sdir, "%s.fa" % tag)
        if not os.path.exists(fa):
            hits = glob.glob(os.path.join(sdir, "*%s*.fa" % tag))
            if len(hits) == 1:
                fa = hits[0]
            else:
                problems.append("%s: no unique .fa for backbone %s" % (name, tag))
                print("  %-28s NO SEQUENCE FILE" % tag)
                continue

        bb = [l for l in open(bb_path)
              if l[:4] == "ATOM" and len(l) > 26 and l[21] == BINDER_CHAIN]
        if not bb:
            problems.append("%s: backbone %s has no chain %s atoms"
                            % (name, tag, BINDER_CHAIN))
            print("  %-28s NO CHAIN %s" % (tag, BINDER_CHAIN))
            continue

        resnums = sorted({int(l[22:26]) for l in bb})
        seqs = read_fasta(fa)
        if not seqs:
            problems.append("%s: %s parsed zero sequences" % (name, tag))
            print("  %-28s ZERO SEQUENCES" % tag)
            continue

        ok = 0
        for k, s in enumerate(seqs, 1):
            if len(s) != len(resnums):
                # A length mismatch means the sequence does not belong to this
                # backbone. Threading it anyway would shift every residue.
                problems.append("%s: %s seq %d is %d aa but backbone has %d residues"
                                % (name, tag, k, len(s), len(resnums)))
                continue
            mapping = dict(zip(resnums, s))
            lines = []
            bad = False
            for l in bb:
                aa = mapping[int(l[22:26])]
                three = AA1TO3.get(aa)
                if three is None:
                    problems.append("%s: %s seq %d has unknown residue %r"
                                    % (name, tag, k, aa))
                    bad = True
                    break
                nl = l[:17] + three + l[20:]
                if len(nl) != len(l):
                    problems.append("%s: %s seq %d line length changed"
                                    % (name, tag, k))
                    bad = True
                    break
                lines.append(nl)
            if bad:
                continue
            out = os.path.join(outdir, "%s_s%d.pdb" % (tag, k))
            with open(out, "w") as fh:
                fh.write("".join(lines))
                fh.write("TER\nEND\n")
            ok += 1
            written += 1
        print("  %-28s %3d residues, %2d/%2d sequences threaded"
              % (tag, len(resnums), ok, len(seqs)))

    return written


def integrity(name, problems):
    """Checks that can actually fail, per the project's own practice."""
    outdir = os.path.join(ROOT, "sc_in", name)
    files = sorted(glob.glob(os.path.join(outdir, "*.pdb")))
    if not files:
        return
    lens = set()
    for f in files:
        for l in open(f):
            if l[:4] == "ATOM":
                lens.add(len(l.rstrip("\n")))
    if len(lens) > 1:
        problems.append("%s: threaded files have mixed ATOM line lengths %s "
                        "- columns are shifted" % (name, sorted(lens)))
        print("  *** MIXED LINE LENGTHS %s" % sorted(lens))

    # A blank residue name is the signature of the 2026-10-28 threading bug.
    names = set()
    for l in open(files[0]):
        if l[:4] == "ATOM":
            names.add(l[17:20])
    blank = [n for n in names if not n.strip()]
    if blank:
        problems.append("%s: blank residue names in threaded output" % name)
        print("  *** BLANK RESIDUE NAMES -- the threading table missed a letter")
    print("  residue names seen in %s: %d distinct, no blanks"
          % (os.path.basename(files[0]), len(names)))


def main(argv):
    table = os.path.join(ROOT, "inputs", "targets.tsv")
    if argv:
        names = argv
    elif os.path.exists(table):
        names = [l.split("\t")[0] for l in open(table) if l.strip()]
    else:
        print("no targets.tsv and no target named on the command line")
        return 1

    problems, total = [], 0
    for name in names:
        n = thread_one(name, problems)
        total += n
        if n:
            integrity(name, problems)

    print("=" * 68)
    print("threaded %d structure(s) into sc_in/<target>/" % total)
    if problems:
        print("\n%d THING(S) TO LOOK AT:" % len(problems))
        for p in problems:
            print("  - %s" % p)
        return 1 if total == 0 else 0
    print("no complaints")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
