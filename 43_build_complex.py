#!/usr/bin/env python3
"""
STAGE 3 input builder: assemble the designed complex for AF2 initial guess,
satisfying both requirements that otherwise fail silently.

Method: /mnt/project-files/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md

Two requirements, both read out of dl_binder_design/af2_initial_guess/predict.py
rather than guessed, and both of which corrupt results with no error if wrong:

  1. THE BINDER MUST BE THE FIRST CHAIN. `binderlen` comes from pose order and
     slices every number downstream. RFdiffusion writes the target first, so
     feeding its output straight in inverts the initial-guess mask and swaps
     every binder/target label while still printing a plausible number.

  2. EVERY RESIDUE NUMBER MUST BE UNIQUE ACROSS THE WHOLE FILE, chain letter
     ignored. The atom reader keys on the number alone, so a collision merges
     atoms into one entry. `input_check` raises on this - job 174604 died here.

So: binder chain B first, renumbered from 201, then target chain A keeping its
crystal numbering.

    cd ~/rfd && python3 43_build_complex.py              # all targets
    cd ~/rfd && python3 43_build_complex.py pd1 --top 3

Takes the threaded, gate-2-passing sequences from sc_in/<name>/ (which carry the
designed sequence on the designed backbone) and merges each with its target
chain from the original RFdiffusion output.

Writes cplx_in/<name>/<design>.pdb and verifies each file before leaving it.
"""

import glob
import os
import sys

ROOT = os.path.expanduser("~/rfd")
BINDER_START = 201


def load(path):
    return [l for l in open(path) if l[:4] == "ATOM" and len(l) > 54]


def renumber(lines, start):
    """Renumber residues consecutively from `start`, preserving residue grouping."""
    out, mapping, nxt = [], {}, start
    for l in lines:
        old = l[22:26]
        if old not in mapping:
            mapping[old] = nxt
            nxt += 1
        out.append(l[:22] + "%4d" % mapping[old] + l[26:])
    return out, mapping


def verify(path, n_binder, n_target):
    """Checks that can actually fail, run on the file as written."""
    problems = []
    lines = load(path)
    order, nums, lens = [], [], set()
    for l in lines:
        ch = l[21]
        if not order or order[-1] != ch:
            order.append(ch)
        nums.append(int(l[22:26]))
        lens.add(len(l.rstrip("\n")))

    if not order or order[0] != "B":
        problems.append("binder (chain B) is not first: chain order %s" % order)
    if len(order) != 2:
        problems.append("expected exactly two chain blocks, got %s" % order)

    # uniqueness is per residue, so compare residue-number sets across chains
    bres = {int(l[22:26]) for l in lines if l[21] == "B"}
    tres = {int(l[22:26]) for l in lines if l[21] == "A"}
    clash = bres & tres
    if clash:
        problems.append("%d residue numbers collide across chains: %s"
                        % (len(clash), sorted(clash)[:8]))
    if len(bres) != n_binder:
        problems.append("binder has %d residues, expected %d" % (len(bres), n_binder))
    if len(tres) != n_target:
        problems.append("target has %d residues, expected %d" % (len(tres), n_target))
    if len(lens) > 1:
        problems.append("mixed ATOM line lengths %s - columns shifted" % sorted(lens))
    return problems


def build_for_target(name, top, problems):
    threaded = sorted(glob.glob(os.path.join(ROOT, "sc_in", name, "*.pdb")))
    if not threaded:
        problems.append("%s: no threaded designs in sc_in/%s" % (name, name))
        return 0

    # keep only designs that passed gate 2, in gate 2's order if available
    pf = os.path.join(ROOT, "gate2", "%s.pass" % name)
    if os.path.exists(pf):
        allowed = set()
        for line in open(pf):
            d = line.strip()
            if d.endswith("_af2pred"):
                d = d[: -len("_af2pred")]
            if d:
                allowed.add(d)
        kept = [p for p in threaded
                if os.path.basename(p)[:-4] in allowed]
        if not kept:
            problems.append("%s: none of the threaded designs are in gate2/%s.pass"
                            % (name, name))
            return 0
        threaded = kept
    else:
        problems.append("%s: no gate2/%s.pass - building WITHOUT the gate 2 filter"
                        % (name, name))

    if top:
        threaded = threaded[:top]

    outdir = os.path.join(ROOT, "cplx_in", name)
    os.makedirs(outdir, exist_ok=True)
    made = 0
    print("-" * 72)
    print("%s: building %d complex input(s)" % (name, len(threaded)))

    for tp in threaded:
        design = os.path.basename(tp)[:-4]
        # the backbone this design came from, which still holds the target chain
        bb = design
        if "_s" in bb:
            head, _, tail = bb.rpartition("_s")
            if tail.isdigit():
                bb = head
        src = os.path.join(ROOT, "outputs", "%s.pdb" % bb)
        if not os.path.exists(src):
            problems.append("%s: backbone %s missing for %s" % (name, src, design))
            print("  %-34s NO BACKBONE FILE" % design)
            continue

        binder = [l for l in load(tp) if l[21] == "B"]
        target = [l for l in load(src) if l[21] == "A"]
        if not binder or not target:
            problems.append("%s: %s missing chain B or chain A" % (name, design))
            print("  %-34s MISSING A CHAIN" % design)
            continue

        nb = len({l[22:26] for l in binder})
        nt = len({l[22:26] for l in target})

        # target keeps its crystal numbering; the binder moves out of its way
        tnums = {int(l[22:26]) for l in target}
        start = BINDER_START
        if max(tnums) >= BINDER_START:
            # crystal numbering reaches into the binder's range, so push further out
            start = max(tnums) + 50
            print("  %-34s target numbering reaches %d, binder starts at %d"
                  % (design, max(tnums), start))
        binder_rn, _ = renumber(binder, start)

        out = os.path.join(outdir, "%s.pdb" % design)
        with open(out, "w") as fh:
            fh.writelines(binder_rn)       # BINDER FIRST
            fh.write("TER\n")
            fh.writelines(target)
            fh.write("TER\nEND\n")

        bad = verify(out, nb, nt)
        if bad:
            for b in bad:
                problems.append("%s/%s: %s" % (name, design, b))
            print("  %-34s WROTE BUT FAILED CHECKS" % design)
            for b in bad:
                print("      %s" % b)
        else:
            print("  %-34s binder %d res (from %d), target %d res  OK"
                  % (design, nb, start, nt))
            made += 1
    return made


def main(argv):
    top = None
    args = []
    i = 0
    while i < len(argv):
        if argv[i] == "--top" and i + 1 < len(argv):
            top = int(argv[i + 1])
            i += 2
            continue
        args.append(argv[i])
        i += 1

    table = os.path.join(ROOT, "inputs", "targets.tsv")
    if args:
        names = args
    elif os.path.exists(table):
        names = [l.split("\t")[0] for l in open(table) if l.strip()]
    else:
        print("no targets.tsv and no target named")
        return 1

    problems, total = [], 0
    for name in names:
        total += build_for_target(name, top, problems)

    print("=" * 72)
    print("wrote %d verified complex input(s) into cplx_in/<target>/" % total)
    print("binder is chain B and comes FIRST; numbering is unique across chains.")
    if problems:
        print("\n%d THING(S) TO LOOK AT:" % len(problems))
        for p in problems:
            print("  - %s" % p)
        return 1 if total == 0 else 0
    print("no complaints")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
