#!/usr/bin/env python3
"""Build the three-arm stage-1 mutation experiment for every batch-1 target.

The experiment, as settled on IL-6 (job 175120): a renamed mutant changes TWO
things at once, because the PDB parser keys atoms to the residue name - rename
TRP to ALA and only ALA's five atoms are read, whatever else is in the file. So
the mutant arm cannot separate "the model read a different sequence" from "the
model saw a smaller residue". The control that CAN is the complement.

    arm  residue name   side-chain atoms   isolates
    wt   TRP            full               baseline
    mut  ALA            stripped           sequence + geometry together
    ctl  TRP            stripped           geometry alone

mut and ctl carry IDENTICAL atoms and differ only in the residue name, so the
difference between them is the sequence one-hot and nothing else. On IL-6 the
control landed on top of the wild type (1.27 vs 1.26, Welch t 0.03) while the
mutant moved 0.69, so there the whole effect was the one-hot.

Picks the residue to mutate from the target's own hotspots, by side-chain size,
and prints what it chose - nothing here is hard-coded per target. Writes
inputs/<target>_<arm>_clean.pdb plus inputs/mut_arms.tsv, which is a drop-in
table for 31_rfd_multi.sbatch:

    cd ~/rfd && python3 50_build_mutation_arms.py
    cd ~/rfd && TABLE=$HOME/rfd/inputs/mut_arms.tsv TAG=m1 \
                NUM_DESIGNS=150 sbatch --array=1-18%4 31_rfd_multi.sbatch

Stdlib only, no GPU, runs on the login node.
"""
import os, sys

HOME = os.path.expanduser("~")
IN = os.path.join(HOME, "rfd", "inputs")
TABLE_IN = os.path.join(IN, "targets.tsv")
TABLE_OUT = os.path.join(IN, "mut_arms.tsv")

KEEP = ("N", "CA", "C", "O", "CB")          # ALA's complete heavy-atom set
# Preference order when several hotspots are bulky. The IL-6 experiment was run
# on a tryptophan, so a TRP keeps this the same experiment rather than a new one.
PREFER = {"TRP": 0, "TYR": 1, "PHE": 2, "ARG": 3, "LYS": 4, "HIS": 5, "MET": 6}
TOO_SMALL = {"GLY", "ALA", "SER", "CYS", "PRO"}


def read_atoms(path):
    return [l for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def residue_key(line):
    return (line[21], line[22:26].strip(), line[26])


def sidechain_atoms(lines):
    return [l for l in lines if l[12:16].strip() not in ("N", "CA", "C", "O", "OXT")]


def pick_residue(atoms, hotspots):
    """The bulkiest hotspot with a complete side chain. Returns (key, resname)."""
    by_res = {}
    for l in atoms:
        by_res.setdefault(residue_key(l), []).append(l)

    cands = []
    for hs in hotspots:
        chain, num = hs[0], hs[1:]
        key = (chain, num, " ")
        if key not in by_res:
            print(f"  ! hotspot {hs} not found in the file - skipped")
            continue
        lines = by_res[key]
        rn = lines[0][17:20].strip()
        n_sc = len(sidechain_atoms(lines))
        if rn in TOO_SMALL or n_sc < 2:
            print(f"  - {hs} is {rn} ({n_sc} side-chain atoms) - too small to strip")
            continue
        cands.append((PREFER.get(rn, 9), -n_sc, hs, key, rn, n_sc))

    if not cands:
        return None
    cands.sort()
    _, _, hs, key, rn, n_sc = cands[0]
    others = ", ".join(f"{c[2]} {c[4]}" for c in cands[1:])
    print(f"  chose {hs} = {rn}, {n_sc} side-chain atoms" + (f"  (also eligible: {others})" if others else ""))
    return key, rn, hs


def write_arm(src_atoms, out, key, newname, strip):
    """newname=None keeps the name. strip=True keeps only KEEP atoms."""
    n_kept = n_dropped = 0
    with open(out, "w") as f:
        for l in src_atoms:
            if residue_key(l) == key:
                if strip and l[12:16].strip() not in KEEP:
                    n_dropped += 1
                    continue
                if newname:
                    # Fixed columns: resname is 18-20 (1-based). Replace in place,
                    # padded to exactly three, never by splitting the line.
                    l = l[:17] + f"{newname:>3}" + l[20:]
                n_kept += 1
            f.write(l)
        f.write("TER\nEND\n")
    return n_kept, n_dropped


def main():
    if not os.path.exists(TABLE_IN):
        sys.exit(f"{TABLE_IN} missing - run 02_prep_targets.py first")

    rows, out_rows = [], []
    for line in open(TABLE_IN):
        if line.strip():
            rows.append(line.rstrip("\n").split("\t"))

    for row in rows:
        name, pdbid, contig, hotspots = row[0], row[1], row[4], row[5]
        src = os.path.join(IN, f"{name}_clean.pdb")
        print(f"\n{name} ({pdbid})")
        if not os.path.exists(src):
            print(f"  ! {src} missing - skipped")
            continue

        atoms = read_atoms(src)
        hs_list = [h.strip() for h in hotspots.split(",") if h.strip()]
        picked = pick_residue(atoms, hs_list)
        if not picked:
            print("  ! no hotspot with a strippable side chain - target skipped")
            continue
        key, rn, hs = picked

        paths = {}
        for arm, newname, strip in (("wt", None, False), ("mut", "ALA", True), ("ctl", None, True)):
            out = os.path.join(IN, f"{name}_{arm}_clean.pdb")
            kept, dropped = write_arm(atoms, out, key, newname, strip)
            paths[arm] = out
            print(f"    {arm:3}  {hs} kept {kept} atom(s), dropped {dropped}")
            out_rows.append([f"{name}_{arm}", pdbid, "", "", contig, hotspots])

        # --- integrity checks, each able to fail ---
        wt, mut, ctl = (read_atoms(paths[a]) for a in ("wt", "mut", "ctl"))

        if len(wt) != len(atoms):
            sys.exit(f"  FAIL {name}: wt arm is not a faithful copy ({len(wt)} vs {len(atoms)})")
        if len(mut) != len(ctl):
            sys.exit(f"  FAIL {name}: mut and ctl must hold identical atoms, got {len(mut)} vs {len(ctl)}")

        # The whole experiment rests on this: mut and ctl differ ONLY in the name.
        diffs = [(a, b) for a, b in zip(mut, ctl) if a != b]
        bad = [(a, b) for a, b in diffs if a[:17] + a[20:] != b[:17] + b[20:]]
        if bad:
            sys.exit(f"  FAIL {name}: mut and ctl differ outside the residue-name columns:\n    {bad[0][0]}    {bad[0][1]}")
        if len(diffs) != len(KEEP):
            print(f"    ! mut/ctl differ on {len(diffs)} lines, expected {len(KEEP)} - check for altlocs")
        names = {a[17:20].strip() for a, _ in diffs}

        # Integrity, not content: one unique ATOM line length means no line was
        # shifted. This is indifferent to what the contents should be, which is
        # why it catches the failure you did not predict - a lookup miss once
        # wrote residue names three characters short here and shifted every
        # coordinate column left by three.
        for arm, lines in (("wt", wt), ("mut", mut), ("ctl", ctl)):
            lens = sorted({len(l.rstrip(chr(10))) for l in lines})
            if len(lens) != 1:
                sys.exit(f"  FAIL {name}/{arm}: ATOM lines have {len(lens)} different lengths {lens} - a column is shifted")
        print(f"    checks passed: {len(diffs)} line(s) differ, name columns only, {sorted(names)} vs {rn}; line lengths uniform")

    if not out_rows:
        sys.exit("\nnothing written - fix the errors above")

    with open(TABLE_OUT, "w") as f:
        for r in out_rows:
            f.write("\t".join(r) + "\n")
    print(f"\nwrote {len(out_rows)} rows to {TABLE_OUT}")
    print("rows, in array order:")
    for i, r in enumerate(out_rows, 1):
        print(f"  {i:2}  {r[0]}")


main()
