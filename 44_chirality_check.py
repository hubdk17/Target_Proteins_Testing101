#!/usr/bin/env python3
"""GATE 0: is every residue an L-amino acid?

Why this exists. A D-amino acid is invisible to every other check in this
pipeline. Amber carries no improper term on the CA centre, so energy, bond
lengths, RMSD, radius of gyration, DSSP and contact counts all read normal on an
inverted residue. Gate 1 (Rg), gate 2 (backbone RMSD) and gate 3
(pae_interaction) cannot see it either. Found in the MD rig: PDBFixer built five
D-residues into the 4ZQK PD-1 85-92 loop, identically in three independent runs.

Measure: the signed volume of (N, C, CB) about CA,

    V = (N-CA) x (C-CA) . (CB-CA)        positive for L, negative for D

Deliberately NOT the N-CA-C-CB improper dihedral, which swings continuously and
misled the original diagnosis. Glycine has no CB and is skipped.

Stdlib only, no GPU, runs on the login node:

    python3 44_chirality_check.py inputs/pd1_clean.pdb
    python3 44_chirality_check.py inputs/*.pdb
    python3 44_chirality_check.py outputs/pd1_h5_*.pdb      # whole batch

Exit status 1 if any D-residue is found, so it can gate a submission:

    python3 44_chirality_check.py inputs/pd1_clean.pdb && sbatch ...
"""
import sys, glob

# A real protein is overwhelmingly L. If the majority of a file comes out
# negative, either the file is mirrored or this sign convention is backwards -
# either way the right response is to stop, not to report every residue as D.
MAJORITY_WARN = 0.5


def residues(path):
    """{(chain, resseq, icode): {atom_name: (x,y,z)}} in file order, ATOM only."""
    out, order = {}, []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        if name not in ("N", "CA", "C", "CB"):
            continue
        key = (line[21], line[22:26].strip(), line[26])
        if key not in out:
            out[key] = {"resname": line[17:20].strip()}
            order.append(key)
        try:
            out[key][name] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            sys.exit(f"{path}: unparseable coordinates on:\n{line.rstrip()}")
    return [(k, out[k]) for k in order]


def signed_volume(r):
    ca = r["CA"]
    n = [r["N"][i] - ca[i] for i in range(3)]
    c = [r["C"][i] - ca[i] for i in range(3)]
    b = [r["CB"][i] - ca[i] for i in range(3)]
    cross = (n[1] * c[2] - n[2] * c[1],
             n[2] * c[0] - n[0] * c[2],
             n[0] * c[1] - n[1] * c[0])
    return sum(cross[i] * b[i] for i in range(3))


def check(path):
    res = residues(path)
    if not res:
        print(f"{path}: NO ATOM records with N/CA/C/CB - wrong file?")
        return 1

    vols, skipped, incomplete = [], 0, []
    for key, r in res:
        if r["resname"] == "GLY":
            skipped += 1
            continue
        if not all(a in r for a in ("N", "CA", "C", "CB")):
            incomplete.append((key, r["resname"], sorted(k for k in r if k != "resname")))
            continue
        vols.append((key, r["resname"], signed_volume(r)))

    if not vols:
        print(f"{path}: no residue had all four of N, CA, C, CB - nothing to test")
        return 1

    d = [v for v in vols if v[2] < 0]
    frac = len(d) / len(vols)
    print(f"{path}")
    print(f"  tested {len(vols)} residues  (glycine skipped: {skipped})")

    rc = 0
    if frac > MAJORITY_WARN:
        print(f"  *** {len(d)} of {len(vols)} residues read NEGATIVE ({frac:.0%}).")
        print("  *** A real protein is almost all L. Either this structure is")
        print("  *** mirrored, or the sign convention here is inverted. Do not")
        print("  *** read this as 'every residue is D' - stop and work out which.")
        rc = 1
    elif d:
        print(f"  *** {len(d)} D-AMINO ACID(S) - this structure is stereochemically wrong:")
        for (ch, num, ic), rn, v in d:
            print(f"        chain {ch} {rn}{num}{ic.strip()}   V = {v:+.3f} A^3")
        print("  *** Nothing downstream can see this. Do not score these designs.")
        rc = 1
    else:
        print("  all L-amino acids")

    # Report the range actually present: a loop missing from the crystal is a
    # gap here, and a gap is the safe state - it is a REBUILT loop that carries
    # the risk. 4ZQK's 85-92 is the known case.
    for ch in dict.fromkeys(k[0] for k, _ in res):
        nums = sorted(int(k[1]) for k, _ in res if k[0] == ch)
        gaps = [(a, b) for a, b in zip(nums, nums[1:]) if b - a > 1]
        g = "  gaps: " + ", ".join(f"{a}->{b}" for a, b in gaps) if gaps else "  no gaps"
        print(f"  chain {ch}: {nums[0]}-{nums[-1]}, {len(nums)} residues{g}")

    if incomplete:
        print(f"  {len(incomplete)} residue(s) missing a backbone/CB atom, untestable:")
        for key, rn, have in incomplete[:10]:
            print(f"        chain {key[0]} {rn}{key[1]}  has only {have}")
    return rc


paths = [p for a in sys.argv[1:] for p in (glob.glob(a) or [a])]
if not paths:
    sys.exit(__doc__)
bad = sum(check(p) for p in paths)
print(f"\n{len(paths)} file(s) checked, {bad} with a problem")
sys.exit(1 if bad else 0)
