#!/usr/bin/env python3
"""
62_blind_pocket.py -- find the binding groove from the receptor ALONE.

Why this exists
---------------
"Design blind" is stronger than "don't open the peptide file". There are
three ways the answer leaks:

  1. A learned generator's weights (RFdiffusion saw the whole 1YCR complex
     in training). We cannot wall that off; we can only decline to call it
     blind.
  2. Our own heads: we computed the F19/W23/L26 hotspots yesterday.
  3. The hotspot list itself: the 22-residue pocket in the gate dissection
     was defined by contacts WITH the peptide. Handing that to a designer
     as "aim here" smuggles the answer back in.

This script closes leak 3. It locates a druggable pocket on the MDM2
receptor using the receptor's own geometry only -- a LIGSITE-style
concavity scan -- with the peptide absent from the input. The pocket it
returns is what a designer is allowed to aim at.

Method (LIGSITE-lite)
---------------------
Lay a grid over the receptor. A grid point is solvent if it is outside
every atom's vdW+probe sphere. For each solvent point, scan outward along
7 directions (3 axes + 4 cubic diagonals); a direction scores a
"protein-solvent-protein" (PSP) event when protein is hit on BOTH sides
within a range. Points with many PSP events sit in a concavity. Cluster
the high-PSP points, take the largest cluster, and map it to the nearest
receptor residues. That residue set is the blind pocket.

Disclosure
----------
The script will, if given the original complex via --disclose-peptide,
print the overlap between the blind pocket and the peptide-contact pocket.
That comparison is a VALIDATION readout printed to stderr, never fed into
the pocket definition. The blind pocket is computed with the peptide
absent, full stop.
"""

import argparse
import math
import sys

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90}
APOLAR = {"C", "S"}


def read_atoms(path, chain=None):
    out = []
    for ln in open(path):
        if not ln.startswith("ATOM"):
            continue
        if chain and ln[21] != chain:
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[0]).upper()
        if el == "H":
            continue
        out.append((ln[21], int(ln[22:26]), ln[17:20].strip(),
                    ln[12:16].strip(), el,
                    (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))))
    return out


def d2(a, b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2


class Grid:
    def __init__(self, pts, cell):
        self.cell = cell
        self.b = {}
        for i, p in enumerate(pts):
            self.b.setdefault(self._i(p), []).append(i)

    def _i(self, p):
        c = self.cell
        return (int(p[0]//c), int(p[1]//c), int(p[2]//c))

    def near(self, p):
        ix, iy, iz = self._i(p)
        out = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    out += self.b.get((ix+dx, iy+dy, iz+dz), [])
        return out


DIRS = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
        (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("receptor", help="receptor-only PDB (no ligand/peptide)")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--spacing", type=float, default=1.0)
    ap.add_argument("--probe", type=float, default=1.4)
    ap.add_argument("--psp-min", type=int, default=5,
                    help="min PSP events (of 7) to call a point pocket")
    ap.add_argument("--scan-max", type=float, default=10.0,
                    help="how far to scan for the far wall, angstrom")
    ap.add_argument("--top", type=int, default=12,
                    help="how many pocket residues to report as hotspots")
    ap.add_argument("--disclose-peptide", default=None,
                    help="original complex PDB; prints overlap to stderr "
                         "as a validation readout, never used in the "
                         "pocket definition")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    atoms = read_atoms(args.receptor, args.chain)
    if not atoms:
        raise SystemExit("no atoms in chain %s of %s"
                         % (args.chain, args.receptor))
    xyz = [a[5] for a in atoms]
    rad = [VDW.get(a[4], 1.7) for a in atoms]
    ag = Grid(xyz, 8.0)

    lo = [min(p[i] for p in xyz) - 6 for i in range(3)]
    hi = [max(p[i] for p in xyz) + 6 for i in range(3)]
    sp = args.spacing

    def occupied(p):
        for j in ag.near(p):
            if d2(p, xyz[j]) < (rad[j] + args.probe) ** 2:
                return True
        return False

    # mark grid occupancy
    nx = int((hi[0]-lo[0])/sp)+1
    ny = int((hi[1]-lo[1])/sp)+1
    nz = int((hi[2]-lo[2])/sp)+1
    occ = {}
    for ix in range(nx):
        x = lo[0]+ix*sp
        for iy in range(ny):
            y = lo[1]+iy*sp
            for iz in range(nz):
                z = lo[2]+iz*sp
                occ[(ix, iy, iz)] = occupied((x, y, z))

    steps = int(args.scan_max/sp)
    pocket_pts = []
    for (ix, iy, iz), o in occ.items():
        if o:
            continue
        psp = 0
        for dx, dy, dz in DIRS:
            hit_pos = hit_neg = False
            for s in range(1, steps+1):
                if occ.get((ix+dx*s, iy+dy*s, iz+dz*s)):
                    hit_pos = True
                    break
            for s in range(1, steps+1):
                if occ.get((ix-dx*s, iy-dy*s, iz-dz*s)):
                    hit_neg = True
                    break
            if hit_pos and hit_neg:
                psp += 1
        if psp >= args.psp_min:
            pocket_pts.append((lo[0]+ix*sp, lo[1]+iy*sp, lo[2]+iz*sp))

    if not pocket_pts:
        raise SystemExit("no pocket found; lower --psp-min")

    # cluster pocket points (single linkage within ~2*spacing)
    pg = Grid(pocket_pts, 2.2*sp)
    seen = [False]*len(pocket_pts)
    clusters = []
    link2 = (2.2*sp)**2
    for i in range(len(pocket_pts)):
        if seen[i]:
            continue
        stack = [i]
        seen[i] = True
        comp = []
        while stack:
            k = stack.pop()
            comp.append(k)
            for j in pg.near(pocket_pts[k]):
                if not seen[j] and d2(pocket_pts[k], pocket_pts[j]) <= link2:
                    seen[j] = True
                    stack.append(j)
        clusters.append(comp)
    clusters.sort(key=len, reverse=True)
    best = [pocket_pts[i] for i in clusters[0]]
    cen = tuple(sum(p[i] for p in best)/len(best) for i in range(3))

    # map pocket points to nearest receptor residues
    res_atoms = {}
    for a in atoms:
        res_atoms.setdefault((a[0], a[1], a[2]), []).append(a)
    score = {}
    apol = {}
    for pp in best:
        for a in atoms:
            if d2(pp, a[5]) <= 4.5**2:
                key = (a[0], a[1], a[2])
                score[key] = score.get(key, 0) + 1
                if a[4] in APOLAR:
                    apol[key] = apol.get(key, 0) + 1
    ranked = sorted(score, key=lambda k: score[k], reverse=True)

    print("blind pocket from receptor geometry only")
    print("  pocket volume (grid): %d points, %.0f A^3"
          % (len(best), len(best)*sp**3))
    print("  pocket centre: %.1f %.1f %.1f" % cen)
    print("  %d lining residues; top %d by pocket-point contact:"
          % (len(ranked), args.top))
    hot = []
    for k in ranked[:args.top]:
        ap_frac = apol.get(k, 0)/score[k]
        hot.append(k)
        print("    %s%-4d %s  pts=%-3d apolar_frac=%.2f"
              % (k[0], k[1], k[2], score[k], ap_frac))
    hotspot_nums = sorted(k[1] for k in hot)
    print("\n  hotspot residues (receptor-derived, blind):")
    print("  " + ",".join("%s%d" % (hot[0][0], n) for n in hotspot_nums))

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("# blind pocket, receptor-only\n")
            fh.write("centre\t%.3f\t%.3f\t%.3f\n" % cen)
            fh.write("volume_A3\t%.1f\n" % (len(best)*sp**3))
            for k in ranked:
                fh.write("res\t%s\t%d\t%s\t%d\t%d\n"
                         % (k[0], k[1], k[2], score[k], apol.get(k, 0)))
        sys.stderr.write("wrote %s\n" % args.out)

    if args.disclose_peptide:
        comp = read_atoms(args.disclose_peptide)
        rec = [a for a in comp if a[0] == args.chain]
        pep = [a for a in comp if a[0] != args.chain]
        contact = set()
        for a in rec:
            for b in pep:
                if d2(a[5], b[5]) <= 4.5**2:
                    contact.add(a[1])
                    break
        blind = set(hotspot_nums)
        peptide_pocket = sorted(contact)
        inter = blind & contact
        sys.stderr.write("\n[DISCLOSURE -- validation only, not used above]\n")
        sys.stderr.write("  peptide-contact pocket (%d res): %s\n"
                         % (len(peptide_pocket),
                            ",".join(str(n) for n in peptide_pocket)))
        sys.stderr.write("  blind top-%d hotspots: %s\n"
                         % (args.top, ",".join(str(n) for n in hotspot_nums)))
        sys.stderr.write("  blind hotspots that ARE real contacts: %d/%d\n"
                         % (len(inter), len(hotspot_nums)))
        sys.stderr.write("  recall of peptide pocket by all %d lining res: "
                         "%d/%d\n"
                         % (len(ranked),
                            len(set(k[1] for k in ranked) & contact),
                            len(contact)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
