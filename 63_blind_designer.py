#!/usr/bin/env python3
"""
63_blind_designer.py -- our own de novo binder designer. Blind by construction.

What makes it blind
-------------------
It has NO trained parameters. Nothing in it was fit to any structure, so it
cannot have memorised a known binder. Its only inputs are:
  * the receptor coordinates (peptide already removed),
  * the blind pocket from 62_blind_pocket.py (receptor geometry only).
It never reads a binder sequence or a bound-peptide pose. The amino acids it
chooses are decided by the groove's shape, not by us and not by data.

The one prior it does use, stated plainly: it assumes the binder is an ideal
alpha-helix. That is a *general* inductive bias (helices bind grooves), not
MDM2-specific knowledge. It is NOT told the natural binder is a helix; it is
told to try helices because that is the simplest rigid secondary structure.
A later version can add strands and loops.

How it works
-----------
1. Build an ideal poly-alanine alpha-helix of length L by NeRF, with CB on
   every residue placed at L-chirality (so our own chirality gate confirms
   it afterwards -- a closed loop).
2. Rigid-body search: translate the helix so its centre sits near the pocket
   centre, then sweep helix-axis orientation and rotation about the axis.
   Score each pose on backbone/CB geometry only:
       + CB atoms whose side chain would point INTO the groove
       - hard clash of any helix atom with the receptor
   Keep the best clash-free poses.
3. Blind sequence design on the best pose: for each residue, measure how
   buried its side-chain direction is in the receptor. Deeply buried -> a
   large hydrophobic sized to the local cavity; solvent-facing -> charged or
   polar. No packing, no rotamers: crude on purpose ("no matter how bad").
4. Write the designed binder as chain B of a complex, ready for the gate
   battery (60) and for AF2, plus a report of the burial/hydrophobicity
   pattern so we can see whether an amphipathic i,i+4,i+7 motif emerged on
   its own.

Usage
-----
  python3 63_blind_designer.py blind/mdm2_receptor.pdb \
      --pocket blind/pocket.tsv --length 14 --out-prefix blind/design \
      --n-designs 5
"""

import argparse
import math
import random

# ideal backbone geometry
BL = {"N_CA": 1.458, "CA_C": 1.525, "C_N": 1.329, "C_O": 1.231, "CA_CB": 1.521}
BA = {"N_CA_C": 111.2, "CA_C_N": 116.2, "C_N_CA": 121.7,
      "CA_C_O": 120.8, "N_CA_CB": 110.5, "C_N_CA_CB": 0}
PHI, PSI, OMEGA = -57.0, -47.0, 180.0

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}

# hydrophobic palette, roughly by size, for buried positions
HYDRO_BY_SIZE = ["A", "V", "L", "I", "M", "F", "Y", "W"]
# polar/charged for exposed positions
POLAR = ["E", "K", "Q", "R", "S", "D", "N", "T", "H"]
# a side-chain reach (CB -> approx side-chain tip), angstrom, by one-letter
REACH = {"A": 0.0, "V": 1.5, "L": 2.6, "I": 2.5, "M": 3.3, "F": 3.4,
         "Y": 4.7, "W": 4.3, "E": 3.6, "K": 4.5, "Q": 3.7, "R": 5.5,
         "S": 1.3, "D": 2.5, "N": 2.5, "T": 1.5, "H": 3.5}
ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL"}


def vsub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def vadd(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])


def vscale(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    return (a[0]/n, a[1]/n, a[2]/n) if n else (0.0, 0.0, 0.0)


def place(a, b, c, bond, angle, dihedral):
    """NeRF: position atom d given a-b-c and internal coords (deg)."""
    ang = math.radians(180.0 - angle)
    dih = math.radians(dihedral)
    d2 = (bond*math.cos(ang), bond*math.cos(dih)*math.sin(ang),
          bond*math.sin(dih)*math.sin(ang))
    bc = unit(vsub(c, b))
    n = unit(cross(vsub(b, a), bc))
    nbc = cross(n, bc)
    M = (bc, nbc, n)
    return (c[0]+M[0][0]*d2[0]+M[1][0]*d2[1]+M[2][0]*d2[2],
            c[1]+M[0][1]*d2[0]+M[1][1]*d2[1]+M[2][1]*d2[2],
            c[2]+M[0][2]*d2[0]+M[1][2]*d2[1]+M[2][2]*d2[2])


def build_helix(n):
    """Ideal poly-A alpha-helix. Returns list of residue dicts N/CA/C/O/CB."""
    # seed first three backbone atoms
    N0 = (0.0, 0.0, 0.0)
    CA0 = (BL["N_CA"], 0.0, 0.0)
    C0 = place((-1.0, 1.0, 0.0), N0, CA0, BL["CA_C"], BA["N_CA_C"], PSI)
    atoms = [("N", 0, N0), ("CA", 0, CA0), ("C", 0, C0)]
    coords = {(0, "N"): N0, (0, "CA"): CA0, (0, "C"): C0}
    for i in range(1, n):
        a = coords[(i-1, "N")]
        b = coords[(i-1, "CA")]
        c = coords[(i-1, "C")]
        # the dihedral each NeRF step sets ends at the atom being placed:
        #   N(i)  sets  N(i-1)-CA(i-1)-C(i-1)-N(i)  = psi(i-1)
        #   CA(i) sets  CA(i-1)-C(i-1)-N(i)-CA(i)   = omega
        #   C(i)  sets  C(i-1)-N(i)-CA(i)-C(i)      = phi(i)
        N = place(a, b, c, BL["C_N"], BA["CA_C_N"], PSI)
        CA = place(b, c, N, BL["N_CA"], BA["C_N_CA"], OMEGA)
        C = place(c, N, CA, BL["CA_C"], BA["N_CA_C"], PHI)
        coords[(i, "N")] = N
        coords[(i, "CA")] = CA
        coords[(i, "C")] = C
    res = []
    for i in range(n):
        N = coords[(i, "N")]
        CA = coords[(i, "CA")]
        C = coords[(i, "C")]
        # CB via NeRF off C-N-CA, L-chirality dihedral -122.6
        CB = place(C, N, CA, BL["CA_CB"], BA["N_CA_CB"], -122.6)
        # carbonyl O: dihedral N(i)-CA(i)-C(i)-O = psi + 180 (trans to next N)
        O = place(N, CA, C, BL["C_O"], BA["CA_C_O"], PSI + 180.0)
        res.append({"N": N, "CA": CA, "C": C, "O": O, "CB": CB})
    return res


def centroid(res, atom="CA"):
    pts = [r[atom] for r in res]
    n = len(pts)
    return tuple(sum(p[i] for p in pts)/n for i in range(3))


def rot_matrix(axis, theta):
    axis = unit(axis)
    x, y, z = axis
    c = math.cos(theta)
    s = math.sin(theta)
    C = 1-c
    return (
        (c+x*x*C, x*y*C-z*s, x*z*C+y*s),
        (y*x*C+z*s, c+y*y*C, y*z*C-x*s),
        (z*x*C-y*s, z*y*C+x*s, c+z*z*C))


def apply_rt(res, R, t):
    out = []
    for r in res:
        nr = {}
        for k, p in r.items():
            nr[k] = (R[0][0]*p[0]+R[0][1]*p[1]+R[0][2]*p[2]+t[0],
                     R[1][0]*p[0]+R[1][1]*p[1]+R[1][2]*p[2]+t[1],
                     R[2][0]*p[0]+R[2][1]*p[1]+R[2][2]*p[2]+t[2])
        out.append(nr)
    return out


def d2f(a, b):
    return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2


class Grid:
    def __init__(self, pts, cell):
        self.cell = cell
        self.b = {}
        for i, p in enumerate(pts):
            self.b.setdefault((int(p[0]//cell), int(p[1]//cell),
                               int(p[2]//cell)), []).append(i)

    def near(self, p):
        c = self.cell
        ix, iy, iz = int(p[0]//c), int(p[1]//c), int(p[2]//c)
        out = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    out += self.b.get((ix+dx, iy+dy, iz+dz), [])
        return out


def read_receptor(path, chain):
    atoms = []
    for ln in open(path):
        if not ln.startswith("ATOM") or ln[21] != chain:
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[0]).upper()
        if el == "H":
            continue
        atoms.append({"el": el, "name": ln[12:16].strip(),
                      "resnum": int(ln[22:26]),
                      "xyz": (float(ln[30:38]), float(ln[38:46]),
                              float(ln[46:54])),
                      "r": VDW.get(el, 1.7)})
    return atoms


def read_pocket(path):
    centre = None
    res = []
    for ln in open(path):
        p = ln.split("\t")
        if p[0] == "centre":
            centre = (float(p[1]), float(p[2]), float(p[3]))
        elif p[0] == "res":
            res.append(int(p[2]))
    return centre, res


def pose_score(helix, rec_xyz, rec_grid, pocket_centre):
    """Backbone/CB geometry score. Returns (clash, cb_in_groove, buried)."""
    clash = 0
    cb_in = 0
    for r in helix:
        for k in ("N", "CA", "C", "O", "CB"):
            p = r[k]
            for j in rec_grid.near(p):
                dd = d2f(p, rec_xyz[j])
                if dd < 2.2**2:
                    clash += 1
        # does this residue's side chain point toward the groove?
        cbdir = unit(vsub(r["CB"], r["CA"]))
        probe = vadd(r["CB"], vscale(cbdir, 2.5))
        near_rec = sum(1 for j in rec_grid.near(probe)
                       if d2f(probe, rec_xyz[j]) < 5.0**2)
        if near_rec >= 3:
            cb_in += 1
    return clash, cb_in


def burial_of(point, rec_xyz, rec_grid, radius=6.0):
    return sum(1 for j in rec_grid.near(point)
               if d2f(point, rec_xyz[j]) < radius*radius)


def design_sequence(helix, rec_xyz, rec_grid):
    """Assign amino acids from burial alone. Returns seq, per-res report."""
    seq = []
    report = []
    for i, r in enumerate(helix):
        cbdir = unit(vsub(r["CB"], r["CA"]))
        # burial sampled where the side-chain tip would reach
        probe_close = vadd(r["CB"], vscale(cbdir, 1.5))
        probe_far = vadd(r["CB"], vscale(cbdir, 3.5))
        b_close = burial_of(probe_close, rec_xyz, rec_grid, 5.0)
        b_far = burial_of(probe_far, rec_xyz, rec_grid, 5.0)
        buried = (b_close + b_far)
        if buried >= 8:
            # pick the largest hydrophobic that fits without clashing
            choice = "A"
            for aa in HYDRO_BY_SIZE:
                tip = vadd(r["CB"], vscale(cbdir, REACH[aa]))
                hard = min((d2f(tip, rec_xyz[j])
                            for j in rec_grid.near(tip)), default=99)
                if hard > 2.6**2:
                    choice = aa
                else:
                    break
            label = "buried"
        elif buried >= 4:
            choice = "L" if b_close >= 3 else "T"
            label = "rim"
        else:
            choice = "E" if i % 2 == 0 else "K"
            label = "exposed"
        seq.append(choice)
        report.append((i, buried, b_close, b_far, choice, label))
    return "".join(seq), report


def write_complex(path, rec_path, chain, helix, seq, binder_chain="B",
                  start=1):
    lines = [ln for ln in open(rec_path)
             if ln.startswith("ATOM") and ln[21] == chain]
    serial = 10000
    out = list(lines)
    for i, (r, aa) in enumerate(zip(helix, seq)):
        rn = ONE_TO_THREE[aa]
        atoms = [("N", r["N"]), ("CA", r["CA"]), ("C", r["C"]),
                 ("O", r["O"])]
        if aa != "G":
            atoms.append(("CB", r["CB"]))
        for name, p in atoms:
            serial += 1
            el = name[0]
            out.append(
                "ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                "          %2s\n"
                % (serial, name, rn, binder_chain, start+i,
                   p[0], p[1], p[2], el))
    out.append("END\n")
    with open(path, "w") as fh:
        fh.writelines(out)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("receptor")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--pocket", required=True)
    ap.add_argument("--length", type=int, default=14)
    ap.add_argument("--n-designs", type=int, default=5)
    ap.add_argument("--orientations", type=int, default=600)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-prefix", default="design")
    args = ap.parse_args(argv)

    rec = read_receptor(args.receptor, args.chain)
    rec_xyz = [a["xyz"] for a in rec]
    rec_grid = Grid(rec_xyz, 6.0)
    centre, pocket_res = read_pocket(args.pocket)
    if centre is None:
        raise SystemExit("pocket file has no centre line")

    rng = random.Random(args.seed)
    helix0 = build_helix(args.length)
    hc = centroid(helix0, "CA")
    helix0 = [{k: vsub(p, hc) for k, p in r.items()} for r in helix0]

    # receptor centroid, to lift the helix OUT of the surface and rest it
    # on the groove rather than through it
    rec_cen = tuple(sum(p[i] for p in rec_xyz)/len(rec_xyz)
                    for i in range(3))
    outward = unit(vsub(centre, rec_cen))

    poses = []
    for _ in range(args.orientations):
        # random helix axis orientation
        ax = unit((rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1)))
        th = rng.uniform(0, 2*math.pi)
        R = rot_matrix(ax, th)
        # lift the helix centroid outward along the surface normal and sweep,
        # plus a small lateral jitter so it can slide along the groove
        for off in (3.0, 4.5, 6.0, 7.5, 9.0):
            lat = (rng.gauss(0, 1.5), rng.gauss(0, 1.5), rng.gauss(0, 1.5))
            t = vadd(vadd(centre, vscale(outward, off)), lat)
            placed = apply_rt(helix0, R, t)
            clash, cb_in = pose_score(placed, rec_xyz, rec_grid, centre)
            if clash == 0 and cb_in >= 2:
                poses.append((cb_in, placed))
    if not poses:
        raise SystemExit("no clash-free pose with groove contact found; "
                         "try --length shorter or more --orientations")
    poses.sort(key=lambda x: x[0], reverse=True)

    print("blind designer: %d clash-free poses, best groove-contact = %d CB"
          % (len(poses), poses[0][0]))
    seen_seq = set()
    made = 0
    for cb_in, placed in poses:
        if made >= args.n_designs:
            break
        seq, report = design_sequence(placed, rec_xyz, rec_grid)
        if seq in seen_seq:
            continue
        seen_seq.add(seq)
        made += 1
        path = "%s_%02d.pdb" % (args.out_prefix, made)
        write_complex(path, args.receptor, args.chain, placed, seq)
        hydro = sum(1 for c in seq if c in "AVLIMFYW")
        print("\ndesign %d  (groove CB=%d)  -> %s" % (made, cb_in, path))
        print("  sequence : %s" % seq)
        print("  hydrophobic fraction: %d/%d = %.2f"
              % (hydro, len(seq), hydro/len(seq)))
        print("  burial pattern (| = buried hydrophobic face):")
        bars = "".join("|" if r[5] == "buried" else
                       ("." if r[5] == "rim" else " ") for r in report)
        print("    " + bars)
        print("    " + seq)
        buried_idx = [r[0] for r in report if r[5] == "buried"]
        if len(buried_idx) >= 2:
            spac = [buried_idx[k+1]-buried_idx[k]
                    for k in range(len(buried_idx)-1)]
            print("  buried-residue spacings: %s" % spac)
    print("\nNext: score every design in one batch, then (only after "
          "freezing) unseal the known binders.")
    print("  python3 60_gate_battery.py %s_*.pdb --target-chain %s "
          "--binder-chain B --out blind_gates.tsv"
          % (args.out_prefix, args.chain))
    return 0


if __name__ == "__main__":
    main()
