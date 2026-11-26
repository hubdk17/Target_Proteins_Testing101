#!/usr/bin/env python3
"""
Check the gate-3 MDM2 hits and write a ChimeraX scene for each.

Reads the AF2-predicted complexes in cplx_out/mdm2/. In an *_af2pred.pdb the
binder is chain A (1..L) and MDM2 is chain B (L+1..L+85): PyRosetta renumbers the
whole pose from 1, so the input's binder-first order makes the binder come back
as chain A, and MDM2's crystal numbering (25-109) is lost. We therefore map the
p53-cleft hotspots we aimed at - MDM2 crystal residues 61,62,67,93,99 - to the
renumbered chain-B positions FROM EACH FILE (binder length read off the file),
verify every mapped residue exists, and refuse to emit a scene if the chain
layout is not exactly binder A then MDM2 B of 85 residues. This is the numbering
trap that once produced a scene colouring the wrong residues, so the check can
fail loudly instead.

The scene colours binder vs MDM2, shows the aimed hotspots as sticks, and
highlights the real interface geometrically (so it cannot be mis-numbered).

CPU, login node, seconds. Default: every design in the newest cp_mdm2_*.sc with
pae_interaction < 10. Or name designs explicitly.
    cd ~/rfd && python3 scripts/46_mdm2_chimerax.py
    cd ~/rfd && python3 scripts/46_mdm2_chimerax.py mdm2_mdm2b1_0_s4 mdm2_mdm2b1_1_s5
"""
import glob
import math
import os
import sys

ROOT = os.path.expanduser("~/rfd")
CPLX = os.path.join(ROOT, "cplx_out", "mdm2")
HOTSPOTS = [61, 62, 67, 93, 99]   # MDM2 crystal numbering, the p53 cleft we aimed at
MDM2_FIRST = 25                   # receptor chain A starts at crystal residue 25
MDM2_LEN = 85                     # 25..109 inclusive
GATE = 10.0                       # pae_interaction gate used to auto-pick hits
CONTACT = 5.0                     # heavy-atom contact cutoff (A)


def atoms(path):
    """[(chain, resnum, (x,y,z)), ...] for heavy atoms (drop hydrogens)."""
    out = []
    for l in open(path):
        if l[:4] != "ATOM" or len(l) < 54:
            continue
        elem = l[76:78].strip() or l[12:16].strip()[:1]
        if elem == "H":
            continue
        try:
            out.append((l[21], int(l[22:26]),
                        (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
        except ValueError:
            continue
    return out


def dist(p, q):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p, q)))


def pick_from_score():
    scs = sorted(glob.glob(os.path.join(ROOT, "cp_mdm2_*.sc")), key=os.path.getmtime)
    if not scs:
        return []
    hdr, rows = None, []
    for l in open(scs[-1]):
        if not l.startswith("SCORE:"):
            continue
        f = l.split()
        if "pae_interaction" in f:
            hdr = f
            continue
        if hdr is None:
            continue
        rec = dict(zip(hdr, f))
        try:
            pae = float(rec["pae_interaction"])
        except (KeyError, ValueError):
            continue
        if pae < GATE:
            rows.append((rec["description"].replace("_af2pred", ""), pae))
    rows.sort(key=lambda r: r[1])
    print("auto-picked %d design(s) under pae_interaction %.0f from %s"
          % (len(rows), GATE, os.path.basename(scs[-1])))
    return [r[0] for r in rows]


def check_one(design):
    path = os.path.join(CPLX, design + "_af2pred.pdb")
    if not os.path.exists(path):
        print("  %s: no %s" % (design, os.path.basename(path)))
        return None
    ats = atoms(path)
    binder = [a for a in ats if a[0] == "A"]
    mdm2 = [a for a in ats if a[0] == "B"]
    bres = sorted({a[1] for a in binder})
    mres = sorted({a[1] for a in mdm2})
    problems = []
    if not binder or not mdm2:
        problems.append("missing chain A or B")
    L = len(bres)
    if bres and bres != list(range(1, L + 1)):
        problems.append("binder chain A is not 1..L (got %d..%d)" % (bres[0], bres[-1]))
    if len(mres) != MDM2_LEN:
        problems.append("MDM2 chain B has %d residues, expected %d" % (len(mres), MDM2_LEN))
    if mres and mres[0] != L + 1:
        problems.append("MDM2 chain B starts at %d, expected %d" % (mres[0], L + 1))
    if problems:
        print("  %s: LAYOUT UNEXPECTED, no scene -> %s" % (design, "; ".join(problems)))
        return None

    # crystal c -> renumbered chain B position: L + (c - 25) + 1
    hot = {c: L + (c - MDM2_FIRST) + 1 for c in HOTSPOTS}
    for c, p in hot.items():
        if p not in mres:
            print("  %s: mapped hotspot %d->B:%d absent, no scene" % (design, c, p))
            return None

    # binder heavy atoms vs each hotspot's heavy atoms
    print("  %s  (binder %d res, MDM2 B:%d..%d)" % (design, L, L + 1, L + MDM2_LEN))
    covered = 0
    for c in HOTSPOTS:
        p = hot[c]
        hatoms = [a[2] for a in mdm2 if a[1] == p]
        nearest = min(dist(b[2], h) for b in binder for h in hatoms)
        tag = "CONTACT" if nearest < CONTACT else ""
        if nearest < CONTACT:
            covered += 1
        print("    MDM2 %-3d (B:%-3d)  nearest binder %.2f A  %s" % (c, p, nearest, tag))
    # total interface: MDM2 residues with any heavy atom within CONTACT of the binder
    iface = set()
    for a in mdm2:
        if any(dist(a[2], b[2]) < CONTACT for b in binder):
            iface.add(a[1])
    print("    -> %d of %d aimed hotspots contacted; %d MDM2 residues at the interface"
          % (covered, len(HOTSPOTS), len(iface)))
    return design, path, hot


def write_cxc(design, hot):
    hotsel = ",".join(str(p) for p in sorted(hot.values()))
    lines = [
        "# %s - open this file from the folder that holds the .pdb" % design,
        "open %s_af2pred.pdb" % design,
        "hide atoms",
        "show cartoon",
        "color /A cornflower blue",        # binder
        "color /B light gray",             # MDM2
        "# aimed p53-cleft hotspots (MDM2 61,62,67,93,99 -> these renumbered ids)",
        "show /B:%s atoms" % hotsel,
        "style /B:%s stick" % hotsel,
        "color /B:%s orange" % hotsel,
        "# the real interface, computed geometrically (cannot be mis-numbered)",
        "select /A & /B :<5",
        "show sel atoms",
        "style sel stick",
        "color sel gold",
        "select clear",
        "view",
    ]
    out = os.path.join(CPLX, design + ".cxc")
    open(out, "w").write("\n".join(lines) + "\n")
    print("    scene: %s" % out)


def main(argv):
    designs = argv or pick_from_score()
    if not designs:
        print("no designs: pass names or run gate 3 first")
        return 1
    print("=" * 64)
    made = 0
    for d in designs:
        r = check_one(d)
        if r:
            write_cxc(r[0], r[2])
            made += 1
    print("=" * 64)
    print("wrote %d scene(s) into %s" % (made, CPLX))
    print("binder = chain A (blue), MDM2 = chain B (grey), hotspots = orange sticks")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
