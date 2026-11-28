#!/usr/bin/env python3
"""
Check the gate-3 hits and write a ChimeraX scene for each. Target-agnostic.

In an *_af2pred.pdb the binder is chain A (1..L) and the target is chain B
(L+1..L+TARGET_LEN): predict.py renumbers the whole pose from 1, so the input's
binder-first order makes the binder come back as chain A, and the target's crystal
numbering is lost. We map the aimed hotspots (crystal numbers) to the renumbered
chain-B positions FROM EACH FILE (binder length read off the file), verify every
mapped residue exists, and refuse a scene if the chain layout is not exactly binder
A then target B of TARGET_LEN residues. This is the numbering trap that once
coloured the wrong residues, so the check fails loudly instead.

Env (defaults = MDM2, so MDM2 use is unchanged):
  CPLXSET      scored set -> cplx_out/<SET> and cp_<SET>_*.sc   (default mdm2)
  HOTSPOTS     comma-separated crystal numbers of the aimed site (default 61,62,67,93,99).
               EMPTY = whole-surface run: no aimed site, just report the interface.
  TARGET_FIRST crystal number of the target's first residue (default 25)
  TARGET_LEN   target residue count (default 85)
  GATE         pae_interaction cutoff for auto-pick (default 10)
Examples:
  cd ~/rfd && python3 scripts/46_mdm2_chimerax.py                                    # MDM2
  cd ~/rfd && CPLXSET=p53pock HOTSPOTS=145,147,151,155,157,220,257 TARGET_FIRST=96 TARGET_LEN=195 \
      python3 scripts/46_mdm2_chimerax.py                                            # p53 Y220C pocket
  cd ~/rfd && CPLXSET=p53wt HOTSPOTS= TARGET_FIRST=96 TARGET_LEN=194 \
      python3 scripts/46_mdm2_chimerax.py                                            # p53 whole surface
"""
import glob
import math
import os
import sys

ROOT = os.path.expanduser("~/rfd")
SET = os.environ.get("CPLXSET", "mdm2")
CPLX = os.path.join(ROOT, "cplx_out", SET)
HOTSPOTS = [int(x) for x in os.environ.get("HOTSPOTS", "61,62,67,93,99").replace(",", " ").split()]
TARGET_FIRST = int(os.environ.get("TARGET_FIRST", "25"))
TARGET_LEN = int(os.environ.get("TARGET_LEN", "85"))
GATE = float(os.environ.get("GATE", "10"))
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
    scs = sorted(glob.glob(os.path.join(ROOT, "cp_%s_*.sc" % SET)), key=os.path.getmtime)
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
    target = [a for a in ats if a[0] == "B"]
    bres = sorted({a[1] for a in binder})
    mres = sorted({a[1] for a in target})
    problems = []
    if not binder or not target:
        problems.append("missing chain A or B")
    L = len(bres)
    if bres and bres != list(range(1, L + 1)):
        problems.append("binder chain A is not 1..L (got %d..%d)" % (bres[0], bres[-1]))
    if len(mres) != TARGET_LEN:
        problems.append("target chain B has %d residues, expected %d" % (len(mres), TARGET_LEN))
    if mres and mres[0] != L + 1:
        problems.append("target chain B starts at %d, expected %d" % (mres[0], L + 1))
    if problems:
        print("  %s: LAYOUT UNEXPECTED, no scene -> %s" % (design, "; ".join(problems)))
        return None

    # crystal c -> renumbered chain B position: L + (c - TARGET_FIRST) + 1
    hot = {c: L + (c - TARGET_FIRST) + 1 for c in HOTSPOTS}
    for c, p in hot.items():
        if p not in mres:
            print("  %s: mapped hotspot %d->B:%d absent, no scene" % (design, c, p))
            return None

    print("  %s  (binder %d res, target B:%d..%d)" % (design, L, L + 1, L + TARGET_LEN))
    covered = 0
    for c in HOTSPOTS:
        p = hot[c]
        hatoms = [a[2] for a in target if a[1] == p]
        nearest = min(dist(b[2], h) for b in binder for h in hatoms)
        tag = "CONTACT" if nearest < CONTACT else ""
        if nearest < CONTACT:
            covered += 1
        print("    target %-4d (B:%-4d)  nearest binder %.2f A  %s" % (c, p, nearest, tag))
    # total interface: target residues with any heavy atom within CONTACT of the binder
    iface = set()
    for a in target:
        if any(dist(a[2], b[2]) < CONTACT for b in binder):
            iface.add(a[1])
    if HOTSPOTS:
        print("    -> %d of %d aimed hotspots contacted; %d target residues at the interface"
              % (covered, len(HOTSPOTS), len(iface)))
    else:
        print("    -> whole-surface run (no aimed site); %d target residues at the interface"
              % len(iface))
    return design, path, hot


def write_cxc(design, hot):
    lines = [
        "# %s - open this file from the folder that holds the .pdb" % design,
        "open %s_af2pred.pdb" % design,
        "hide atoms",
        "show cartoon",
        "color /A cornflower blue",        # binder
        "color /B light gray",             # target
    ]
    if hot:
        hotsel = ",".join(str(p) for p in sorted(hot.values()))
        lines += [
            "# aimed hotspots (crystal -> these renumbered ids)",
            "show /B:%s atoms" % hotsel,
            "style /B:%s stick" % hotsel,
            "color /B:%s orange" % hotsel,
        ]
    lines += [
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
    print("set=%s target B first=%d len=%d; hotspots=%s"
          % (SET, TARGET_FIRST, TARGET_LEN, HOTSPOTS or "(none, whole-surface)"))
    made = 0
    for d in designs:
        r = check_one(d)
        if r:
            write_cxc(r[0], r[2])
            made += 1
    print("=" * 64)
    print("wrote %d scene(s) into %s" % (made, CPLX))
    print("binder = chain A (blue), target = chain B (grey)"
          + ("; hotspots = orange sticks" if HOTSPOTS else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
