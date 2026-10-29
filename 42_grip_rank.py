#!/usr/bin/env python3
"""
RANKING (not a gate): among designs that passed gate 2, rank by grip.

Method: /mnt/project-files/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md
  CA-CA within 8 A  = grip (real contact)
  CA-CA within 10 A = footprint (includes grazing)
  plus coverage of the intended epitope, i.e. the hotspots we aimed at.

Grip is a property of the BACKBONE, not of the sequence, so backbones are ranked
and only those with at least one sequence under 2 A are considered. Ranking
before gating is how this project once selected unfoldable binders; this script
refuses to rank anything that has not passed.

CPU, login node, seconds.
    cd ~/rfd && python3 42_grip_rank.py              # all targets
    cd ~/rfd && python3 42_grip_rank.py pd1

Contacts are normalised by binder length, because an unnormalised count ranks
partly by how long the chain is.
"""

import os
import sys

ROOT = os.path.expanduser("~/rfd")
GRIP_CUT = 8.0
FOOT_CUT = 10.0


def read_cas(path):
    """Return {chain: [(resnum, (x,y,z)), ...]} for CA atoms."""
    out = {}
    with open(path) as fh:
        for l in fh:
            if l[:4] != "ATOM" or len(l) < 54:
                continue
            if l[12:16] != " CA ":
                continue
            ch = l[21]
            try:
                num = int(l[22:26])
                xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError:
                continue
            out.setdefault(ch, []).append((num, xyz))
    return out


def d2(p, q):
    return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2


def backbones_that_passed(name):
    """Backbone names with at least one sequence under the gate 2 cut."""
    pf = os.path.join(ROOT, "gate2", "%s.pass" % name)
    if not os.path.exists(pf):
        return None
    seen = {}
    for line in open(pf):
        desc = line.strip()
        if not desc:
            continue
        bb = desc
        for suffix in ("_af2pred",):
            if bb.endswith(suffix):
                bb = bb[: -len(suffix)]
        # strip the _s<k> sequence index
        if "_s" in bb:
            head, _, tail = bb.rpartition("_s")
            if tail.isdigit():
                bb = head
        seen[bb] = seen.get(bb, 0) + 1
    return seen


def hotspots_for(name):
    table = os.path.join(ROOT, "inputs", "targets.tsv")
    if not os.path.exists(table):
        return []
    for line in open(table):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6 and f[0] == name:
            out = []
            for h in f[5].split(","):
                h = h.strip()
                if len(h) > 1 and h[1:].isdigit():
                    out.append(int(h[1:]))
            return out
    return []


def measure(path, hotspots):
    cas = read_cas(path)
    if "A" not in cas or "B" not in cas:
        return None
    target, binder = cas["A"], cas["B"]
    n = len(binder)
    grip = foot = 0
    hot_hit = set()
    g2, f2 = GRIP_CUT ** 2, FOOT_CUT ** 2
    for bnum, bxyz in binder:
        for tnum, txyz in target:
            dd = d2(bxyz, txyz)
            if dd <= f2:
                foot += 1
                if dd <= g2:
                    grip += 1
                    if tnum in hotspots:
                        hot_hit.add(tnum)
    return {
        "n": n,
        "grip": grip,
        "foot": foot,
        "grip_per_res": grip / float(n) if n else 0.0,
        "hot_hit": sorted(hot_hit),
    }


def main(argv):
    table = os.path.join(ROOT, "inputs", "targets.tsv")
    if argv:
        names = argv
    elif os.path.exists(table):
        names = [l.split("\t")[0] for l in open(table) if l.strip()]
    else:
        print("no targets.tsv and no target named")
        return 1

    for name in names:
        print("=" * 76)
        passed = backbones_that_passed(name)
        if passed is None:
            print("%s: no gate2/%s.pass - run 41_gate2_report.sh first" % (name, name))
            continue
        if not passed:
            print("%s: nothing passed gate 2, nothing to rank" % name)
            continue
        hotspots = hotspots_for(name)
        print("%s: %d backbone(s) passed gate 2   hotspots aimed at: %s"
              % (name, len(passed), hotspots or "unknown"))

        rows = []
        for bb, nseq in sorted(passed.items()):
            path = os.path.join(ROOT, "outputs", "%s.pdb" % bb)
            if not os.path.exists(path):
                print("  %-28s backbone file missing: %s" % (bb, path))
                continue
            m = measure(path, hotspots)
            if m is None:
                print("  %-28s needs both chain A and chain B" % bb)
                continue
            m["bb"] = bb
            m["nseq"] = nseq
            rows.append(m)

        if not rows:
            continue
        rows.sort(key=lambda r: -r["grip_per_res"])
        print("  %-26s %4s %5s %6s %6s %s"
              % ("backbone", "n", "seqs", "grip", "grip/n", "hotspots contacted"))
        for r in rows:
            cov = ""
            if hotspots:
                cov = "%d/%d  %s" % (len(r["hot_hit"]), len(hotspots),
                                     ",".join(str(h) for h in r["hot_hit"]))
            print("  %-26s %4d %5d %6d %6.3f  %s"
                  % (r["bb"], r["n"], r["nseq"], r["grip"], r["grip_per_res"], cov))

        best = rows[0]
        print("  -> best grip: %s (%.3f per residue)" % (best["bb"], best["grip_per_res"]))
        if hotspots and not best["hot_hit"]:
            print("  *** the best-grip design touches NONE of the hotspots -- it is")
            print("      binding somewhere else on the target, which is not the aim")

    print("=" * 76)
    print("grip is a RANKING among designs that already passed both gates.")
    print("Carry the top few into stage 3 (43_build_complex.py).")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
