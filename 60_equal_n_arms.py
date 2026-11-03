#!/usr/bin/env python3
"""Truncate the three arms of a group to EQUAL n before stages 2 and 3.

The arms of one group come back with different pool sizes (job 176028's pcsk9
arms landed 147 / 141 / 144 of 150) and gate 1 then passes a different number
from each. Running stages 2 and 3 on unequal pools makes every later comparison
a comparison of sample sizes as well as of arms, and the bigger arm wins extra
passes for free. So: take the smallest passing count in the group and cut every
arm in the group to exactly that many.

    cd ~/rfd && TAG=m1 bash 32_gate1_compactness.sh casp3_wt casp3_mut casp3_ctl
    cd ~/rfd && TAG=m1 python3 60_equal_n_arms.py casp3 pcsk9

Reads   gate1/<group>_<arm>.pass     written by 32_gate1_compactness.sh
Writes  gate1/<group>_<arm>.eqn      the truncated list, one path per line
        inputs/arms_<TAG>.tsv        array table for 33/62, one arm per line
        gate1/equal_n_<TAG>.tsv      the record: pool, pass, kept, per arm

Selection is the lowest design index first, which is deterministic and
independent of gate-1 rho - picking the best-rho subset would select on a
quantity that differs between arms and quietly bias the comparison. Nothing is
deleted; the .pass files stay as written.

The full-pool and full-pass counts are recorded before truncation, so gate-1
pass RATE is still read from the whole pool. Only stages 2 and 3 see the cut.

Stdlib only, CPU, login node, instant. Compile it before running so a mangled
paste fails loudly:

    cd ~/rfd && python3 -m py_compile 60_equal_n_arms.py && echo INTACT || echo MANGLED
"""
import os
import re
import sys

ROOT = os.path.expanduser("~/rfd")
TAG = os.environ.get("TAG", "m1")
ARMS = ("wt", "mut", "ctl")


def design_index(path):
    """The design number RFdiffusion put at the end of the filename.

    outputs/<group>_<arm>_<TAG>_<i>.pdb -> i. Sorting on the string would put
    design 10 before design 2, which is not wrong but is not reproducible across
    a top-up either; the integer is.
    """
    m = re.search(r"_(\d+)\.pdb$", os.path.basename(path))
    return int(m.group(1)) if m else -1


def read_pass(group, arm):
    path = os.path.join(ROOT, "gate1", "%s_%s.pass" % (group, arm))
    if not os.path.exists(path):
        return None, path
    out = []
    for line in open(path):
        p = line.strip()
        if not p:
            continue
        if not os.path.exists(p):
            print("  ! listed but missing on disk, dropped: %s" % p)
            continue
        out.append(p)
    return out, path


def pool_count(group, arm):
    """Designs that exist for this arm, whatever gate 1 made of them."""
    d = os.path.join(ROOT, "outputs")
    pre = "%s_%s_%s_" % (group, arm, TAG)
    if not os.path.isdir(d):
        return 0
    return sum(1 for f in os.listdir(d) if f.startswith(pre) and f.endswith(".pdb"))


def main():
    groups = sys.argv[1:]
    if not groups:
        sys.exit("name the groups to process, e.g.: TAG=%s python3 60_equal_n_arms.py casp3 pcsk9" % TAG)

    os.makedirs(os.path.join(ROOT, "gate1"), exist_ok=True)
    os.makedirs(os.path.join(ROOT, "inputs"), exist_ok=True)

    record = []
    table = []
    failed = []

    for group in groups:
        print("%s  (TAG=%s)" % (group, TAG))
        lists = {}
        ok = True
        for arm in ARMS:
            got, path = read_pass(group, arm)
            if got is None:
                print("  ! no %s - run 32_gate1_compactness.sh for %s_%s first" % (path, group, arm))
                ok = False
                continue
            if not got:
                print("  ! %s is empty - nothing passed gate 1 for %s_%s" % (path, group, arm))
                ok = False
                continue
            lists[arm] = sorted(got, key=design_index)
        if not ok:
            failed.append(group)
            print()
            continue

        keep = min(len(v) for v in lists.values())
        for arm in ARMS:
            got = lists[arm]
            kept = got[:keep]
            out = os.path.join(ROOT, "gate1", "%s_%s.eqn" % (group, arm))
            with open(out, "w") as fh:
                fh.write("".join(p + "\n" for p in kept))
            # Verify what ARRIVED on disk, not what was meant to be written.
            back = [l.strip() for l in open(out) if l.strip()]
            assert len(back) == keep, "wrote %d paths to %s but it holds %d" % (keep, out, len(back))
            pool = pool_count(group, arm)
            print("  %-3s pool %4d   gate1 pass %4d   kept %4d   dropped %4d"
                  % (arm, pool, len(got), keep, len(got) - keep))
            record.append((group, arm, pool, len(got), keep))
            table.append("%s_%s" % (group, arm))

        spread = max(len(v) for v in lists.values()) - keep
        if spread > 0.5 * keep:
            print("  ! the arms differ by %d of %d - truncating this far throws away most of an arm."
                  " Top the short arm up with a stage-1 resubmit (same TAG tops up) before scoring."
                  % (spread, keep))
        print()

    if not record:
        sys.exit("nothing to write - fix the groups above first")

    rec_path = os.path.join(ROOT, "gate1", "equal_n_%s.tsv" % TAG)
    with open(rec_path, "w") as fh:
        fh.write("group\tarm\tpool\tgate1_pass\tkept\n")
        for row in record:
            fh.write("%s\t%s\t%d\t%d\t%d\n" % row)

    tsv_path = os.path.join(ROOT, "inputs", "arms_%s.tsv" % TAG)
    with open(tsv_path, "w") as fh:
        fh.write("".join(n + "\n" for n in table))

    print("equal-n lists: gate1/<group>_<arm>.eqn")
    print("record:        %s" % rec_path)
    print("array table:   %s   (%d lines -> --array=1-%d)" % (tsv_path, len(table), len(table)))
    if failed:
        print("LEFT OUT: %s" % " ".join(failed))
    print()
    print("next: cd ~/rfd && TABLE=%s PASSEXT=eqn TAG=%s sbatch --array=1-%d%%2 \\" % (tsv_path, TAG, len(table)))
    print("        --exclude=dgxanode03 33_mpnn_multi.sbatch")


if __name__ == "__main__":
    main()
