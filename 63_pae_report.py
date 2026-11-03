#!/usr/bin/env python3
"""Read out the AF2 complex scores for a group's three arms.

    cd ~/rfd && TAG=m1 python3 63_pae_report.py casp3 pcsk9
    cd ~/rfd && TAG=m1 CUT=10 python3 63_pae_report.py casp3

Reads  af2_<group>_<arm>_<TAG>_<stamp>.sc  written by 62_af2_arms.sbatch
       (the newest per arm, and it says which file it read)

Three things this script refuses to do, each of which has produced a wrong
number on this project before:

  1. COLUMNS ARE READ BY NAME, never by position. predict.py's column order is
     not a promise, and a positional read of a reordered file returns a real
     number from the wrong column with nothing to show it is wrong.

  2. THE GATE IS A FLOOR, NOT A RANKING. pae_interaction below CUT (default 10)
     is the pass test; among designs that fail it the value does not order them
     meaningfully. So this prints pass COUNTS per arm, and a mean with its
     standard error, and does not print a best-design table.

  3. A SINGLE DESIGN IS NOT A MEASUREMENT. AF2 here is not deterministic even on
     one pinned node and the noise is heavy-tailed - about one design in ten
     flips either side of the gate on a re-run, and roughly 3 units are needed
     before two single designs differ. Only a mean over about 9 designs has a
     floor near 0.5 units, so any arm with n < 9 is flagged and no arm-vs-arm
     difference smaller than its own standard errors is called.

Also refuses a score file that holds a duplicate design name, which is what an
append-only file written twice under one name looks like.

Stdlib only, CPU, login node, instant.

    cd ~/rfd && python3 -m py_compile 63_pae_report.py && echo INTACT || echo MANGLED
"""
import glob
import math
import os
import sys

ROOT = os.path.expanduser("~/rfd")
TAG = os.environ.get("TAG", "m1")
CUT = float(os.environ.get("CUT", "10"))
COL = os.environ.get("COL", "pae_interaction")
ARMS = ("wt", "mut", "ctl")
MIN_N = 9           # below this a mean has no usable floor
RESOLVE = 3.0       # units needed to separate two single designs


def is_number(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def read_sc(path, col):
    """(values, names, problems) for one column, read by NAME.

    The file is Rosetta score format: every data line starts with SCORE:, and
    the first such line is the header. `grep -c '^SCORE:'` therefore overcounts
    the designs by exactly one.
    """
    header, rows, problems = None, [], []
    for line in open(path):
        if not line.startswith("SCORE:"):
            continue
        fields = line.split()[1:]
        if header is None and fields and not is_number(fields[0]):
            header = fields
            continue
        if header is None:
            problems.append("data before the header - is this a score file?")
            return [], [], problems
        rows.append(fields)

    if header is None:
        problems.append("no header line")
        return [], [], problems
    if col not in header:
        problems.append("no column named %r; the file has: %s" % (col, " ".join(header)))
        return [], [], problems

    ci = header.index(col)
    di = header.index("description") if "description" in header else None
    vals, names = [], []
    for f in rows:
        if len(f) != len(header):
            problems.append("a row has %d fields against %d column names - row skipped"
                            % (len(f), len(header)))
            continue
        if not is_number(f[ci]):
            problems.append("non-numeric %s = %r - row skipped" % (col, f[ci]))
            continue
        vals.append(float(f[ci]))
        names.append(f[di] if di is not None else "?")

    if di is not None:
        dupes = len(names) - len(set(names))
        if dupes:
            # Refuse the file rather than warn and read it anyway: a warning
            # above a printed mean gets read as a mean.
            problems.append("DUPLICATE design names (%d) - this file was appended to twice."
                            " ARM SKIPPED. Re-run that arm with a fresh -scorefilename." % dupes)
            return [], names, problems
    return vals, names, problems


def stats(v):
    n = len(v)
    if n == 0:
        return (0, 0.0, 0.0, 0.0)
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else 0.0
    return (n, m, sd, sd / math.sqrt(n))


def newest_sc(group, arm):
    pat = os.path.join(ROOT, "af2_%s_%s_%s_*.sc" % (group, arm, TAG))
    found = sorted(glob.glob(pat))
    if not found:
        return None, 0
    return found[-1], len(found)


def main():
    groups = sys.argv[1:]
    if not groups:
        sys.exit("name the groups, e.g.: TAG=%s python3 63_pae_report.py casp3 pcsk9" % TAG)

    print("TAG=%s   column %s   gate: below %.1f (a FLOOR, not a ranking)\n" % (TAG, COL, CUT))

    for group in groups:
        print("%s" % group)
        per_arm = {}
        for arm in ARMS:
            path, nfiles = newest_sc(group, arm)
            if path is None:
                print("  %-3s no score file af2_%s_%s_%s_*.sc - has 62_af2_arms.sbatch run?"
                      % (arm, group, arm, TAG))
                continue
            vals, _, problems = read_sc(path, COL)
            for p in problems:
                print("  %-3s ! %s" % (arm, p))
            if not vals:
                continue
            if nfiles > 1:
                print("  %-3s (%d score files match; read the newest: %s)"
                      % (arm, nfiles, os.path.basename(path)))
            n, m, sd, sem = stats(vals)
            npass = sum(1 for x in vals if x < CUT)
            per_arm[arm] = vals
            flag = ""
            if n < MIN_N:
                flag = "   * n < %d, the mean has no usable floor here" % MIN_N
            print("  %-3s n=%-4d pass %3d (%4.1f%%)   mean %6.2f +/- %4.2f (sd %5.2f)   min %6.2f%s"
                  % (arm, n, npass, 100.0 * npass / n, m, sem, sd, min(vals), flag))
            print("      file: %s" % os.path.basename(path))

        # Arm-vs-arm, only where the numbers can carry it.
        if len(per_arm) == len(ARMS):
            ns = [len(per_arm[a]) for a in ARMS]
            if len(set(ns)) != 1:
                print("      ! the arms are n=%s - UNEQUAL. Run 60_equal_n_arms.py and re-score;"
                      " a pass-count comparison across unequal n compares sample sizes too."
                      % ",".join(str(x) for x in ns))
            base = stats(per_arm["wt"])
            for arm in ("mut", "ctl"):
                cur = stats(per_arm[arm])
                d = cur[1] - base[1]
                comb = math.sqrt(cur[3] ** 2 + base[3] ** 2)
                if comb == 0:
                    continue
                if cur[0] < MIN_N or base[0] < MIN_N:
                    verdict = "NOT callable, n < %d in one arm" % MIN_N
                else:
                    verdict = "resolved" if abs(d) > 2 * comb else "NOT resolved"
                print("      %s - wt: %+6.2f units (+/- %4.2f combined sem)  ->  %s"
                      % (arm, d, comb, verdict))
            print("      a difference under %.1f units is inside single-design noise;"
                  " only these means carry it, and only at equal n." % RESOLVE)
        print()


if __name__ == "__main__":
    main()
