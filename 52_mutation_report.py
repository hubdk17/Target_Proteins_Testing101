#!/usr/bin/env python3
"""Read out the three-arm stage-1 mutation experiment.

Metric, as used on IL-6 job 175120: how much binder mass the model places
against the mutated target residue - the number of binder CA atoms within 8 A of
that residue's CA, per design.

**Measured to the residue's CA, deliberately.** The wild-type arm has a full side
chain and the other two are stripped to five atoms, so a "nearest any atom of the
residue" metric would differ between arms by construction and measure our own
edit rather than the model's behaviour. The CA is present and identical in all
three arms.

Verifies the residue NAME in every single design before using it, so a
mislabelled arm cannot be read as a result.

    cd ~/rfd && python3 52_mutation_report.py            # all targets, TAG=m1
    cd ~/rfd && TAG=m1 python3 52_mutation_report.py pcsk9
"""
import os, sys, glob, math

HOME = os.path.expanduser("~")
RFD = os.path.join(HOME, "rfd")
TAG = os.environ.get("TAG", "m1")
NEAR, WIDE = 8.0, 10.0
ARMS = ("wt", "mut", "ctl")
EXPECT_NAME = {"mut": "ALA"}        # ctl and wt must keep the original name


def target_and_binder(path):
    """(target CAs by (chain,resnum) -> (name,xyz), [binder CA xyz]).

    RFdiffusion writes the fixed target as chain A with crystal numbering and the
    diffused binder as chain B, backbone-only glycine.
    """
    tgt, binder = {}, []
    for l in open(path):
        if not l.startswith("ATOM") or l[12:16].strip() != "CA":
            continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if l[21] == "B":
            binder.append(xyz)
        else:
            tgt[(l[21], l[22:26].strip())] = (l[17:20].strip(), xyz)
    return tgt, binder


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def stats(v):
    n = len(v)
    if n == 0:
        return (0, 0.0, 0.0, 0.0)
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else 0.0
    return (n, m, sd, sd / math.sqrt(n) if n else 0.0)


def welch(a, b):
    na, ma, sa, _ = stats(a)
    nb, mb, sb, _ = stats(b)
    den = math.sqrt(sa * sa / na + sb * sb / nb) if na and nb else 0.0
    return (ma - mb) / den if den else 0.0


def resnum_from_table():
    """Which residue each target's arms were built around, read back from the
    arm files rather than recomputed - the builder's choice is the record."""
    out = {}
    for p in glob.glob(os.path.join(RFD, "inputs", "*_mut_clean.pdb")):
        name = os.path.basename(p)[: -len("_mut_clean.pdb")]
        wt = os.path.join(RFD, "inputs", f"{name}_wt_clean.pdb")
        if not os.path.exists(wt):
            continue
        mut_res = {(l[21], l[22:26].strip()) for l in open(p)
                   if l.startswith("ATOM") and l[17:20].strip() == "ALA"}
        wt_res = {(l[21], l[22:26].strip()) for l in open(wt)
                  if l.startswith("ATOM") and l[17:20].strip() == "ALA"}
        only = mut_res - wt_res
        if len(only) == 1:
            key = only.pop()
            orig = {l[17:20].strip() for l in open(wt)
                    if l.startswith("ATOM") and (l[21], l[22:26].strip()) == key}
            out[name] = (key, orig.pop() if len(orig) == 1 else "?")
        else:
            print(f"  ! {name}: could not identify the mutated residue ({len(only)} candidates)")
    return out


targets = sys.argv[1:] or sorted(resnum_from_table())
picked = resnum_from_table()

print(f"TAG={TAG}   metric: binder CAs within {NEAR:.0f} A of the residue's CA\n")
for name in targets:
    if name not in picked:
        print(f"{name}: no arm files in inputs/ - run 50_build_mutation_arms.py")
        continue
    key, orig_name = picked[name]
    print(f"{name}  residue {key[0]}{key[1]}, wild type {orig_name}")

    res = {}
    for arm in ARMS:
        want = EXPECT_NAME.get(arm, orig_name)
        files = sorted(glob.glob(os.path.join(RFD, "outputs", f"{name}_{arm}_{TAG}_*.pdb")))
        near, wide, nearest, badname, nobinder = [], [], [], 0, 0
        for f in files:
            tgt, binder = target_and_binder(f)
            if key not in tgt:
                badname += 1
                continue
            got, ca = tgt[key]
            if got != want:
                badname += 1
                continue
            if not binder:
                nobinder += 1
                continue
            d = [dist(ca, b) for b in binder]
            near.append(sum(1 for x in d if x <= NEAR))
            wide.append(sum(1 for x in d if x <= WIDE))
            nearest.append(min(d))
        res[arm] = near
        n, m, sd, sem = stats(near)
        _, mw, _, _ = stats(wide)
        _, mn, _, _ = stats(nearest)
        flag = ""
        if badname:
            flag += f"   *** {badname} design(s) had the wrong residue name, EXCLUDED - arm may be mislabelled"
        if nobinder:
            flag += f"   ({nobinder} with no chain B)"
        print(f"  {arm:3} {want:3} n={n:4}  within {NEAR:.0f}A {m:5.2f} +/- {sd:4.2f} (sem {sem:4.3f})"
              f"   within {WIDE:.0f}A {mw:5.2f}   nearest {mn:5.2f} A{flag}")

    if all(len(res[a]) > 1 for a in ARMS):
        t_wm = welch(res["wt"], res["mut"])
        t_cw = welch(res["ctl"], res["wt"])
        t_cm = welch(res["ctl"], res["mut"])
        _, m_wt, _, _ = stats(res["wt"])
        _, m_mut, _, _ = stats(res["mut"])
        _, m_ctl, _, _ = stats(res["ctl"])
        span = m_mut - m_wt
        print(f"      Welch t:  wt vs mut {t_wm:+.2f}   ctl vs wt {t_cw:+.2f}   ctl vs mut {t_cm:+.2f}")

        # Only attribute an effect that exists. Dividing by a span that is itself
        # noise gives a confident-looking percentage with nothing behind it - on
        # IL-6 the span was 0.69 at t=6.15, which is what made the attribution
        # meaningful there.
        if abs(t_wm) < 2.0:
            print(f"      no wt-vs-mut effect to attribute (|t| = {abs(t_wm):.2f} < 2)."
                  f" There is nothing here for the control to explain; the % figure"
                  f" would be a ratio of two noise terms. Either this residue does not"
                  f" steer the binder, or n is too small.")
        else:
            pos = (m_ctl - m_wt) / span * 100
            verdict = ("sequence one-hot" if abs(pos) < 35 else
                       "side-chain geometry" if pos > 65 else
                       "MIXED, read it carefully")
            print(f"      the control sits {pos:.0f}% of the way from wt to mut  ->  {verdict}")
    else:
        print("      not enough designs in every arm yet")
    print()
