#!/usr/bin/env python3
"""
Convert BoltzGen inverse-folded designs into binder-first complex PDBs for our
gate-3 (AF2 initial-guess pae_interaction) - the SAME scoring we ran on
RFdiffusion, so the generators are judged on one independent metric.

BoltzGen output `<run>/intermediate_designs_inverse_folded/<spec>_NN.cif`:
  chain A = the de novo binder (designed, variable length)
  chain B = the fixed target

We write the binder FIRST as chain B and the target as chain A, with residue
numbers unique across chains (the binder is auto-offset clear of the target's
max), so predict.py's binder-first requirement holds and the site-check maps
target residues the same way after AF2 renumbers the pose. Uses the DESIGNED cif,
NOT refold_cif/ (the Boltz-2 refold), to keep AF2 independent of Boltz.

Which run, set name, and TARGET (env vars; defaults = MDM2, so the MDM2 runs are
unchanged):
  SET           complex-set name -> cplx_in/<SET>/ and inputs/<SET>.tsv   (default mdm2bg)
  BG_GLOB       glob under ~/rfd/outputs for the run dir; newest match wins (default bg_mdm2_bindany_*)
  TARGET_PREFIX first residues of the target chain, to identify/validate it (default MDM2 ETLVRPKPLLLK)
  TARGET_LEN    expected target residue count (default 85 = MDM2)
  TARGET_START  crystal-style start number for the target in the output (default 25 = MDM2)
Examples:
  source ~/rfd/venv_boltzgen/bin/activate
  cd ~/rfd && python3 scripts/47_bg_to_complex.py                                   # MDM2 unaimed (defaults)
  cd ~/rfd && SET=p53wt   BG_GLOB='bg_p53dbd_bindany_*' TARGET_PREFIX=SVPSQKTYQGSY \
      TARGET_LEN=194 TARGET_START=96 python3 scripts/47_bg_to_complex.py            # p53 DBD, wild-type core
  cd ~/rfd && SET=p53pock BG_GLOB='bg_p53dbd_pocket_*'  TARGET_PREFIX=SVPSQKTYQGSY \
      TARGET_LEN=195 TARGET_START=96 python3 scripts/47_bg_to_complex.py            # p53 DBD, Y220C core
"""
import glob
import os
import sys

import gemmi

ROOT = os.path.expanduser("~/rfd")
SET = os.environ.get("SET", "mdm2bg")
BG_GLOB = os.environ.get("BG_GLOB", "bg_mdm2_bindany_*")
TARGET_PREFIX = os.environ.get("TARGET_PREFIX", "ETLVRPKPLLLK")
TARGET_LEN = int(os.environ.get("TARGET_LEN", "85"))
TARGET_START = int(os.environ.get("TARGET_START", "25"))
OUT = os.path.join(ROOT, "cplx_in", SET)
# Binder numbered clear of the target's max, so numbers are unique across chains
# (predict.py's input_check requires it). MDM2 target max 109 -> 201; p53 290 -> 301.
_TMAX = TARGET_START + TARGET_LEN - 1
BINDER_START = (_TMAX // 100 + 1) * 100 + 1


def newest_run():
    dirs = [d for d in glob.glob(os.path.join(ROOT, "outputs", BG_GLOB)) if os.path.isdir(d)]
    return max(dirs, key=os.path.getmtime) if dirs else None


def seq(chain):
    return gemmi.one_letter_code([r.name for r in chain])


def copy_chain(src, cid, start):
    ch = gemmi.Chain(cid)
    n = start
    for res in src:
        nr = gemmi.Residue()
        nr.name = res.name
        nr.seqid = gemmi.SeqId(n, ' ')
        for at in res:
            if at.element.atomic_number == 1:   # skip hydrogens
                continue
            na = gemmi.Atom()
            na.name = at.name
            na.element = at.element
            na.pos = at.pos
            na.occ = at.occ
            na.b_iso = at.b_iso
            nr.add_atom(na)
        ch.add_residue(nr)
        n += 1
    return ch


def convert(cif, problems):
    idx = os.path.basename(cif).split("_")[-1].replace(".cif", "")
    st = gemmi.read_structure(cif)
    m = st[0]
    have = [c.name for c in m]
    if "A" not in have or "B" not in have:
        problems.append("%s: chains not A/B (%s)" % (idx, have))
        return 0
    binder, target = m["A"], m["B"]
    # safety: the target chain must match TARGET_PREFIX. If A looks like the target
    # and B doesn't, they are swapped relative to the expected layout.
    if not seq(target).startswith(TARGET_PREFIX) and seq(binder).startswith(TARGET_PREFIX):
        binder, target = target, binder
    if not seq(target).startswith(TARGET_PREFIX):
        problems.append("%s: neither chain matches target prefix %s" % (idx, TARGET_PREFIX))
        return 0
    nb, nt = len(binder), len(target)
    if nt != TARGET_LEN:
        problems.append("%s: target chain has %d res, expected %d" % (idx, nt, TARGET_LEN))
        return 0

    out = gemmi.Structure()
    out.add_model(gemmi.Model("1"))
    om = out[0]
    om.add_chain(copy_chain(binder, "B", BINDER_START))   # binder FIRST
    om.add_chain(copy_chain(target, "A", TARGET_START))   # target second
    out.setup_entities()
    path = os.path.join(OUT, "%s_%s.pdb" % (SET, idx))
    out.write_pdb(path)

    degen = seq(binder).startswith(TARGET_PREFIX)
    print("  %s_%s  binder %d res (B:%d..), target %d res (A:%d..)  OK%s"
          % (SET, idx, nb, BINDER_START, nt, TARGET_START,
             "   <- degenerate: binder seq ~target" if degen else ""))
    return 1


def main():
    run = newest_run()
    if not run:
        print("no boltzgen run dir matching %s under %s/outputs" % (BG_GLOB, ROOT))
        return 1
    src = os.path.join(run, "intermediate_designs_inverse_folded")
    cifs = sorted(glob.glob(os.path.join(src, "*.cif")))
    if not cifs:
        print("no design cifs in %s" % src)
        return 1
    os.makedirs(OUT, exist_ok=True)
    print("run: %s" % run)
    print("set=%s target_prefix=%s target_len=%d target_start=%d binder_start=%d"
          % (SET, TARGET_PREFIX, TARGET_LEN, TARGET_START, BINDER_START))
    print("converting %d BoltzGen designs -> %s" % (len(cifs), OUT))
    print("-" * 64)
    problems, made = [], 0
    for c in cifs:
        made += convert(c, problems)
    tsv = os.path.join(ROOT, "inputs", "%s.tsv" % SET)
    os.makedirs(os.path.dirname(tsv), exist_ok=True)
    with open(tsv, "w") as fh:
        fh.write("%s\tBoltzGen\t-\t-\t-\t-\n" % SET)
    print("-" * 64)
    print("wrote %d complex input(s); table %s" % (made, tsv))
    if problems:
        print("PROBLEMS:")
        for p in problems:
            print("  " + p)
    tag = "bg" if SET == "mdm2bg" else SET.replace("mdm2", "bg")
    print("next: TABLE=%s TAG=%s sbatch --exclude=dgxanode03,dgxanode04 --array=1 scripts/45_gate3_complex.sbatch"
          % (tsv, tag))
    return 0


if __name__ == "__main__":
    sys.exit(main())
