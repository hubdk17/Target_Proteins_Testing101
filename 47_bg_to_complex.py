#!/usr/bin/env python3
"""
Convert BoltzGen inverse-folded designs into binder-first complex PDBs for our
gate-3 (AF2 initial-guess pae_interaction) - the SAME scoring we ran on
RFdiffusion, so the generators are judged on one independent metric.

BoltzGen output `<run>/intermediate_designs_inverse_folded/<spec>_NN.cif`:
  chain A = the de novo binder (designed, variable length 60-100)
  chain B = MDM2 (fixed target, 85 res, seq ETLVRPKPLLLK...)

We write the binder FIRST as chain B (renumbered from 201) and MDM2 as chain A
(renumbered from 25), residue numbers unique across chains - matching
43_build_complex, so predict.py's binder-first requirement holds and
46_mdm2_chimerax maps the hotspots the same way after AF2 renumbers the pose.
Uses the DESIGNED cif, NOT refold_cif/ (the Boltz-2 refold), to keep AF2
independent of Boltz's own structure prediction.

Which run to convert, and what to call the complex set (env vars):
  SET      complex-set name -> cplx_in/<SET>/ and inputs/<SET>.tsv  (default mdm2bg)
  BG_GLOB  glob under ~/rfd/outputs for the boltzgen run dir; newest match wins
           (default bg_mdm2_bindany_*)
Examples:
  source ~/rfd/venv_boltzgen/bin/activate        # gemmi lives here
  # unaimed baseline (defaults):
  cd ~/rfd && python3 scripts/47_bg_to_complex.py
  # aimed p53-cleft run:
  cd ~/rfd && SET=mdm2cleft BG_GLOB='bg_mdm2_cleft_*' python3 scripts/47_bg_to_complex.py
"""
import glob
import os
import sys

import gemmi

ROOT = os.path.expanduser("~/rfd")
SET = os.environ.get("SET", "mdm2bg")
BG_GLOB = os.environ.get("BG_GLOB", "bg_mdm2_bindany_*")
OUT = os.path.join(ROOT, "cplx_in", SET)
MDM2_PREFIX = "ETLVRPKPLLLK"   # MDM2 p53-domain start, to confirm the target chain
BINDER_START = 201
TARGET_START = 25


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
    # safety: the target chain must be MDM2. If A looks like MDM2 and B doesn't,
    # they are swapped relative to the expected layout.
    if not seq(target).startswith(MDM2_PREFIX) and seq(binder).startswith(MDM2_PREFIX):
        binder, target = target, binder
    if not seq(target).startswith(MDM2_PREFIX):
        problems.append("%s: neither chain matches MDM2 prefix" % idx)
        return 0
    nb, nt = len(binder), len(target)
    if nt != 85:
        problems.append("%s: MDM2 chain has %d res, expected 85" % (idx, nt))
        return 0

    out = gemmi.Structure()
    out.add_model(gemmi.Model("1"))
    om = out[0]
    om.add_chain(copy_chain(binder, "B", BINDER_START))   # binder FIRST
    om.add_chain(copy_chain(target, "A", TARGET_START))   # MDM2 second
    out.setup_entities()
    path = os.path.join(OUT, "%s_%s.pdb" % (SET, idx))
    out.write_pdb(path)

    degen = seq(binder).startswith(MDM2_PREFIX)
    print("  %s_%s  binder %d res (B:201..), MDM2 %d res (A:25..)  OK%s"
          % (SET, idx, nb, nt, "   <- degenerate: binder seq ~MDM2" if degen else ""))
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
    print("converting %d BoltzGen designs (set=%s) -> %s" % (len(cifs), SET, OUT))
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
