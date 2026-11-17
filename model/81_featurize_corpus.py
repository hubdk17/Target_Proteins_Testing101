#!/usr/bin/env python3
"""
81_featurize_corpus.py -- turn the blind corpus into a target-conditioned
interface-chemistry table, one row per interface residue.

The point (ranked claim 3)
--------------------------
Stage 1 of the pipeline is chemically empty on the target side: a hotspot is
one bit per residue, no chemistry. The capacity we say is missing is
TARGET-CONDITIONED interface features. This table is exactly that: for every
residue sitting in an interface, it records what that residue IS (the label
the model will learn to produce) alongside the CHEMISTRY IT FACES on the
partner -- how buried it gets, and how many donors / acceptors / apolar /
charged partner atoms are within reach.

A model trained on this learns P(residue identity | partner-pocket chemistry):
"given this local target environment, what amino acid belongs here." That is
the conditioning stage 1 lacks, learned from real interfaces, with MDM2/MDMX/
p53 already held out by 80_fetch_corpus.py.

No GPU and NO network here -- it reads the already-downloaded corpus -- so run
it as a plain CPU Slurm job (it is SASA-heavy; ~1-5 s per structure):

  sbatch --cpus-per-task=4 --mem=8G --time=06:00:00 \
    --job-name=featurize --output=logs/%x_%j.out --error=logs/%x_%j.err \
    --wrap 'python3 model/81_featurize_corpus.py ~/rfd/corpus_v1/raw \
            --out ~/rfd/corpus_v1/interface_features.tsv'

Reuses the VALIDATED SASA / parser from 60_gate_battery.py (loaded by path, so
the numbers match the gate battery exactly -- no second implementation to drift).
"""

import argparse
import glob
import importlib.util
import math
import os
import sys

# -- load the validated gate-battery code by path (its name starts with a digit,
#    so it cannot be imported normally) ---------------------------------------
_GB = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "60_gate_battery.py")
_spec = importlib.util.spec_from_file_location("gate_battery", _GB)
gb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gb)

POLAR = gb.POLAR_ELEMENTS
CHARGED_POS = gb.CHARGED_POS
CHARGED_NEG = gb.CHARGED_NEG
THREE2ONE = gb.THREE_TO_ONE


def partner_role(atom):
    """donor / acceptor / both / None for a partner polar atom (heavy-atom)."""
    el, nm, rn = atom.element, atom.name, atom.resname
    if el not in ("N", "O"):
        return None
    if nm == "N":
        return "donor"
    if nm in ("O", "OXT"):
        return "acceptor"
    if el == "O":
        return "both" if nm in ("OG", "OG1", "OH") else "acceptor"
    if rn == "HIS":
        return "both"
    return "donor"


def protein_chains(atoms, min_res=10):
    res_by_chain = {}
    for a in atoms:
        if a.resname in THREE2ONE:
            res_by_chain.setdefault(a.chain, set()).add(a.resnum)
    return {c for c, rs in res_by_chain.items() if len(rs) >= min_res}


def best_pair(atoms, chains, cutoff=5.0):
    """The chain pair with the most inter-chain heavy-atom contacts."""
    by = {c: [a for a in atoms if a.chain == c and a.resname in THREE2ONE]
          for c in chains}
    cl = sorted(chains)
    best, bestn = None, 0
    cut2 = cutoff * cutoff
    for i in range(len(cl)):
        for j in range(i + 1, len(cl)):
            a1, a2 = by[cl[i]], by[cl[j]]
            if not a1 or not a2:
                continue
            g = gb.Grid(a2, 2.0 * cutoff)
            n = 0
            for a in a1:
                for k in g.near(a.xyz):
                    if gb.dist2(a.xyz, a2[k].xyz) <= cut2:
                        n += 1
                        break
            if n > bestn:
                bestn, best = n, (cl[i], cl[j])
    return best, bestn


def residue_sasa_buried(side_atoms, pair_atoms, n_points):
    """Per-residue buried apolar/polar area for one side of the interface.

    buried = SASA(side alone) - SASA(side within the pair complex), summed over
    each residue's atoms. Uses the same Shrake-Rupley as the gate battery.
    """
    s_alone = gb.sasa(side_atoms, n_points)
    s_cplx = gb.sasa(pair_atoms, n_points)
    idx = {id(a): k for k, a in enumerate(pair_atoms)}
    out = {}
    for a, sa in zip(side_atoms, s_alone):
        sc = s_cplx[idx[id(a)]]
        d = sa - sc
        if d <= 0:
            continue
        key = (a.chain, a.resnum, a.resname)
        ap, po = out.setdefault(key, [0.0, 0.0])
        if gb.is_apolar(a):
            out[key][0] = ap + d
        else:
            out[key][1] = po + d
    return out


def featurize_structure(path, n_points, contact_cut=5.0, hb_cut=4.5):
    atoms = gb.read_pdb(path)
    chains = protein_chains(atoms)
    if len(chains) < 2:
        return []
    pair, ncon = best_pair(atoms, chains, contact_cut)
    if pair is None or ncon < 10:
        return []
    c1, c2 = pair
    a1 = [a for a in atoms if a.chain == c1 and a.resname in THREE2ONE]
    a2 = [a for a in atoms if a.chain == c2 and a.resname in THREE2ONE]
    both = a1 + a2
    rows = []
    # featurize BOTH sides: each is a legitimate "residue facing a partner"
    for side, partner in ((a1, a2), (a2, a1)):
        buried = residue_sasa_buried(side, both, n_points)
        pgrid = gb.Grid(partner, 2.0 * contact_cut)
        # group this side into residues
        res = {}
        for a in side:
            res.setdefault((a.chain, a.resnum, a.resname), []).append(a)
        for key, ratoms in res.items():
            ba, bp = buried.get(key, (0.0, 0.0))
            if ba + bp <= 0:
                continue  # not actually at the interface
            don = acc = apol = pos = neg = ncontact = 0
            seenp = set()
            for a in ratoms:
                for k in pgrid.near(a.xyz):
                    b = partner[k]
                    d2 = gb.dist2(a.xyz, b.xyz)
                    if d2 > contact_cut * contact_cut:
                        continue
                    if id(b) not in seenp:
                        seenp.add(id(b))
                        if gb.is_apolar(b):
                            apol += 1
                        if (b.resname, b.name) in CHARGED_POS:
                            pos += 1
                        if (b.resname, b.name) in CHARGED_NEG:
                            neg += 1
                        if d2 <= hb_cut * hb_cut:
                            role = partner_role(b)
                            if role in ("donor", "both"):
                                don += 1
                            if role in ("acceptor", "both"):
                                acc += 1
                    ncontact += 1
            rows.append({
                "structure": os.path.basename(path),
                "chain": key[0], "resnum": key[1],
                "resname": key[2], "aa": THREE2ONE.get(key[2], "X"),
                "buried_apolar": round(ba, 2), "buried_polar": round(bp, 2),
                "partner_apolar": apol,
                "partner_donor": don, "partner_acceptor": acc,
                "partner_pos": pos, "partner_neg": neg,
                "n_partner_contacts": ncontact,
            })
    return rows


COLS = ["structure", "chain", "resnum", "resname", "aa",
        "buried_apolar", "buried_polar", "partner_apolar",
        "partner_donor", "partner_acceptor", "partner_pos", "partner_neg",
        "n_partner_contacts"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", help="corpus raw/ directory of .pdb files")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sasa-points", type=int, default=96)
    ap.add_argument("--shard", default="1/1",
                    help="i/N to process only the i-th of N chunks "
                         "(for a Slurm array); 1-indexed")
    args = ap.parse_args(argv)
    if os.path.exists(args.out):
        raise SystemExit("%s exists; pass a fresh --out (no silent append)."
                         % args.out)

    files = sorted(glob.glob(os.path.join(args.raw, "*.pdb")))
    i, n = (int(x) for x in args.shard.split("/"))
    files = files[(i - 1)::n]
    print("featurizing %d structures (shard %s), SASA points=%d"
          % (len(files), args.shard, args.sasa_points))

    nrows = nfail = 0
    with open(args.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for k, p in enumerate(files, 1):
            try:
                rows = featurize_structure(p, args.sasa_points)
            except Exception as e:
                nfail += 1
                sys.stderr.write("FAIL %s: %s\n" % (os.path.basename(p), e))
                continue
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in COLS) + "\n")
            nrows += len(rows)
            if k % 25 == 0:
                print("  %d/%d structures, %d interface residues"
                      % (k, len(files), nrows))
    print("\nDONE  %d interface-residue rows from %d structures (%d failed)"
          % (nrows, len(files), nfail))
    print("wrote %s" % args.out)
    print("Each row = one interface residue + the partner chemistry it faces. "
          "Next: train a model to predict 'aa' from the partner_* features.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
