#!/usr/bin/env python3
"""
98_sample_backbone.py -- step 3 of the from-scratch BACKBONE generator.

Load a trained 96 checkpoint and generate backbones: integrate the SE(3) flow
from noise (t=0) to data (t=1), then rebuild idealized backbone atoms from each
residue frame and write PDB files. The generated structures feed step 4
(ProteinMPNN -> AF2 designability) and can be opened directly in ChimeraX.

Honest-geometry note
--------------------
We place N, CA, C, O with IDEAL *intra*-residue geometry in each generated
frame. The peptide bond that is NOT fixed this way is the *inter*-residue
C_i -- N_{i+1} bond: its length depends entirely on whether the model placed
consecutive frames at a sane relative pose. So this script reports the mean
C-N peptide-bond length as a cheap geometry signal -- a real backbone sits near
1.33 A; a model that hasn't learned chain geometry will be far off. That is a
pre-check, not a substitute for the designability gate in step 4.

Usage (in the enroot container, GPU optional):
  python model/98_sample_backbone.py --ckpt ~/rfd/corpus_v1/bbflow.pt \
      --n 20 --length 80 --steps 100 --out-dir ~/rfd/corpus_v1/bb_samples
  # smoke (no checkpoint, random weights -- checks the geometry/IO path):
  python model/98_sample_backbone.py --random-init --n 2 --length 60 \
      --out-dir /tmp/bb_smoke
"""

import argparse
import importlib.util
import os
import sys

import numpy as np
import torch

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "bbflow", os.path.join(_HERE, "96_flow_backbone.py"))
M = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M)

# idealized atom positions in the residue frame (columns e1=C-CA dir, e2 N-side).
# N-CA 1.458 at N-CA-C angle 111 deg; C-CA 1.525 along e1; carbonyl O 1.231,
# CA-C-O 120.8 deg, trans to N.
LOCAL = {
    "N":  np.array([-0.5225, 1.3606, 0.0], np.float32),
    "CA": np.array([0.0, 0.0, 0.0], np.float32),
    "C":  np.array([1.5250, 0.0, 0.0], np.float32),
    "O":  np.array([2.1550, -1.0580, 0.0], np.float32),
}
ORDER = ["N", "CA", "C", "O"]


def rebuild(x, R):
    """x (L,3) CA coords (angstrom), R (L,3,3) frames -> {atom: (L,3)}."""
    out = {}
    for a in ORDER:
        off = torch.tensor(LOCAL[a], device=x.device)
        out[a] = x + torch.einsum('lkm,m->lk', R, off)      # CA + R @ local
    return out


def write_pdb(atoms, path):
    L = atoms["CA"].shape[0]
    serial = 1
    with open(path, "w") as fh:
        for i in range(L):
            for a in ORDER:
                p = atoms[a][i]
                nf = (" " + a).ljust(4)
                fh.write("ATOM  %5d %s GLY A%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                         "           %s\n"
                         % (serial, nf, i + 1, p[0], p[1], p[2], a[0]))
                serial += 1
        fh.write("TER\nEND\n")


def peptide_bond_stats(atoms):
    """Mean/std of the inter-residue C_i - N_{i+1} bond length (ideal ~1.33)."""
    C = atoms["C"][:-1]; N = atoms["N"][1:]
    d = torch.linalg.norm(N - C, dim=-1)
    return d.mean().item(), d.std().item()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", default=None)
    ap.add_argument("--random-init", action="store_true",
                    help="skip the checkpoint, use random weights (IO smoke)")
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--length", type=int, default=80)
    ap.add_argument("--steps", type=int, default=100)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    torch.manual_seed(args.seed)
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    scale = M.SCALE
    prior, alpha, kappa = "iid", 1.0, 0.2
    if args.random_init:
        net = M.FrameNet(128, 4).to(dev).eval()
        print("random-init (smoke): no checkpoint loaded")
    else:
        if not args.ckpt:
            ap.error("--ckpt is required unless --random-init")
        ck = torch.load(args.ckpt, map_location=dev)
        net = M.FrameNet(ck["d_model"], ck["layers"]).to(dev).eval()
        net.load_state_dict(ck["state"])
        scale = ck.get("scale", M.SCALE)
        prior = ck.get("prior", "iid")                       # match the training prior
        alpha = ck.get("alpha", 1.0); kappa = ck.get("kappa", 0.2)
        print("loaded", args.ckpt, "| prior:", prior)

    os.makedirs(args.out_dir, exist_ok=True)
    x, R = M.sample_latent(net, args.n, args.length, args.steps, dev,
                           prior=prior, alpha=alpha, kappa=kappa)
    x = x * scale                                            # back to angstrom
    cn_all = []
    for k in range(args.n):
        atoms = rebuild(x[k], R[k])
        path = os.path.join(args.out_dir, "sample_%03d.pdb" % k)
        write_pdb({a: v.cpu().numpy() for a, v in atoms.items()}, path)
        m, s = peptide_bond_stats(atoms)
        cn_all.append(m)
    cn = float(np.mean(cn_all))
    print("wrote %d backbones (len %d) -> %s"
          % (args.n, args.length, args.out_dir))
    print("mean C-N peptide bond: %.2f A  (ideal ~1.33; far off = chain "
          "geometry not learned)" % cn)
    return 0


if __name__ == "__main__":
    sys.exit(main())
