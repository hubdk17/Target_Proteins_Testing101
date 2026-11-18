#!/usr/bin/env python3
"""
99_designability.py -- step 4 scoring core: the self-consistency / designability
metric for generated backbones.

Designability pipeline (field standard, e.g. FrameDiff/RFdiffusion papers):
  generated backbone --ProteinMPNN--> N sequences --AF2--> N predicted structures
  --> scRMSD = min over the N of CA-RMSD(prediction, generated backbone).
A backbone is DESIGNABLE if scRMSD < 2.0 A: some sequence folds back into the
shape we drew. The fraction of designable backbones is the headline, alongside
diversity (recovery alone undersells a generator).

This file is the SCORING step only -- pure geometry (Kabsch superposition), so
it is testable without a GPU. The ProteinMPNN and AF2 runs reuse the already-
proven cluster scripts (33_mpnn_multi pattern on venv/ProteinMPNN; 40_gate2
predict.py -force_monomer on venv_af2); they are wired in 100_designability.sbatch
once a checkpoint exists to test against.

Usage:
  # one backbone vs its AF2 predictions:
  python model/99_designability.py --ref gen/sample_000.pdb --preds af2_out/sample_000
  # a whole set (one subdir of predictions per backbone):
  python model/99_designability.py --ref-dir gen --preds-dir af2_out --cutoff 2.0
"""

import argparse
import glob
import os
import sys

import numpy as np


def ca_coords(path):
    """CA coordinates in residue order from a PDB (first model, one chain ok)."""
    xs = []
    for ln in open(path):
        if ln.startswith("ENDMDL"):
            break
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA" and ln[16] in " A":
            try:
                xs.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            except ValueError:
                pass
    return np.asarray(xs, dtype=np.float64)


def kabsch_rmsd(P, Q):
    """Minimal CA-RMSD between P and Q (L,3) over optimal rigid superposition."""
    if P.shape != Q.shape or len(P) == 0:
        return float("nan")
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    diff = (R @ Pc.T).T - Qc
    return float(np.sqrt((diff ** 2).sum() / len(P)))


def score_one(ref_pdb, preds_dir):
    ref = ca_coords(ref_pdb)
    best = float("inf"); bestf = None
    preds = sorted(glob.glob(os.path.join(preds_dir, "*.pdb")))
    for p in preds:
        q = ca_coords(p)
        if q.shape == ref.shape:
            r = kabsch_rmsd(q, ref)
            if r < best:
                best, bestf = r, os.path.basename(p)
    return best, bestf, len(preds)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref")
    ap.add_argument("--preds")
    ap.add_argument("--ref-dir")
    ap.add_argument("--preds-dir")
    ap.add_argument("--cutoff", type=float, default=2.0)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest()

    if args.ref and args.preds:
        best, bf, n = score_one(args.ref, args.preds)
        print("scRMSD %.3f A  (best of %d: %s)  -> %s"
              % (best, n, bf, "DESIGNABLE" if best < args.cutoff else "no"))
        return 0

    if args.ref_dir and args.preds_dir:
        refs = sorted(glob.glob(os.path.join(args.ref_dir, "*.pdb")))
        n_ok = n_tot = 0
        scs = []
        for ref in refs:
            stem = os.path.splitext(os.path.basename(ref))[0]
            pdir = os.path.join(args.preds_dir, stem)
            if not os.path.isdir(pdir):
                continue
            best, bf, n = score_one(ref, pdir)
            n_tot += 1
            scs.append(best)
            ok = best < args.cutoff
            n_ok += ok
            print("%-20s scRMSD %.3f  %s" % (stem, best, "DESIGNABLE" if ok else ""))
        if n_tot:
            scs = np.array(scs)
            print("\n=== designability ===")
            print("backbones scored : %d" % n_tot)
            print("designable (<%.1f): %d  (%.1f%%)"
                  % (args.cutoff, n_ok, 100.0 * n_ok / n_tot))
            print("scRMSD median    : %.3f A" % np.median(scs))
        return 0

    ap.error("give --ref/--preds or --ref-dir/--preds-dir (or --selftest)")


def selftest():
    """Kabsch must be invariant to rigid motion and recover a known noise level."""
    rng = np.random.default_rng(0)
    P = rng.normal(size=(80, 3)) * 10
    # exact rigid copy -> ~0
    A = rng.normal(size=(3, 3)); Qr, _ = np.linalg.qr(A)
    if np.linalg.det(Qr) < 0:
        Qr[:, 0] = -Qr[:, 0]
    Pr = (Qr @ P.T).T + rng.normal(size=3) * 5
    r0 = kabsch_rmsd(P, Pr)
    # add known gaussian noise sigma -> rmsd ~ sqrt(3)*sigma for small sigma
    sigma = 0.5
    Pn = (Qr @ (P + rng.normal(size=P.shape) * sigma).T).T
    rn = kabsch_rmsd(P, Pn)
    print("rigid-motion rmsd : %.2e  (expect ~0)" % r0)
    print("noise sigma=%.2f -> rmsd %.3f  (expect ~%.3f)"
          % (sigma, rn, np.sqrt(3) * sigma))
    ok = r0 < 1e-6 and abs(rn - np.sqrt(3) * sigma) < 0.2
    print("selftest:", "PASS" if ok else "FAIL")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
