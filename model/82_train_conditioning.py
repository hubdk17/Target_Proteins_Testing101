#!/usr/bin/env python3
"""
82_train_conditioning.py -- the first learned model: does the TARGET chemistry
at an interface position predict which residue belongs there?

This is the honest first rung toward a generative binder model (ranked claim 3:
target-conditioned interface chemistry). It is NOT a generator. It learns

        P( residue identity | the chemistry the position faces on the partner )

from 81_featurize_corpus.py's table, so we find out whether the signal a bigger
flow model would depend on is even there -- before spending weeks (and an A100)
building the flow model. If this cannot beat the amino-acid-frequency baseline,
no 24M-parameter model will.

Honesty rails
-------------
* Split is BY STRUCTURE, never by row: every row from one complex goes wholly to
  train or wholly to test, so a near-duplicate interface cannot sit on both sides
  and inflate the score. (The corpus is not yet sequence-deduped -- this split is
  the cheap guard until it is.)
* Features are the TARGET-side chemistry the position faces: partner apolar /
  donor / acceptor / +charge / -charge counts and contact count. Burial
  (buried_apolar/polar) is OFF by default because it partly encodes the answer
  (hydrophobics bury more); turn it on with --with-burial to see how much it adds.
* Two baselines, both reported: always-guess-the-commonest-residue, and
  guess-weighted-by-training-frequency. The model has to beat both.
* Exact 20-way residue accuracy is hard and not the claim; the claim is that
  chemistry predicts residue CHEMISTRY, so the chemical-GROUP accuracy
  (hydrophobic / polar / +/- / special) is the headline number.
* The learned weights are printed, so you can see WHAT it learned (apolar->big
  hydrophobics, -charge partner -> Lys/Arg, etc.) rather than trusting a number.

Usage (CPU, ~minutes; numpy only -- run in the pytorch container to be safe):
  python 82_train_conditioning.py ~/rfd/corpus_v1/interface_features.tsv
"""

import argparse
import sys
import numpy as np

AAS = list("ACDEFGHIKLMNPQRSTVWY")
AA_IDX = {a: i for i, a in enumerate(AAS)}
GROUP = {  # chemical group of each residue
    "A": "hydrophobic", "V": "hydrophobic", "L": "hydrophobic",
    "I": "hydrophobic", "M": "hydrophobic", "F": "hydrophobic",
    "W": "hydrophobic", "Y": "hydrophobic", "C": "hydrophobic",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar", "H": "polar",
    "K": "pos", "R": "pos", "D": "neg", "E": "neg",
    "G": "special", "P": "special",
}
PARTNER_FEATS = ["partner_apolar", "partner_donor", "partner_acceptor",
                 "partner_pos", "partner_neg", "n_partner_contacts"]
BURIAL_FEATS = ["buried_apolar", "buried_polar"]


def load(path, feats):
    header = None
    X, y, struct = [], [], []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = {c: i for i, c in enumerate(p)}
                need = feats + ["aa", "structure"]
                miss = [c for c in need if c not in header]
                if miss:
                    raise SystemExit("missing columns: %s" % miss)
                continue
            aa = p[header["aa"]]
            if aa not in AA_IDX:
                continue
            try:
                X.append([float(p[header[c]]) for c in feats])
            except (ValueError, IndexError):
                continue
            y.append(AA_IDX[aa])
            struct.append(p[header["structure"]])
    return np.array(X, float), np.array(y, int), np.array(struct)


def softmax(z):
    z = z - z.max(1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(1, keepdims=True)


def train_softmax(X, y, n_class, epochs=400, lr=0.5, l2=1e-4, seed=0):
    rng = np.random.default_rng(seed)
    n, d = X.shape
    W = np.zeros((d, n_class))
    b = np.zeros(n_class)
    Y = np.eye(n_class)[y]
    for _ in range(epochs):
        P = softmax(X @ W + b)
        gW = X.T @ (P - Y) / n + l2 * W
        gb = (P - Y).mean(0)
        W -= lr * gW
        b -= lr * gb
    return W, b


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--with-burial", action="store_true")
    ap.add_argument("--test-frac", type=float, default=0.2)
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--save", default=None, help="save weights to .npz")
    args = ap.parse_args(argv)

    feats = PARTNER_FEATS + (BURIAL_FEATS if args.with_burial else [])
    X, y, struct = load(args.table, feats)
    print("loaded %d interface residues, %d features (%s)"
          % (len(y), len(feats), ", ".join(feats)))
    if len(y) < 1000:
        print("WARNING: very few rows; results will be noisy.")

    # split BY STRUCTURE
    uniq = np.array(sorted(set(struct)))
    rng = np.random.default_rng(args.seed)
    rng.shuffle(uniq)
    n_test = max(1, int(len(uniq) * args.test_frac))
    test_set = set(uniq[:n_test])
    te = np.array([s in test_set for s in struct])
    tr = ~te
    print("structures: %d train / %d test;  rows: %d train / %d test"
          % (len(uniq) - n_test, n_test, tr.sum(), te.sum()))

    # standardise on train stats only
    mu, sd = X[tr].mean(0), X[tr].std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd

    W, b = train_softmax(Xs[tr], y[tr], len(AAS),
                         epochs=args.epochs, seed=args.seed)
    P = softmax(Xs[te] @ W + b)
    pred = P.argmax(1)
    yte = y[te]

    # baselines (from TRAIN)
    counts = np.bincount(y[tr], minlength=len(AAS))
    freq = counts / counts.sum()
    maj = counts.argmax()
    acc_model = (pred == yte).mean()
    acc_maj = (yte == maj).mean()
    # expected accuracy of frequency-weighted guessing = sum f_i^2 (on test dist)
    te_counts = np.bincount(yte, minlength=len(AAS))
    te_freq = te_counts / te_counts.sum()
    acc_weighted = float((freq * te_freq).sum())  # guess~train freq, eval~test

    # chemical-group accuracy (the headline)
    grp = np.array([list(GROUP)[0]] * 0)  # placeholder
    g_pred = np.array([GROUP[AAS[i]] for i in pred])
    g_true = np.array([GROUP[AAS[i]] for i in yte])
    acc_group = (g_pred == g_true).mean()
    # group baseline: most common group in train
    from collections import Counter
    gtr = Counter(GROUP[AAS[i]] for i in y[tr])
    maj_group = gtr.most_common(1)[0][0]
    acc_group_base = (g_true == maj_group).mean()

    print("\n================ RESULTS (held-out structures) ================")
    print("exact 20-way residue accuracy:")
    print("  model                  : %.3f" % acc_model)
    print("  baseline most-common   : %.3f  (always '%s')" % (acc_maj, AAS[maj]))
    print("  baseline freq-weighted : %.3f" % acc_weighted)
    print("chemical-GROUP accuracy (hydrophobic/polar/+/-/special) -- headline:")
    print("  model                  : %.3f" % acc_group)
    print("  baseline most-common   : %.3f  (always '%s')"
          % (acc_group_base, maj_group))
    lift = acc_group - acc_group_base
    print("  --> group lift over baseline: %+.3f %s"
          % (lift, "(signal present)" if lift > 0.02 else "(NO signal)"))

    # what did it learn? top residues pushed by each partner feature
    print("\n================ WHAT IT LEARNED (weights) ================")
    print("(for each target-chemistry feature, the residues it most pushes "
          "toward / away)")
    for fi, fname in enumerate(feats):
        order = np.argsort(W[fi])
        up = " ".join("%s(+%.2f)" % (AAS[k], W[fi, k]) for k in order[::-1][:4])
        dn = " ".join("%s(%.2f)" % (AAS[k], W[fi, k]) for k in order[:3])
        print("  %-20s  toward: %s   | away: %s" % (fname, up, dn))

    if args.save:
        np.savez(args.save, W=W, b=b, mu=mu, sd=sd, feats=feats, aas=AAS)
        print("\nsaved model to %s" % args.save)
    print("\nRead: group lift > 0 means target chemistry DOES predict residue "
          "chemistry -- the signal a generative model would need is real.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
