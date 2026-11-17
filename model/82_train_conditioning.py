#!/usr/bin/env python3
"""
82_train_conditioning.py -- does the TARGET chemistry at an interface position
predict which residue belongs there? (ranked claim 3)

Learns  P( residue | chemistry the position faces on the partner )  from
81_featurize_corpus.py's table. It is NOT a generator -- it is the cheap test of
whether the signal a bigger flow model would depend on is even there.

v2 note (2026-11-17): the first version's headline was "chemical-group accuracy"
vs an always-guess-the-commonest-group baseline. That baseline is ~0.41 because
41% of interface residues are hydrophobic, so a model that correctly diversifies
to catch charged/polar positions can LOSE on top-1 group while clearly winning on
exact accuracy. That metric was fooled by class imbalance. This version reports
imbalance-robust metrics instead:
  * exact top-1 and top-3 accuracy (vs frequency baselines)
  * BALANCED accuracy = mean per-residue recall (majority-guess scores ~1/20)
  * per-chemical-group recall, so you SEE whether it catches +/-/polar positions
and with --compare it runs partner-chemistry-only AND +burial side by side.

Honesty rails: split BY STRUCTURE (never by row); weights printed so you can read
WHAT it learned; burial OFF by default (it partly encodes the answer).

Run inside the enroot container (login python3 has no numpy):
  srun --cpus-per-task=4 --mem=16G \
    --container-image=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh \
    --container-mounts=/home/dakshcsed:/workspace --container-workdir=/workspace \
    python /workspace/rfd/scripts/model/82_train_conditioning.py \
      /workspace/rfd/corpus_v1/interface_features.tsv --compare
"""

import argparse
import sys
import numpy as np

AAS = list("ACDEFGHIKLMNPQRSTVWY")
AA_IDX = {a: i for i, a in enumerate(AAS)}
GROUP = {"A": "hyd", "V": "hyd", "L": "hyd", "I": "hyd", "M": "hyd",
         "F": "hyd", "W": "hyd", "Y": "hyd", "C": "hyd",
         "S": "pol", "T": "pol", "N": "pol", "Q": "pol", "H": "pol",
         "K": "pos", "R": "pos", "D": "neg", "E": "neg",
         "G": "sp", "P": "sp"}
GRP_OF = np.array([GROUP[a] for a in AAS])
PARTNER_FEATS = ["partner_apolar", "partner_donor", "partner_acceptor",
                 "partner_pos", "partner_neg", "n_partner_contacts"]
BURIAL_FEATS = ["buried_apolar", "buried_polar"]
ALLCOLS = PARTNER_FEATS + BURIAL_FEATS


def load(path):
    header = None
    cols = {c: [] for c in ALLCOLS}
    y, struct = [], []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = {c: i for i, c in enumerate(p)}
                miss = [c for c in ALLCOLS + ["aa", "structure"]
                        if c not in header]
                if miss:
                    raise SystemExit("missing columns: %s" % miss)
                continue
            aa = p[header["aa"]]
            if aa not in AA_IDX:
                continue
            try:
                vals = {c: float(p[header[c]]) for c in ALLCOLS}
            except (ValueError, IndexError):
                continue
            for c in ALLCOLS:
                cols[c].append(vals[c])
            y.append(AA_IDX[aa])
            struct.append(p[header["structure"]])
    X = {c: np.array(v, float) for c, v in cols.items()}
    return X, np.array(y, int), np.array(struct)


def softmax(z):
    z = z - z.max(1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(1, keepdims=True)


def train_softmax(X, y, n_class, epochs=400, lr=0.5, l2=1e-4):
    n, d = X.shape
    W = np.zeros((d, n_class))
    b = np.zeros(n_class)
    Y = np.eye(n_class)[y]
    for _ in range(epochs):
        P = softmax(X @ W + b)
        W -= lr * (X.T @ (P - Y) / n + l2 * W)
        b -= lr * (P - Y).mean(0)
    return W, b


def balanced_acc(y_true, y_pred, n_class):
    recs = []
    for c in range(n_class):
        m = y_true == c
        if m.any():
            recs.append((y_pred[m] == c).mean())
    return float(np.mean(recs))


def run(Xall, y, tr, te, feats, label, epochs, save=None):
    X = np.column_stack([Xall[c] for c in feats])
    mu, sd = X[tr].mean(0), X[tr].std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    W, b = train_softmax(Xs[tr], y[tr], len(AAS), epochs=epochs)
    P = softmax(Xs[te] @ W + b)
    yte = y[te]
    pred = P.argmax(1)
    top3 = np.argsort(P, 1)[:, -3:]

    counts = np.bincount(y[tr], minlength=len(AAS))
    freq = counts / counts.sum()
    maj = counts.argmax()
    top3_freq = np.argsort(freq)[-3:]

    acc = (pred == yte).mean()
    acc_maj = (yte == maj).mean()
    acc3 = np.mean([yt in row for yt, row in zip(yte, top3)])
    acc3_base = sum(freq[k] for k in top3_freq)  # approx on test ~ train freq
    bacc = balanced_acc(yte, pred, len(AAS))
    bacc_base = 1.0 / len(AAS)

    print("\n================ %s ================" % label)
    print("features: %s" % ", ".join(feats))
    print("  exact top-1 : model %.3f  | baseline(always %s) %.3f"
          % (acc, AAS[maj], acc_maj))
    print("  exact top-3 : model %.3f  | baseline(freq top-3)  %.3f"
          % (acc3, acc3_base))
    print("  BALANCED acc: model %.3f  | baseline(majority)    %.3f   "
          "<-- imbalance-robust headline" % (bacc, bacc_base))
    # per-group recall
    gt = GRP_OF[yte]
    gp = GRP_OF[pred]
    print("  per-group recall (did it call +/-/polar positions right?):")
    for g in ["hyd", "pol", "pos", "neg", "sp"]:
        m = gt == g
        if m.any():
            print("     %-3s  recall %.3f   (n=%d)" % (g, (gp[m] == g).mean(),
                                                       m.sum()))
    verdict = "SIGNAL" if (acc > acc_maj and bacc > 2 * bacc_base) else "weak/none"
    print("  --> %s (beats exact baseline AND balanced acc >> chance)" % verdict)

    print("  weights (target feature -> residues it pushes toward | away):")
    for fi, fname in enumerate(feats):
        order = np.argsort(W[fi])
        up = " ".join("%s(+%.2f)" % (AAS[k], W[fi, k]) for k in order[::-1][:4])
        dn = " ".join("%s(%.2f)" % (AAS[k], W[fi, k]) for k in order[:2])
        print("     %-20s toward: %-28s | away: %s" % (fname, up, dn))
    if save:
        np.savez(save, W=W, b=b, mu=mu, sd=sd, feats=feats, aas=AAS)
        print("  saved -> %s" % save)
    return acc, bacc


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--with-burial", action="store_true")
    ap.add_argument("--compare", action="store_true",
                    help="train partner-only AND +burial, side by side")
    ap.add_argument("--test-frac", type=float, default=0.2)
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--save", default=None)
    args = ap.parse_args(argv)

    Xall, y, struct = load(args.table)
    print("loaded %d interface residues from %d structures"
          % (len(y), len(set(struct))))

    uniq = np.array(sorted(set(struct)))
    rng = np.random.default_rng(args.seed)
    rng.shuffle(uniq)
    n_test = max(1, int(len(uniq) * args.test_frac))
    test_set = set(uniq[:n_test])
    te = np.array([s in test_set for s in struct])
    tr = ~te
    print("split BY STRUCTURE: %d train / %d test structures; %d / %d rows"
          % (len(uniq) - n_test, n_test, tr.sum(), te.sum()))

    if args.compare:
        run(Xall, y, tr, te, PARTNER_FEATS, "PARTNER CHEMISTRY ONLY",
            args.epochs)
        run(Xall, y, tr, te, PARTNER_FEATS + BURIAL_FEATS,
            "PARTNER CHEMISTRY + BURIAL", args.epochs, save=args.save)
    else:
        feats = PARTNER_FEATS + (BURIAL_FEATS if args.with_burial else [])
        run(Xall, y, tr, te, feats, "MODEL", args.epochs, save=args.save)

    print("\nRead: BALANCED accuracy >> 0.05 (chance) and > baseline means the "
          "target chemistry genuinely predicts residue identity, class imbalance "
          "aside. Compare partner-only vs +burial to see how much geometry adds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
