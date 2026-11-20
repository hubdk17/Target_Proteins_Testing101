#!/usr/bin/env python3
"""
84_train_mlp.py -- our first real neural model: a small MLP that learns
P( residue | the target chemistry the position faces ), trained on GPU.

This is the next rung above 82_train_conditioning.py (a linear softmax). Same
blind corpus, same honest metrics, same by-structure split -- the ONLY thing
that changes is model capacity. That makes this a clean, falsifiable experiment:

  does a nonlinear model pull MORE of the target-chemistry signal out of the
  same features than the linear baseline?

  * if YES  -> capacity helps; scaling toward the generative model is justified.
  * if NO   -> the FEATURES are the bottleneck, not the model, and the right
               move is richer features (enrich 81), not a bigger net.

Either answer is useful and neither can be faked, because the test set is held
out BY STRUCTURE and the baselines are printed next to every number.

This is the sequence head of the eventual joint generator, trained standalone.
It is NOT a generator yet.

Honesty rails (shared with 82 -- loaded from it, one source of truth):
  * split BY STRUCTURE, never by row
  * balanced accuracy (mean per-residue recall, chance ~0.05) is the headline,
    because 41% of interface residues are hydrophobic and top-1 is imbalance-fooled
  * per-chemical-group recall, so you SEE if it catches +/-/polar positions
  * burial OFF by default (it partly encodes the answer); --compare runs both

Run inside the enroot container on a GPU (see 85_train_mlp.sbatch). Falls back
to CPU with --cpu; the dataset is small so CPU also finishes in minutes.
"""

import argparse
import importlib.util
import os
import sys

import numpy as np

# Reuse 82's loader, label maps and metrics verbatim so the two models are
# strictly comparable and there is a single definition of "interface residue".
_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "train_conditioning", os.path.join(_HERE, "82_train_conditioning.py"))
tc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tc)

AAS = tc.AAS
GRP_OF = tc.GRP_OF
PARTNER_FEATS = tc.PARTNER_FEATS
BURIAL_FEATS = tc.BURIAL_FEATS
balanced_acc = tc.balanced_acc


def build_mlp(torch, nn, d_in, n_class, hidden, dropout):
    layers = []
    prev = d_in
    for h in hidden:
        layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
        prev = h
    layers += [nn.Linear(prev, n_class)]
    return nn.Sequential(*layers)


def run(torch, Xall, y, tr, te, feats, label, args, save=None):
    nn = torch.nn
    dev = torch.device("cpu" if args.cpu or not torch.cuda.is_available()
                       else "cuda")

    X = np.column_stack([Xall[c] for c in feats]).astype(np.float32)
    mu, sd = X[tr].mean(0), X[tr].std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd

    Xtr = torch.tensor(Xs[tr], device=dev)
    ytr = torch.tensor(y[tr], device=dev, dtype=torch.long)
    Xte = torch.tensor(Xs[te], device=dev)
    yte_np = y[te]

    # Class-weighted loss so the net is pushed to call the rare charged/polar
    # residues, not just park on hydrophobic. This is why we read BALANCED acc.
    counts = np.bincount(y[tr], minlength=len(AAS)).astype(np.float32)
    w = counts.sum() / (len(AAS) * np.clip(counts, 1, None))
    class_w = torch.tensor(w, device=dev)

    torch.manual_seed(args.seed)
    net = build_mlp(torch, nn, X.shape[1], len(AAS),
                    [int(h) for h in args.hidden.split(",") if h], args.dropout)
    net.to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=args.lr, weight_decay=args.l2)
    lossf = nn.CrossEntropyLoss(weight=class_w)

    n = Xtr.shape[0]
    for ep in range(args.epochs):
        net.train()
        perm = torch.randperm(n, device=dev)
        for i in range(0, n, args.batch):
            idx = perm[i:i + args.batch]
            opt.zero_grad()
            loss = lossf(net(Xtr[idx]), ytr[idx])
            loss.backward()
            opt.step()
        if ep == 0 or (ep + 1) % max(1, args.epochs // 5) == 0:
            print("  epoch %4d/%d  train loss %.4f"
                  % (ep + 1, args.epochs, float(loss)))

    net.eval()
    with torch.no_grad():
        logits = net(Xte)
        P = torch.softmax(logits, 1).cpu().numpy()
    pred = P.argmax(1)
    top3 = np.argsort(P, 1)[:, -3:]

    cnt = np.bincount(y[tr], minlength=len(AAS))
    freq = cnt / cnt.sum()
    maj = cnt.argmax()
    top3_freq = np.argsort(freq)[-3:]

    acc = (pred == yte_np).mean()
    acc_maj = (yte_np == maj).mean()
    acc3 = np.mean([yt in row for yt, row in zip(yte_np, top3)])
    acc3_base = sum(freq[k] for k in top3_freq)
    bacc = balanced_acc(yte_np, pred, len(AAS))

    print("\n================ %s (MLP %s) ================"
          % (label, args.hidden))
    print("features: %s   device: %s" % (", ".join(feats), dev))
    print("  exact top-1 : model %.3f  | baseline(always %s) %.3f"
          % (acc, AAS[maj], acc_maj))
    print("  exact top-3 : model %.3f  | baseline(freq top-3)  %.3f"
          % (acc3, acc3_base))
    print("  BALANCED acc: model %.3f  | chance %.3f   "
          "<-- imbalance-robust headline" % (bacc, 1.0 / len(AAS)))
    gt = GRP_OF[yte_np]
    gp = GRP_OF[pred]
    print("  per-group recall (did it call +/-/polar positions right?):")
    for g in ["hyd", "pol", "pos", "neg", "sp"]:
        m = gt == g
        if m.any():
            print("     %-3s  recall %.3f   (n=%d)" % (g, (gp[m] == g).mean(),
                                                       m.sum()))
    verdict = ("SIGNAL" if (acc > acc_maj and bacc > 2.0 / len(AAS))
               else "weak/none")
    print("  --> %s. Compare BALANCED acc to 82 (linear) on the SAME split: "
          "higher => capacity helps; same => features are the bottleneck."
          % verdict)
    if save:
        np.savez(save, mu=mu, sd=sd, feats=feats, aas=AAS,
                 hidden=args.hidden)
        torch.save(net.state_dict(), save + ".pt")
        print("  saved -> %s(.npz/.pt)" % save)
    return acc, bacc


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--hidden", default="128,64",
                    help="comma-separated hidden layer sizes")
    ap.add_argument("--dropout", type=float, default=0.2)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--batch", type=int, default=4096)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--l2", type=float, default=1e-5)
    ap.add_argument("--test-frac", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--with-burial", action="store_true")
    ap.add_argument("--compare", action="store_true",
                    help="train partner-only AND +burial, side by side")
    ap.add_argument("--save", default=None)
    ap.add_argument("--cpu", action="store_true")
    args = ap.parse_args(argv)

    import torch  # imported here so --help works without torch installed
    print("torch %s  cuda available: %s" % (torch.__version__,
                                            torch.cuda.is_available()))

    Xall, y, struct = tc.load(args.table)
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
        run(torch, Xall, y, tr, te, PARTNER_FEATS,
            "PARTNER CHEMISTRY ONLY", args)
        run(torch, Xall, y, tr, te, PARTNER_FEATS + BURIAL_FEATS,
            "PARTNER CHEMISTRY + BURIAL", args, save=args.save)
    else:
        feats = PARTNER_FEATS + (BURIAL_FEATS if args.with_burial else [])
        run(torch, Xall, y, tr, te, feats, "MODEL", args, save=args.save)

    print("\nRead: this is a capacity test on fixed features. If BALANCED acc "
          "beats 82's linear number on the same split, a bigger model is worth "
          "it; if not, enrich the features (81) before scaling.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
