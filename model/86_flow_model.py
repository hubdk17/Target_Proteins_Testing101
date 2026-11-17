#!/usr/bin/env python3
"""
86_flow_model.py -- our own generative FLOW-MATCHING model for de novo binder
interface sequences, conditioned on target chemistry. Trained on the blind PDB
corpus (MDM2/MDMX/p53 held out by 80/83).

What it is
----------
This is the generative version of 82/84. Those answered "does target chemistry
predict the right residue?" (yes, confirmed). This GENERATES a whole binder
interface sequence given the target pocket it must sit against -- flow matching
on the amino-acid simplex, conditioned per position on the environment-only
features from 81, with a transformer over the interface so the generated
residues are mutually COHERENT, not independent argmaxes.

It realises ranked claim 3 generatively: target-conditioned interface design
that optimises for what chemistry the pocket presents, which stage 1 of the
pipeline never sees.

Scope choice (said plainly): this generates SEQUENCE, not backbone. For MDM2 the
fold is a single ideal helix -- the backbone is nearly solved by physics, and
every bit of learnable signal we have measured lives in sequence-vs-chemistry.
A backbone SE(3) flow (our RFdiffusion) is the heavier next model; this is the
one our data actually supports today. One word pivots us to backbone/joint.

Method (conditional flow matching, x1-prediction parameterisation)
-----------------------------------------------------------------
Per interface position i with target-chemistry features c_i and true residue
one-hot x1_i on the 20-simplex:
  * draw x0_i ~ Dirichlet(1)  (uniform on the simplex)
  * draw one t ~ U(0,1) for the whole interface
  * linear interpolant   x_t = (1-t) x0 + t x1
  * net predicts clean residue logits from (x_t, c, t); loss = CE(logits, x1)
Sampling integrates dx/dt = (softmax(logits) - x_t)/(1-t) from noise to t->1.

Honesty rails (shared discipline with 82/84)
  * split BY STRUCTURE, never by row (an interface is one structure+chain)
  * class-weighted loss + BALANCED recovery (mean per-residue recall) headline,
    because 41% of interface residues are hydrophobic
  * per-chemical-group recall, so you SEE if it designs +/-/polar positions
  * default features are the columns already on disk; --enriched uses 81's new
    distance-weighted / directional columns (re-featurize first)

Run inside the enroot pytorch container on a GPU -- see 87_train_flow.sbatch.
"""

import argparse
import os
import sys

import numpy as np
import torch
import torch.nn as nn

AAS = list("ACDEFGHIKLMNPQRSTVWY")
AA_IDX = {a: i for i, a in enumerate(AAS)}
GROUP = {"A": "hyd", "V": "hyd", "L": "hyd", "I": "hyd", "M": "hyd",
         "F": "hyd", "W": "hyd", "Y": "hyd", "C": "hyd",
         "S": "pol", "T": "pol", "N": "pol", "Q": "pol", "H": "pol",
         "K": "pos", "R": "pos", "D": "neg", "E": "neg",
         "G": "sp", "P": "sp"}
GRP_OF = np.array([GROUP[a] for a in AAS])

# features already present in every interface_features.tsv
BASE_FEATS = ["partner_apolar", "partner_donor", "partner_acceptor",
              "partner_pos", "partner_neg", "n_partner_contacts",
              "buried_apolar", "buried_polar"]
# 81's enriched columns (distance-weighted + directional + local hydrophobicity)
ENRICHED_FEATS = ["w_apolar", "w_donor", "w_acceptor", "w_pos", "w_neg",
                  "cb_apolar", "cb_polar", "cb_pos", "cb_neg", "partner_kd"]


def load_groups(path, feats):
    """One group per (structure, chain) = one interface = one training example."""
    header = None
    groups = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = {c: i for i, c in enumerate(p)}
                miss = [c for c in feats + ["aa", "structure", "chain"]
                        if c not in header]
                if miss:
                    raise SystemExit(
                        "missing columns: %s -- re-featurize with the current "
                        "81_featurize_corpus.py (it adds the enriched columns)."
                        % miss)
                continue
            aa = p[header["aa"]]
            if aa not in AA_IDX:
                continue
            try:
                fv = [float(p[header[c]]) for c in feats]
            except (ValueError, IndexError):
                continue
            key = (p[header["structure"]], p[header["chain"]])
            g = groups.setdefault(key, ([], []))
            g[0].append(fv)
            g[1].append(AA_IDX[aa])
    out = []
    for (struct, _chain), (fv, y) in groups.items():
        out.append((struct, np.asarray(fv, np.float32), np.asarray(y, np.int64)))
    return out


def batch_pad(groups, idx, dev):
    """Pad a list of interfaces to the longest in the batch; return mask."""
    sub = [groups[i] for i in idx]
    L = max(g[1].shape[0] for g in sub)
    F = sub[0][1].shape[1]
    B = len(sub)
    C = np.zeros((B, L, F), np.float32)
    Y = np.zeros((B, L), np.int64)
    M = np.zeros((B, L), bool)
    for k, (_s, fv, y) in enumerate(sub):
        n = y.shape[0]
        C[k, :n] = fv
        Y[k, :n] = y
        M[k, :n] = True
    return (torch.tensor(C, device=dev), torch.tensor(Y, device=dev),
            torch.tensor(M, device=dev))


def dirichlet_uniform(shape, dev):
    """x0 ~ Dirichlet(1) = uniform on the simplex (via normalised Exp(1))."""
    g = -torch.log(torch.rand(shape, device=dev).clamp_min(1e-9))
    return g / g.sum(-1, keepdim=True)


class FlowNet(nn.Module):
    def __init__(self, d_feat, n_aa=20, d=128, layers=3, heads=4, drop=0.1):
        super().__init__()
        self.inp = nn.Linear(n_aa + d_feat + 1, d)
        enc = nn.TransformerEncoderLayer(d, heads, 4 * d, drop,
                                         batch_first=True)
        self.tr = nn.TransformerEncoder(enc, layers)
        self.out = nn.Linear(d, n_aa)

    def forward(self, xt, c, t, mask):
        # xt: B,L,20 | c: B,L,F | t: B,1,1 | mask: B,L (True=valid)
        tt = t.expand(xt.size(0), xt.size(1), 1)
        h = self.inp(torch.cat([xt, c, tt], -1))
        h = self.tr(h, src_key_padding_mask=~mask)  # transformer: True = pad
        return self.out(h)


def train(net, groups, tr_idx, mu, sd, args, dev, class_w):
    opt = torch.optim.Adam(net.parameters(), lr=args.lr,
                           weight_decay=args.l2)
    lossf = nn.CrossEntropyLoss(weight=class_w)
    rng = np.random.default_rng(args.seed)
    n_aa = len(AAS)
    for ep in range(args.epochs):
        net.train()
        order = rng.permutation(tr_idx)
        tot = nb = 0.0
        for i in range(0, len(order), args.batch):
            idx = order[i:i + args.batch]
            C, Y, M = batch_pad(groups, idx, dev)
            C = (C - mu) / sd
            x1 = torch.eye(n_aa, device=dev)[Y]          # B,L,20
            x0 = dirichlet_uniform(x1.shape, dev)
            t = torch.rand(C.size(0), 1, 1, device=dev)
            xt = (1 - t) * x0 + t * x1
            logits = net(xt, C, t, M)
            loss = lossf(logits[M], Y[M])
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss)
            nb += 1
        if ep == 0 or (ep + 1) % max(1, args.epochs // 10) == 0:
            print("  epoch %4d/%d  train loss %.4f"
                  % (ep + 1, args.epochs, tot / max(nb, 1)))


@torch.no_grad()
def sample_recovery(net, groups, te_idx, mu, sd, args, dev, steps=50):
    """Integrate the flow from noise and compare designed vs native residue."""
    net.eval()
    n_aa = len(AAS)
    all_true, all_pred = [], []
    for i in range(0, len(te_idx), args.batch):
        idx = te_idx[i:i + args.batch]
        C, Y, M = batch_pad(groups, idx, dev)
        Cn = (C - mu) / sd
        xt = dirichlet_uniform((C.size(0), C.size(1), n_aa), dev)
        for s in range(steps):
            t = torch.full((C.size(0), 1, 1), s / steps, device=dev)
            logits = net(xt, Cn, t, M)
            p1 = torch.softmax(logits, -1)
            dt = 1.0 / steps
            xt = xt + dt * (p1 - xt) / max(1e-3, 1 - s / steps)
        pred = xt.argmax(-1)
        all_true.append(Y[M].cpu().numpy())
        all_pred.append(pred[M].cpu().numpy())
    return np.concatenate(all_true), np.concatenate(all_pred)


def report(y_true, y_pred, freq_counts):
    n_aa = len(AAS)
    rec = (y_pred == y_true).mean()
    freq = freq_counts / freq_counts.sum()
    maj = freq_counts.argmax()
    base = (y_true == maj).mean()
    recs = []
    for c in range(n_aa):
        m = y_true == c
        if m.any():
            recs.append((y_pred[m] == c).mean())
    bacc = float(np.mean(recs))
    print("\n================ FLOW MODEL -- held-out recovery ================")
    print("  sequence recovery : model %.3f  | baseline(always %s) %.3f"
          % (rec, AAS[maj], base))
    print("  BALANCED recovery : model %.3f  | chance %.3f   "
          "<-- imbalance-robust headline" % (bacc, 1.0 / n_aa))
    gt, gp = GRP_OF[y_true], GRP_OF[y_pred]
    print("  per-group recall (did it DESIGN +/-/polar positions right?):")
    for g in ["hyd", "pol", "pos", "neg", "sp"]:
        m = gt == g
        if m.any():
            print("     %-3s  recall %.3f   (n=%d)" % (g, (gp[m] == g).mean(),
                                                       m.sum()))
    print("  --> compare BALANCED recovery to 84 (MLP) / 82 (linear) on the "
          "same corpus: a generator that couples positions should match or beat "
          "the per-position classifier, and give DIVERSE valid designs.")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--enriched", action="store_true",
                    help="use 81's distance-weighted/directional columns too")
    ap.add_argument("--d-model", type=int, default=128)
    ap.add_argument("--layers", type=int, default=3)
    ap.add_argument("--heads", type=int, default=4)
    ap.add_argument("--dropout", type=float, default=0.1)
    ap.add_argument("--epochs", type=int, default=120)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--lr", type=float, default=5e-4)
    ap.add_argument("--l2", type=float, default=1e-5)
    ap.add_argument("--test-frac", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--save", default=None)
    ap.add_argument("--cpu", action="store_true")
    args = ap.parse_args(argv)

    feats = BASE_FEATS + (ENRICHED_FEATS if args.enriched else [])
    dev = torch.device("cpu" if args.cpu or not torch.cuda.is_available()
                       else "cuda")
    print("torch %s  device %s  features(%d): %s"
          % (torch.__version__, dev, len(feats), ", ".join(feats)))

    groups = load_groups(args.table, feats)
    print("loaded %d interfaces (one per structure+chain) from %d structures"
          % (len(groups), len(set(g[0] for g in groups))))

    # split BY STRUCTURE
    structs = sorted(set(g[0] for g in groups))
    rng = np.random.default_rng(args.seed)
    rng.shuffle(structs)
    n_test = max(1, int(len(structs) * args.test_frac))
    test_structs = set(structs[:n_test])
    tr_idx = np.array([i for i, g in enumerate(groups)
                       if g[0] not in test_structs])
    te_idx = np.array([i for i, g in enumerate(groups)
                       if g[0] in test_structs])
    print("split BY STRUCTURE: %d train / %d test interfaces"
          % (len(tr_idx), len(te_idx)))

    # standardise features on TRAIN positions only
    trC = np.concatenate([groups[i][1] for i in tr_idx], 0)
    mu = torch.tensor(trC.mean(0), device=dev)
    sd = torch.tensor(trC.std(0), device=dev)
    sd[sd == 0] = 1.0

    # class weights + frequency baseline from TRAIN labels
    tr_y = np.concatenate([groups[i][2] for i in tr_idx])
    counts = np.bincount(tr_y, minlength=len(AAS)).astype(np.float32)
    w = counts.sum() / (len(AAS) * np.clip(counts, 1, None))
    class_w = torch.tensor(w, device=dev)

    torch.manual_seed(args.seed)
    net = FlowNet(len(feats), len(AAS), args.d_model, args.layers,
                  args.heads, args.dropout).to(dev)
    nparam = sum(p.numel() for p in net.parameters())
    print("FlowNet: d=%d layers=%d heads=%d  params=%d"
          % (args.d_model, args.layers, args.heads, nparam))

    train(net, groups, tr_idx, mu, sd, args, dev, class_w)
    y_true, y_pred = sample_recovery(net, groups, te_idx, mu, sd, args, dev)
    report(y_true, y_pred, counts)

    if args.save:
        torch.save({"state": net.state_dict(), "feats": feats,
                    "mu": mu.cpu().numpy(), "sd": sd.cpu().numpy(),
                    "aas": AAS, "cfg": vars(args)}, args.save)
        print("\nsaved -> %s" % args.save)
    print("\nThis is a generator: it SAMPLES interface sequences for a given "
          "target pocket. Recovery on held-out structures is the blind check; "
          "lab validation comes later, as planned.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
