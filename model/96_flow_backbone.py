#!/usr/bin/env python3
"""
96_flow_backbone.py -- step 2 of the from-scratch BACKBONE generator.

An SE(3) flow-matching model over backbone frames (FrameDiff / FrameFlow
recipe, implemented from scratch). Each residue is a rigid frame (R in SO(3),
t in R^3) built from its N, CA, C atoms. Training flows random frames into
protein frames; sampling (step 3) runs the flow backward from noise.

The flow (per residue, conditioned on the whole chain)
------------------------------------------------------
  translations : straight-line / OT path  x_t = (1-t) x0 + t x1,   x0 ~ N(0, I)
  rotations    : SO(3) geodesic            R_t = R0 exp(t log(R0^T R1)),
                 R0 ~ uniform SO(3)
The network does x1-PREDICTION: given (x_t, R_t, t) it predicts the clean frame
(x1, R1). This is the formulation whose Euler sampler we validated in numpy to
1e-15 with a perfect oracle, so if training drives the prediction toward the
data the sampler is correct by construction.

Equivariance
------------
The network is an IPA-lite frame transformer. It reads ONLY SE(3)-invariant
pairwise features -- the construction from the quaternion note:
  d_ij = R_i^T (x_j - x_i),   |x_j - x_i|,   R_ij = R_i^T R_j
and emits per-residue updates in each residue's LOCAL frame, which map back to
global coordinates equivariantly. A global rotation/translation of the input
therefore transforms the output the same way, and nothing depends on absolute
position. (Verify this empirically in step 3's eval.)

Run the CHEAP smoke test IN THE CONTAINER before any GPU job (seconds, CPU):
  python model/96_flow_backbone.py --smoke
Real training (GPU, enroot container) is driven by 97_train_backbone.sbatch.
"""

import argparse
import math
import os
import sys

import numpy as np
import torch
import torch.nn as nn

# =========================== SO(3) utilities ===============================
# All batched over arbitrary leading dims; vectors are (...,3), rots (...,3,3).

def skew(w):
    o = torch.zeros(w.shape[:-1] + (3, 3), dtype=w.dtype, device=w.device)
    o[..., 0, 1] = -w[..., 2]; o[..., 0, 2] = w[..., 1]
    o[..., 1, 0] = w[..., 2];  o[..., 1, 2] = -w[..., 0]
    o[..., 2, 0] = -w[..., 1]; o[..., 2, 1] = w[..., 0]
    return o


def so3_exp(w):
    th = torch.linalg.norm(w, dim=-1, keepdim=True).clamp_min(1e-8)      # (...,1)
    K = skew(w / th)
    th = th[..., None]                                                   # (...,1,1)
    I = torch.eye(3, dtype=w.dtype, device=w.device).expand(K.shape)
    return I + torch.sin(th) * K + (1 - torch.cos(th)) * (K @ K)


def so3_log(R):
    tr = R[..., 0, 0] + R[..., 1, 1] + R[..., 2, 2]
    c = ((tr - 1) / 2).clamp(-1.0, 1.0)
    th = torch.arccos(c)[..., None]                                      # (...,1)
    v = torch.stack([R[..., 2, 1] - R[..., 1, 2],
                     R[..., 0, 2] - R[..., 2, 0],
                     R[..., 1, 0] - R[..., 0, 1]], dim=-1)               # (...,3)
    small = th < 1e-4
    denom = (2 * torch.sin(th)).clamp_min(1e-8)
    return torch.where(small, v / 2, th * v / denom)


def rand_so3(shape, device, dtype=torch.float32):
    A = torch.randn(shape + (3, 3), device=device, dtype=dtype)
    Q, Rr = torch.linalg.qr(A)
    Q = Q * torch.sign(torch.diagonal(Rr, dim1=-2, dim2=-1))[..., None, :]
    det = torch.linalg.det(Q)                                           # +-1
    Q[..., :, 0] = Q[..., :, 0] * det[..., None]                        # force det +1
    return Q


def frames_from_coords(N, CA, C):
    """R has columns [e1,e2,e3]; e1 along C-CA, e2 in the N side of the plane."""
    e1 = torch.nn.functional.normalize(C - CA, dim=-1)
    v2 = N - CA
    u2 = v2 - (v2 * e1).sum(-1, keepdim=True) * e1
    e2 = torch.nn.functional.normalize(u2, dim=-1)
    e3 = torch.cross(e1, e2, dim=-1)
    return torch.stack([e1, e2, e3], dim=-1)                            # (...,3,3)


def gram_schmidt_6d(v6):
    """6D rotation representation -> proper rotation (Zhou et al. 2019)."""
    a1, a2 = v6[..., :3], v6[..., 3:]
    e1 = torch.nn.functional.normalize(a1, dim=-1)
    u2 = a2 - (a2 * e1).sum(-1, keepdim=True) * e1
    e2 = torch.nn.functional.normalize(u2, dim=-1)
    e3 = torch.cross(e1, e2, dim=-1)
    return torch.stack([e1, e2, e3], dim=-1)


# =========================== the network ===================================

class FrameLayer(nn.Module):
    def __init__(self, d):
        super().__init__()
        self.edge = nn.Sequential(nn.Linear(2 * d + 13, d), nn.ReLU(),
                                  nn.Linear(d, d))
        self.upd = nn.Sequential(nn.Linear(2 * d, d), nn.ReLU(), nn.Linear(d, d))
        self.norm = nn.LayerNorm(d)
        self.out = nn.Linear(d, 6)      # local dt (3) + so(3) dw (3)

    def forward(self, h, x, R, mask):
        B, L, _ = h.shape
        diff = x.unsqueeze(1) - x.unsqueeze(2)                 # [b,i,j] = xj - xi
        Rt = R.transpose(-1, -2)
        dloc = torch.einsum('bikm,bijm->bijk', Rt, diff)       # R_i^T (xj-xi)
        dist = torch.linalg.norm(diff, dim=-1, keepdim=True)
        Rij = torch.einsum('bimk,bjmn->bijkn', R, R).reshape(B, L, L, 9)
        edge_geo = torch.cat([dloc, dist, Rij], dim=-1)        # (B,L,L,13)
        hi = h.unsqueeze(2).expand(B, L, L, -1)
        hj = h.unsqueeze(1).expand(B, L, L, -1)
        msg = self.edge(torch.cat([hi, hj, edge_geo], dim=-1))  # (B,L,L,d)
        km = mask[:, None, :, None]                             # key mask over j
        msg = msg * km
        agg = msg.sum(2) / km.sum(2).clamp_min(1.0)             # masked mean
        h = self.norm(h + self.upd(torch.cat([h, agg], dim=-1)))
        d6 = self.out(h)                                        # (B,L,6)
        dt_loc, dw = d6[..., :3], d6[..., 3:]
        x = x + torch.einsum('bikm,bim->bik', R, dt_loc) * mask[..., None]
        R = R @ so3_exp(dw * mask[..., None])
        return h, x, R


class FrameNet(nn.Module):
    def __init__(self, d=128, layers=4):
        super().__init__()
        self.d = d
        self.tproj = nn.Linear(1, d)
        self.h0 = nn.Parameter(torch.zeros(d))
        self.layers = nn.ModuleList(FrameLayer(d) for _ in range(layers))

    def forward(self, x, R, t, mask):
        B, L, _ = x.shape
        h = self.h0[None, None, :].expand(B, L, -1).clone()
        h = h + self.tproj(t[:, None, None].expand(B, L, 1))
        h = h * mask[..., None]
        for lyr in self.layers:
            h, x, R = lyr(h, x, R, mask)
        return x, R                                            # pred x1, pred R1


# =========================== data ==========================================

SCALE = 10.0   # translations are stored in angstrom/SCALE so noise ~ N(0,1)


def load_backbones(path, min_len, max_len, device, limit=None):
    rows = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            rows.setdefault(p[0], []).append([float(v) for v in p[2:11]])
    chains = []
    for name, rs in rows.items():
        if len(rs) < min_len:
            continue
        a = np.asarray(rs[:max_len], dtype=np.float32)          # (L,9)
        N, CA, C = a[:, 0:3], a[:, 3:6], a[:, 6:9]
        t = torch.tensor(N); ca = torch.tensor(CA); c = torch.tensor(C)
        R = frames_from_coords(t, ca, c)
        x = (ca - ca.mean(0, keepdim=True)) / SCALE             # centered, scaled
        chains.append((x.to(device), R.to(device)))
        if limit and len(chains) >= limit:
            break
    return chains


def make_batch(chains, idx, device):
    L = max(chains[i][0].shape[0] for i in idx)
    B = len(idx)
    x1 = torch.zeros(B, L, 3, device=device)
    R1 = torch.eye(3, device=device).expand(B, L, 3, 3).clone()
    mask = torch.zeros(B, L, device=device)
    for b, i in enumerate(idx):
        xi, Ri = chains[i]
        n = xi.shape[0]
        x1[b, :n] = xi; R1[b, :n] = Ri; mask[b, :n] = 1.0
    return x1, R1, mask


# =========================== train / loss ==================================

# =================== chain-correlated (colored) noise prior ==================
# Instead of drawing the flow's base noise as independent N(0,I) per residue, draw
# it from N(0, Sigma) with Sigma = P (Lap+kI)^(-alpha) P over the COVALENT chain
# graph (path graph), normalised to unit mean marginal variance. Low (smooth,
# collective) chain modes carry most of the noise, so a sampled starting cloud is
# already a connected string and the model only has to refine a chain -- the
# principled fix for the under-sharpened translation field (C-N ~20 A with iid).
# The covalent graph is always known, so nothing about the fold leaks (a contact
# graph would). prior="iid" reproduces the original behaviour exactly.
_PRIOR_CACHE = {}

def colored_sqrt(L, alpha, kappa, device, dtype):
    key = (L, round(float(alpha), 4), round(float(kappa), 6))
    C = _PRIOR_CACHE.get(key)
    if C is None:
        lap = torch.zeros(L, L, dtype=torch.float64)
        for i in range(L - 1):
            lap[i, i] += 1; lap[i + 1, i + 1] += 1
            lap[i, i + 1] -= 1; lap[i + 1, i] -= 1
        evals, evecs = torch.linalg.eigh(lap)              # ascending; evals[0]~0
        w = (evals + kappa).clamp_min(1e-9) ** (-alpha / 2.0)
        w[0] = 0.0                                          # P: drop const (COM) mode
        Chalf = (evecs * w) @ evecs.T                      # Sigma^{1/2}
        # Normalise to the DATA's consecutive-CA spacing (3.8 A / SCALE = 0.38 units)
        # so the starting cloud is already a chain at the right local scale, not a
        # gas (iid) nor an arbitrary smooth curve. alpha/kappa then set only the
        # SHAPE of correlation beyond neighbours (higher alpha = more collective).
        Sig = Chalf @ Chalf.T
        di = Sig.diagonal()
        e2 = (di[:-1] + di[1:] - 2.0 * Sig.diagonal(1)).mean() * 3.0   # E|delta|^2 (3D)
        target = (3.8 / SCALE) ** 2
        Chalf = Chalf * (target / e2.clamp_min(1e-12)).sqrt()
        C = Chalf.to(dtype=dtype)
        _PRIOR_CACHE[key] = C
    return C.to(device=device)

def draw_prior(shape, prior, L, alpha, kappa, device, dtype=torch.float32):
    z = torch.randn(shape, device=device, dtype=dtype)      # (B,L,3)
    if prior == "colored" and L > 1:
        C = colored_sqrt(L, alpha, kappa, device, z.dtype)
        z = torch.einsum("lm,bmd->bld", C, z)
    return z


def flow_loss(net, x1, R1, mask, w_rot, w_dist=1.0, prior="iid", alpha=1.0, kappa=0.2,
              w_dgram=0.0, dgram_w=12):
    B, L, _ = x1.shape
    t = torch.rand(B, device=x1.device)
    tx = t[:, None, None]
    x0 = draw_prior(x1.shape, prior, L, alpha, kappa, x1.device, x1.dtype)
    R0 = rand_so3((B, L), x1.device, x1.dtype)
    xt = (1 - tx) * x0 + tx * x1
    tangent = so3_log(R0.transpose(-1, -2) @ R1)                # body-frame log
    Rt = R0 @ so3_exp(t[:, None, None] * tangent)
    px1, pR1 = net(xt, Rt, t, mask)
    m = mask[..., None]
    denom = mask.sum().clamp_min(1.0)
    lt = (((px1 - x1) * m) ** 2).sum() / (denom * 3)
    lr = (((pR1 - R1) * m[..., None]) ** 2).sum() / (denom * 9)
    # auxiliary chain-connectivity loss: match the CONSECUTIVE-CA distance of the
    # prediction to the data. Without this, the plain x1-MSE underfits the one
    # thing a valid backbone needs most -- neighbouring residues ~3.8 A apart --
    # and the sampler produces loose residue clouds (measured C-N ~20 A).
    ld = px1.new_zeros(())
    if w_dist > 0 and L > 1:
        dp = torch.linalg.norm(px1[:, 1:] - px1[:, :-1], dim=-1)   # (B,L-1)
        dt = torch.linalg.norm(x1[:, 1:] - x1[:, :-1], dim=-1)
        nm = mask[:, 1:] * mask[:, :-1]
        ld = (((dp - dt) * nm) ** 2).sum() / nm.sum().clamp_min(1.0)
    # distogram loss: match the LOCAL-band pairwise CA-CA distances (|i-j| <= dgram_w).
    # Consecutive spacing (ld) only fixes connectivity; the local distance pattern is
    # what encodes secondary structure (i,i+2/3/4) and local packing. Restricted to a
    # band because the full O(L^2) matrix is dominated by long-range pairs, which chase
    # global extent and blow the structure up (CPU-tested: full-matrix diverged).
    ldg = px1.new_zeros(())
    if w_dgram > 0 and L > 1:
        dpp = torch.cdist(px1, px1)                              # (B,L,L)
        dtt = torch.cdist(x1, x1)
        ii = torch.arange(L, device=x1.device)
        band = (ii[None, :] - ii[:, None]).abs() <= dgram_w     # (L,L) local window
        pm = mask[:, :, None] * mask[:, None, :] * band
        ldg = (((dpp - dtt) * pm) ** 2).sum() / pm.sum().clamp_min(1.0)
    return (lt + w_rot * lr + w_dist * ld + w_dgram * ldg,
            lt.item(), lr.item(), float(ld), float(ldg))


@torch.no_grad()
def sample_latent(net, n, L, steps, device, prior="iid", alpha=1.0, kappa=0.2):
    """Integrate noise->data in frame space (no PDB; step 3 writes structures)."""
    x = draw_prior((n, L, 3), prior, L, alpha, kappa, device)
    R = rand_so3((n, L), device)
    mask = torch.ones(n, L, device=device)
    dt = 1.0 / steps
    for k in range(steps):
        t = torch.full((n,), k * dt, device=device)
        rem = max(1 - k * dt, 1e-3)
        px1, pR1 = net(x, R, t, mask)
        x = x + dt * (px1 - x) / rem
        tang = so3_log(R.transpose(-1, -2) @ pR1) / rem
        R = R @ so3_exp(dt * tang)
    return x, R


def train(args):
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device:", dev)
    chains = load_backbones(args.data, args.min_len, args.max_len, dev)
    print("chains loaded:", len(chains))
    net = FrameNet(args.d_model, args.layers).to(dev)
    print("params:", sum(p.numel() for p in net.parameters()))
    opt = torch.optim.Adam(net.parameters(), lr=args.lr)
    n = len(chains)
    for ep in range(args.epochs):
        perm = torch.randperm(n).tolist()
        tot = tt = tr = td = tg = 0.0; nb = 0
        for s in range(0, n, args.batch):
            idx = perm[s:s + args.batch]
            x1, R1, mask = make_batch(chains, idx, dev)
            loss, lt, lr, ld, ldg = flow_loss(net, x1, R1, mask, args.w_rot, args.w_dist,
                                              args.prior, args.alpha, args.kappa, args.w_dgram)
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(net.parameters(), 1.0)
            opt.step()
            tot += loss.item(); tt += lt; tr += lr; td += ld; tg += ldg; nb += 1
        if ep % args.log_every == 0 or ep == args.epochs - 1:
            print("epoch %4d  loss %.4f  trans %.4f  rot %.4f  dist %.4f  dgram %.4f"
                  % (ep, tot / nb, tt / nb, tr / nb, td / nb, tg / nb))
    if args.save:
        torch.save({"state": net.state_dict(), "d_model": args.d_model,
                    "layers": args.layers, "scale": SCALE,
                    "prior": args.prior, "alpha": args.alpha, "kappa": args.kappa},
                   args.save)
        print("saved", args.save)


# =========================== smoke test ====================================

def smoke():
    """Tiny CPU end-to-end check: shapes, a few training steps with loss going
    down, and a latent sampling pass with no NaNs. Run IN THE CONTAINER before
    any GPU job -- it catches wiring bugs in seconds without a GPU."""
    torch.manual_seed(0)
    dev = torch.device("cpu")
    # synthetic chains: random frames + random-walk CA
    chains = []
    for _ in range(12):
        L = int(torch.randint(20, 40, (1,)))
        x = torch.cumsum(torch.randn(L, 3) * 0.1, 0)
        R = rand_so3((L,), dev)
        chains.append((x, R))
    net = FrameNet(d=32, layers=3)
    opt = torch.optim.Adam(net.parameters(), lr=1e-3)
    first = last = None
    for step in range(30):
        idx = list(range(len(chains)))
        x1, R1, mask = make_batch(chains, idx, dev)
        loss, lt, lr, ld, ldg = flow_loss(net, x1, R1, mask, 1.0, 1.0, w_dgram=1.0)
        opt.zero_grad(); loss.backward(); opt.step()
        if step == 0: first = loss.item()
        last = loss.item()
    xs, Rs = sample_latent(net, 2, 25, 20, dev)
    ortho = (Rs.transpose(-1, -2) @ Rs
             - torch.eye(3)).abs().max().item()
    print("smoke: loss %.4f -> %.4f   (down=%s)"
          % (first, last, last < first))
    print("smoke: sample x finite=%s R finite=%s  max|RtR-I|=%.2e"
          % (torch.isfinite(xs).all().item(),
             torch.isfinite(Rs).all().item(), ortho))
    ok = (last < first and torch.isfinite(xs).all()
          and torch.isfinite(Rs).all() and ortho < 1e-3)
    print("smoke:", "PASS" if ok else "FAIL")
    return 0 if ok else 2


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=128)
    ap.add_argument("--d-model", type=int, default=128)
    ap.add_argument("--layers", type=int, default=4)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--w-rot", type=float, default=1.0)
    ap.add_argument("--w-dist", type=float, default=1.0,
                    help="weight on the consecutive-CA distance (chain) loss")
    ap.add_argument("--w-dgram", type=float, default=0.0,
                    help="weight on the full pairwise-CA distogram loss (SS/packing/compactness)")
    ap.add_argument("--prior", choices=["iid", "colored"], default="iid",
                    help="base noise: iid N(0,I) or chain-correlated (covalent graph)")
    ap.add_argument("--alpha", type=float, default=2.0, help="colored prior smoothness exponent")
    ap.add_argument("--kappa", type=float, default=0.1, help="colored prior graph-Laplacian shift")
    ap.add_argument("--log-every", type=int, default=10)
    ap.add_argument("--save", default=None)
    args = ap.parse_args(argv)
    if args.smoke:
        return smoke()
    if not args.data:
        ap.error("--data is required unless --smoke")
    train(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
