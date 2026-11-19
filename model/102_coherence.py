#!/usr/bin/env python3
"""Structural-coherence scorecard for generated backbones (stdlib only).
Goes beyond C-N connectivity: bond geometry, steric clashes, secondary-structure
content, and compactness -- the things that separate a real fold from a connected
but incoherent chain. Run on a directory of sampled PDBs.
  python3 102_coherence.py <dir-of-pdbs>
"""
import sys, os, math, glob

IDEAL = {"N_CA": 1.458, "CA_C": 1.525, "C_N": 1.329, "CA_CA": 3.803}

def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def dot(a, b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def nrm(a): return math.sqrt(dot(a, a))
def dist(a, b): return nrm(sub(a, b))

def dih(p0, p1, p2, p3):
    b0 = sub(p0, p1); b1 = sub(p2, p1); b2 = sub(p3, p2)
    n = nrm(b1)
    if n == 0: return None
    u = (b1[0]/n, b1[1]/n, b1[2]/n)
    v = sub(b0, tuple(dot(b0, u)*c for c in u))
    w = sub(b2, tuple(dot(b2, u)*c for c in u))
    return math.degrees(math.atan2(dot(cross(u, v), w), dot(v, w)))

def parse(path):
    order = []; res = {}
    for ln in open(path):
        if not ln.startswith("ATOM"): continue
        atom = ln[12:16].strip(); c = ln[21]
        try:
            rs = int(ln[22:26]); x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
        except ValueError: continue
        key = (c, rs)
        if key not in res: res[key] = {}; order.append(key)
        res[key][atom] = (x, y, z)
    return order, res

def score(path):
    order, res = parse(path)
    cas = [res[k]["CA"] for k in order if "CA" in res[k]]
    n = len(cas)
    cn = []; caca = []
    for i in range(len(order)-1):
        a = res[order[i]]; b = res[order[i+1]]
        if "C" in a and "N" in b: cn.append(dist(a["C"], b["N"]))
        if "CA" in a and "CA" in b: caca.append(dist(a["CA"], b["CA"]))
    # steric clashes: non-adjacent CA pairs closer than 3.5 A
    clash = 0
    for i in range(n):
        for j in range(i+2, n):
            if dist(cas[i], cas[j]) < 3.5: clash += 1
    # secondary structure by phi/psi
    H = E = L = 0; nss = 0
    for i in range(1, len(order)-1):
        rp = res[order[i-1]]; ri = res[order[i]]; rn = res[order[i+1]]
        if not all(k in ri for k in ("N", "CA", "C")) or "C" not in rp or "N" not in rn: continue
        phi = dih(rp["C"], ri["N"], ri["CA"], ri["C"]); psi = dih(ri["N"], ri["CA"], ri["C"], rn["N"])
        if phi is None or psi is None: continue
        nss += 1
        if -145 <= phi <= -30 and -70 <= psi <= 20: H += 1
        elif -180 <= phi <= -40 and (90 <= psi <= 180 or -180 <= psi <= -150): E += 1
        else: L += 1
    # compactness
    cx = tuple(sum(p[d] for p in cas)/n for d in range(3))
    rg = math.sqrt(sum(dist(p, cx)**2 for p in cas)/n)
    gate = 2.2*(n**0.38)
    cnm = sum(cn)/len(cn) if cn else float("nan")
    cam = sum(caca)/len(caca) if caca else float("nan")
    return dict(n=n, cn=cnm, caca=cam, clash=clash,
                H=100*H/nss if nss else 0, E=100*E/nss if nss else 0,
                L=100*L/nss if nss else 0, rgr=rg/gate if gate else 0)

files = sorted(glob.glob(os.path.join(sys.argv[1], "*.pdb")))
if not files:
    print("no PDBs in", sys.argv[1]); sys.exit(0)
agg = {}
print(f"--- coherence scorecard: {sys.argv[1]} ({len(files)} structures) ---")
for f in files:
    s = score(f)
    for k, v in s.items(): agg[k] = agg.get(k, 0) + v
m = {k: v/len(files) for k, v in agg.items()}
print(f"  C-N {m['cn']:.2f} (ideal 1.33)   CA-CA {m['caca']:.2f} (ideal 3.80)   "
      f"clashes/struct {m['clash']:.1f}")
print(f"  SS: helix {m['H']:.0f}%  strand {m['E']:.0f}%  loop {m['L']:.0f}%   "
      f"Rg/compact-gate {m['rgr']:.2f} (<1.3 good)")
ok = (abs(m['cn']-1.33) < 0.3 and abs(m['caca']-3.80) < 0.6 and m['clash'] < 2
      and (m['H']+m['E']) > 40 and m['rgr'] < 1.4)
print(f"  COHERENT-ish: {ok}  (connected + clash-free + has SS + compact)")
