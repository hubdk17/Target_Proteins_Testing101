#!/bin/bash
# Build an AMBER system from a two-chain complex PDB. Runs INSIDE the nvhpc
# container. Usage (inside container): bash 70_md_build.sh <complex.pdb> <outdir>
# Produces <outdir>/sys.prmtop, sys.inpcrd, sys_built.pdb, and runs the
# chirality check on the built system BEFORE any GPU time is spent, because
# tleap rebuilds hydrogens and can introduce D-residues that every downstream
# gate is blind to.
set -euo pipefail
PDB="$1"; OUT="$2"; mkdir -p "$OUT"
source /workspace/rfd/scripts/md/_amber_env.sh
# clean: add TER between chains, strip waters/alt locs, standardise
pdb4amber -i "$PDB" -o "$OUT/clean.pdb" --nohyd --dry
cat > "$OUT/build.in" <<TLEAP
source leaprc.protein.ff14SB
source leaprc.water.tip3p
mol = loadpdb $OUT/clean.pdb
check mol
addions mol Na+ 0
addions mol Cl- 0
solvateoct mol TIP3PBOX 10.0
saveamberparm mol $OUT/sys.prmtop $OUT/sys.inpcrd
savepdb mol $OUT/sys_built.pdb
quit
TLEAP
tleap -f "$OUT/build.in" | tee "$OUT/tleap.log"
echo "=== chirality check on built system (must report 0 D-residues) ==="
# self-contained: N-C-CA-CB improper ~ -122 deg for L, positive for D
# (calibrated on 1YCR: 94/94 L-residues at -122.3 +/- 3.0 deg)
# `if ! ...; then exit 1` so the gate genuinely aborts the job (a swallowed
# check is exactly the "check that cannot fail" trap we keep hitting).
if ! python3 - "$OUT/sys_built.pdb" <<'PYEOF'
import sys, math
res={}
for ln in open(sys.argv[1]):
    if not ln.startswith("ATOM"): continue
    k=(ln[21],int(ln[22:26])); res.setdefault(k,[ln[17:20].strip(),{}])[1][ln[12:16].strip()]=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
def sub(a,b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2])
def crs(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dt(a,b): return sum(x*y for x,y in zip(a,b))
def imp(p1,p2,p3,p4):
    b1,b2,b3=sub(p2,p1),sub(p3,p2),sub(p4,p3); n1,n2=crs(b1,b2),crs(b2,b3)
    nb=math.sqrt(dt(b2,b2));
    if nb==0: return None
    m=crs(n1,(b2[0]/nb,b2[1]/nb,b2[2]/nb)); return math.degrees(math.atan2(dt(m,n2),dt(n1,n2)))
nd=0; n=0; bad=[]
for k,(rn,a) in res.items():
    if rn=="GLY" or not {"N","CA","C","CB"}<=a.keys(): continue
    v=imp(a["N"],a["C"],a["CA"],a["CB"]);
    if v is None: continue
    n+=1
    if v>0: nd+=1; bad.append("%s%s%d"%(k[0],rn,k[1]))
print("chiral residues checked:",n," D-residues:",nd, (" "+",".join(bad[:10]) if bad else ""))
sys.exit(1 if nd else 0)
PYEOF
then
  :
else
  echo "D-RESIDUES FOUND in built system -- aborting before MD"
  exit 1
fi

echo "=== chain-separation check (the two chains must NOT be fused) ==="
# ref_p53.pdb and the design complexes carry no TER record between chain A and
# chain B. pdb4amber normally inserts one, but if it does not, tleap happily
# makes a peptide bond from the last target residue to the first binder
# residue -- a single covalent chain that can never come apart. Every
# downstream read (binder RMSD, interface contacts) would still look sane,
# so this is gated explicitly rather than assumed.
if ! python3 - "$OUT/sys_built.pdb" <<'PYEOF'
import sys, math
# ordered protein residues: (chain, resseq) -> {atom: xyz}
order=[]; res={}
skip={"WAT","HOH","Na+","Cl-","NA","CL"}
for ln in open(sys.argv[1]):
    if not ln.startswith("ATOM"): continue
    rn=ln[17:20].strip()
    if rn in skip: continue
    k=(ln[21], int(ln[22:26]))
    if k not in res: res[k]={}; order.append(k)
    res[k][ln[12:16].strip()]=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
breaks=[]
for a,b in zip(order, order[1:]):
    pa, pb = res[a].get("C"), res[b].get("N")
    if pa is None or pb is None: continue
    d=math.dist(pa,pb)
    if d>2.0: breaks.append((order.index(b)+1, round(d,2)))
print("protein residues:", len(order), " chain breaks found:", len(breaks),
      "at built-residue index(es):", [b[0] for b in breaks[:6]])
sys.exit(0 if breaks else 1)
PYEOF
then
  :
else
  echo "NO CHAIN BREAK in built system -- target and binder were fused into one"
  echo "chain by tleap. Aborting before MD. Fix: insert a TER record between"
  echo "the chains in the input PDB and rebuild."
  exit 1
fi
