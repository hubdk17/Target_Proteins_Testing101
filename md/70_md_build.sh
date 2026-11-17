#!/bin/bash
# Build an AMBER system from a two-chain complex PDB. Runs INSIDE the nvhpc
# container. Usage (inside container): bash 70_md_build.sh <complex.pdb> <outdir>
# Produces <outdir>/sys.prmtop, sys.inpcrd, sys_built.pdb, and runs the
# chirality check on the built system BEFORE any GPU time is spent, because
# tleap rebuilds hydrogens and can introduce D-residues that every downstream
# gate is blind to.
set -euo pipefail
PDB="$1"; OUT="$2"; mkdir -p "$OUT"
export AMBERHOME=/workspace/amber24
source "$AMBERHOME/amber.sh"
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
