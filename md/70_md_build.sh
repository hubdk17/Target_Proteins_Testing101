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
python3 /workspace/rfd/scripts/44_chirality_check.py "$OUT/sys_built.pdb" \
  || echo "CHIRALITY CHECK SCRIPT MISSING OR FAILED -- do not run MD until resolved"
