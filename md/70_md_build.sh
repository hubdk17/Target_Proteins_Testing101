#!/bin/bash
# Build an AMBER system from a two-chain complex PDB. Runs INSIDE the nvhpc
# container. Usage (inside container): bash 70_md_build.sh <complex.pdb> <outdir>
# Produces <outdir>/sys.prmtop, sys.inpcrd, sys_built.pdb, and runs both build
# gates on the built system BEFORE any GPU time is spent.
#
# NO PYTHON ANYWHERE. The nvhpc container ships no python at all (no
# /usr/bin/python3, nothing on PATH, and the AMBER 24 archive bundles no
# interpreter either), so `pdb4amber` -- a #!/usr/bin/python3 script -- exits
# 126 on its first line and the build dies before touching the structure.
# Measured 2026-11-19/20, job 188654. The cleaning and the gates are therefore
# awk: 71_clean_pdb.awk and 71_check_built.awk.
set -euo pipefail
PDB="$1"; OUT="$2"; mkdir -p "$OUT"
# Derive the scripts root from this script's own location, so the clone can
# live anywhere (~/rfd/scripts, ~/rfd/scripts_md, ...) without edits.
S="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$S/md/_amber_env.sh"

command -v tleap >/dev/null || { echo "ERROR: tleap not found on PATH under $AMBERHOME"; exit 1; }

echo "=== clean (awk: keep ATOM, drop H/waters/altlocs, TER at every chain change) ==="
awk -f "$S/md/71_clean_pdb.awk" "$PDB" > "$OUT/clean.pdb"
echo "  input ATOM records : $(grep -c '^ATOM' "$PDB")"
echo "  cleaned ATOM records: $(grep -c '^ATOM' "$OUT/clean.pdb")  TER records: $(grep -c '^TER' "$OUT/clean.pdb")"
# A TER between the chains is what stops tleap fusing target and binder. If the
# cleaner produced none, stop here rather than build a complex that can never
# come apart.
[ "$(grep -c '^TER' "$OUT/clean.pdb")" -ge 2 ] || { echo "ERROR: fewer than 2 TER records -- chains not separated"; exit 1; }

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
[ -f "$OUT/sys.prmtop" ] || { echo "ERROR: tleap produced no prmtop -- see $OUT/tleap.log"; exit 1; }

echo "=== build gates on the built system (chirality + chain separation) ==="
awk -f "$S/md/71_check_built.awk" "$OUT/sys_built.pdb"
