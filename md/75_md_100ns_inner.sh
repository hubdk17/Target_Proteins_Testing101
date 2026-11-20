#!/bin/bash
# 100 ns production MD of a two-chain complex, in 10 ns chunks.
# Runs INSIDE the nvhpc container; called by 75_md_100ns.sbatch.
#   bash 75_md_100ns_inner.sh <complex.pdb> <outdir> [nchunks=10]
#
# Resumable and idempotent: build, min/heat/equil and every finished chunk are
# skipped on a re-run, so re-submitting the identical command after a crash,
# a timeout or a node failure continues where it stopped. A chunk is only
# counted as finished once BOTH its .nc and .rst have been renamed off .tmp,
# so a half-written chunk is redone rather than silently treated as 10 ns.
set -euo pipefail
PDB="$1"; OUT="$2"; NCHUNK="${3:-10}"
S="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$S/md/_amber_env.sh"
cd "$S"
mkdir -p "$OUT"

echo "=== GPU ==="; nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

if [ ! -f "$OUT/sys.prmtop" ]; then
  echo "=== build (tleap + chirality + chain-break gates) ==="
  bash md/70_md_build.sh "$S/$PDB" "$S/$OUT"
else
  echo "[skip] build: $OUT/sys.prmtop exists"
fi

cd "$S/$OUT"
P="$S/md"

if [ ! -f equil.rst ]; then
  echo "=== minimisation (also the pmemd.cuda / GPU-load test) ==="
  pmemd.cuda -O -i "$P/min.in"   -p sys.prmtop -c sys.inpcrd -ref sys.inpcrd -o min.out   -r min.rst   -x min.nc
  echo "=== heat 0->300 K ==="
  pmemd.cuda -O -i "$P/heat.in"  -p sys.prmtop -c min.rst    -ref min.rst    -o heat.out  -r heat.rst  -x heat.nc
  echo "=== equilibration 500 ps NPT ==="
  pmemd.cuda -O -i "$P/equil.in" -p sys.prmtop -c heat.rst   -ref heat.rst   -o equil.out -r equil.rst -x equil.nc
else
  echo "[skip] min/heat/equil: equil.rst exists"
fi

prev=equil.rst
for i in $(seq 1 "$NCHUNK"); do
  tag=$(printf "%02d" "$i")
  if [ -f "prod_$tag.rst" ] && [ -f "prod_$tag.nc" ]; then
    echo "[skip] chunk $tag already done"; prev="prod_$tag.rst"; continue
  fi
  rm -f "prod_$tag.rst.tmp" "prod_$tag.nc.tmp"
  t0=$SECONDS
  echo "=== production chunk $tag / $NCHUNK (10 ns) from $prev ==="
  pmemd.cuda -O -i "$P/prod.in" -p sys.prmtop -c "$prev" \
             -o "prod_$tag.out" -r "prod_$tag.rst.tmp" -x "prod_$tag.nc.tmp"
  mv "prod_$tag.nc.tmp"  "prod_$tag.nc"
  mv "prod_$tag.rst.tmp" "prod_$tag.rst"
  el=$((SECONDS-t0))
  echo "[chunk $tag] 10 ns in ${el}s  (~$(awk -v e="$el" 'BEGIN{printf "%.0f", 10*86400/e}') ns/day)"
  prev="prod_$tag.rst"
done

done_ns=$(( $(ls prod_??.nc 2>/dev/null | wc -l) * 10 ))
echo "=== DONE: ${done_ns} ns of production in $(pwd) ==="
