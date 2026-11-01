#!/usr/bin/env bash
# STEP 2 -- run this ON THE LOGIN NODE (bcmnode1). No GPU, no container.
#   bash 01_fetch.sh
#
# Downloads what the GPU job needs into ~/rfd so the compute node never has to:
#   - RFdiffusion source
#   - the two model checkpoints we use (~4 GB)
#   - trypsin (PDB 2PTN), cleaned, plus its residue numbering
# Safe to re-run: it skips anything already present.

set -euo pipefail

RFD_ROOT="${RFD_ROOT:-$HOME/rfd}"
REPO="$RFD_ROOT/RFdiffusion"
MODELS="$RFD_ROOT/models"
INPUTS="$RFD_ROOT/inputs"
mkdir -p "$MODELS" "$INPUTS" "$RFD_ROOT/logs" "$RFD_ROOT/outputs"

echo "=== 1/3 RFdiffusion source ==="
if [ -d "$REPO/.git" ]; then echo "already cloned"; else
  git clone https://github.com/RosettaCommons/RFdiffusion.git "$REPO"
fi

echo "=== 2/3 model weights ==="
# The repo's download_models.sh pulls all ten checkpoints (~25 GB). We need two.
BASE=http://files.ipd.uw.edu/pub/RFdiffusion
get() {
  if [ -s "$MODELS/$2" ]; then echo "  have $2"; else
    echo "  getting $2"
    wget -q --show-progress -O "$MODELS/$2.part" "$1/$2" && mv "$MODELS/$2.part" "$MODELS/$2"
  fi
}
get $BASE/6f5902ac237024bdd0c176cb93063dc4 Base_ckpt.pt
get $BASE/e29311f6f1bf1af907f9ef9f44b8328b Complex_base_ckpt.pt
ls -lh "$MODELS"

echo "=== 3/3 trypsin target (PDB 2PTN, bovine beta-trypsin) ==="
if [ ! -s "$INPUTS/2PTN.pdb" ]; then
  wget -q -O "$INPUTS/2PTN.pdb" https://files.rcsb.org/download/2PTN.pdb
fi
# Protein atoms only: drops waters, the Ca2+ ion and the benzamidine ligand.
grep -E '^(ATOM|TER)' "$INPUTS/2PTN.pdb" > "$INPUTS/trypsin_clean.pdb"
echo "END" >> "$INPUTS/trypsin_clean.pdb"

echo
echo "--- residue numbering (needed for the binder contig) ---"
# 2PTN uses chymotrypsinogen numbering: it does not start at 1 and has gaps,
# so the contig string must use these numbers, not 1..N.
echo -n "residue count: "; grep -c ' CA ' "$INPUTS/trypsin_clean.pdb"
echo -n "first residue: "; grep ' CA ' "$INPUTS/trypsin_clean.pdb" | head -1 | cut -c22-27
echo -n "last residue:  "; grep ' CA ' "$INPUTS/trypsin_clean.pdb" | tail -1 | cut -c22-27
echo "gaps in numbering:"
grep ' CA ' "$INPUTS/trypsin_clean.pdb" | cut -c23-26 \
  | awk '{n=$1+0; if (p && n-p>1) printf "  %d -> %d\n", p, n; p=n} END{if(!f) print "  (see above; none printed means none found)"}'

echo
echo "=== step 2 done.  RFD_ROOT = $RFD_ROOT ==="
