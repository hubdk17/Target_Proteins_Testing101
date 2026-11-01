#!/usr/bin/env bash
# IL-6 STEP 1 -- run this ON THE LOGIN NODE (bcmnode1). No GPU, no container.
#   bash 01_fetch_il6.sh
#
# Does three things:
#   1. checks the RFdiffusion environment you already built is still there
#   2. downloads human IL-6 (PDB 1ALU) and strips it to protein atoms only
#   3. prints the chains and residue numbering, which is what the binder
#      contig string has to be built from
#
# Safe to re-run. Nothing here touches a GPU and nothing here is submitted
# to Slurm, so it costs you no allocation.

set -euo pipefail

RFD_ROOT="${RFD_ROOT:-$HOME/rfd}"
INPUTS="$RFD_ROOT/inputs"
mkdir -p "$INPUTS" "$RFD_ROOT/logs" "$RFD_ROOT/outputs"

echo "=== 1/3 checking the environment from the trypsin run ==="
for f in "$RFD_ROOT/venv/.complete" \
         "$RFD_ROOT/models/Complex_base_ckpt.pt" \
         "$RFD_ROOT/RFdiffusion/scripts/run_inference.py"; do
  if [ -e "$f" ]; then echo "  ok      $f"
  else                 echo "  MISSING $f"; fi
done
echo "  (if anything says MISSING, stop and tell me - do not run the rest)"

echo
echo "=== 2/3 downloading IL-6 (PDB 1ALU, human interleukin-6) ==="
if [ -s "$INPUTS/1ALU.pdb" ]; then
  echo "  already have 1ALU.pdb"
else
  wget -q --show-progress -O "$INPUTS/1ALU.pdb" https://files.rcsb.org/download/1ALU.pdb
fi

# Protein atoms only: drop waters, ions and any ligand. Also drop the second
# copy of any atom that was modelled in two alternate positions (altloc B+),
# keeping the first, because RFdiffusion wants one atom per position.
awk '
  /^ATOM/  { alt = substr($0,17,1); if (alt == " " || alt == "A") print; next }
  /^TER/   { print }
' "$INPUTS/1ALU.pdb" > "$INPUTS/il6_clean.pdb"
echo "END" >> "$INPUTS/il6_clean.pdb"
echo "  wrote $INPUTS/il6_clean.pdb"

echo
echo "=== 3/3 what is actually in the file ==="
echo "--- chains, and how many residues each has ---"
awk '/^ATOM/ && substr($0,13,4) == " CA " { print substr($0,22,1) }' \
  "$INPUTS/il6_clean.pdb" | sort | uniq -c | sed 's/^/  chain /'

echo "--- first and last residue number per chain ---"
awk '/^ATOM/ && substr($0,13,4) == " CA " {
       c = substr($0,22,1); n = substr($0,23,4) + 0;
       if (!(c in lo)) lo[c] = n; hi[c] = n;
     }
     END { for (c in lo) printf "  chain %s: %d to %d\n", c, lo[c], hi[c] }' \
  "$INPUTS/il6_clean.pdb"

echo "--- breaks in the numbering (missing residues) ---"
awk '/^ATOM/ && substr($0,13,4) == " CA " {
       c = substr($0,22,1); n = substr($0,23,4) + 0;
       if (c == pc && n - pn > 1) { printf "  chain %s: %d -> %d  (%d missing)\n", c, pn, n, n-pn-1; found++ }
       pc = c; pn = n;
     }
     END { if (!found) print "  none - the chain is continuous" }' \
  "$INPUTS/il6_clean.pdb"

echo
echo "=== step 1 done. Paste everything from '=== 1/3' down to here back to me. ==="
