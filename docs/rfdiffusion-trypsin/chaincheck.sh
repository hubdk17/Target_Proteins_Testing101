#!/usr/bin/env bash
# Sanity-check any PDB: per-chain residue count, numbering range, and real
# chain breaks (consecutive CA atoms further apart than a peptide bond allows).
#
#   bash chaincheck.sh ~/rfd/outputs/uncond_50aa_0.pdb
#   bash chaincheck.sh ~/rfd/inputs/trypsin_clean.pdb
#
# Consecutive CA atoms sit ~3.8 A apart, so >4.5 A is a genuine break. This is
# how to tell a gap in NUMBERING (harmless, common in crystal structures that
# follow a reference numbering) from a gap in the CHAIN (which a contig string
# must be split across).

[ $# -eq 1 ] || { echo "usage: bash chaincheck.sh <file.pdb>"; exit 1; }
[ -f "$1" ] || { echo "no such file: $1"; exit 1; }

awk '
/^ATOM/ && substr($0,13,4)==" CA " {
  ch = substr($0,22,1); n = substr($0,23,4)+0
  x = substr($0,31,8)+0; y = substr($0,39,8)+0; z = substr($0,47,8)+0
  cnt[ch]++
  if (!(ch in firstres)) { firstres[ch] = n; order[++nch] = ch }
  lastres[ch] = n
  if (pch == ch) {
    d = sqrt((x-px)^2 + (y-py)^2 + (z-pz)^2)
    if (d > 4.5) { printf "  CHAIN BREAK %s: %d -> %d  (%.2f A)\n", ch, pn, n, d; breaks++ }
    else if (n - pn > 1) { printf "  numbering gap only %s: %d -> %d  (%.2f A, continuous)\n", ch, pn, n, d; gaps++ }
  }
  pch = ch; pn = n; px = x; py = y; pz = z
}
END {
  for (i = 1; i <= nch; i++) {
    c = order[i]
    printf "chain %s: %d residues, numbered %d..%d\n", c, cnt[c], firstres[c], lastres[c]
  }
  printf "%d real chain break(s), %d numbering-only gap(s)\n", breaks+0, gaps+0
  if (breaks+0 == 0) print "-> chain is continuous; a single contig range is correct"
  else print "-> split the contig at the real break(s)"
}' "$1"
