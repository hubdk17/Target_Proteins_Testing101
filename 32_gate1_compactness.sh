#!/bin/bash
# GATE 1, compactness, applied per target across the whole batch.
#
# rho = Rg(binder CAs) / (2.2 * N^0.38); keep rho < 1.30.
# Formula and reasoning: /mnt/project-files/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md
#
# CPU only, login node, seconds.
#   cd ~/rfd && bash 32_gate1_compactness.sh            # all targets in targets.tsv
#   cd ~/rfd && bash 32_gate1_compactness.sh pd1        # one target
#
# Writes, per target:  gate1/<name>.tsv   full table, every design
#                      gate1/<name>.pass  the survivor paths, one per line
# Prints a one-line summary per target so a bad target is visible immediately.

set -u
ROOT="$HOME/rfd"
TABLE="$ROOT/inputs/targets.tsv"
TAG="${TAG:-b1}"
CUT="${CUT:-1.30}"

test -s "$TABLE" || { echo "$TABLE missing - run 02_prep_targets.py"; exit 1; }
mkdir -p "$ROOT/gate1"

if [ $# -gt 0 ]; then NAMES="$*"; else NAMES=$(cut -f1 "$TABLE"); fi

printf "%-8s %6s %6s %7s %7s %s\n" target designs pass rate medrho note

for NAME in $NAMES; do
  OUT="$ROOT/gate1/${NAME}.tsv"
  PASS="$ROOT/gate1/${NAME}.pass"
  : > "$OUT"; : > "$PASS"

  FILES=$(ls "$ROOT/outputs/${NAME}_${TAG}"_*.pdb 2>/dev/null)
  if [ -z "$FILES" ]; then
    printf "%-8s %6s %6s %7s %7s %s\n" "$NAME" 0 0 - - "NO OUTPUT - did stage 1 run?"
    continue
  fi

  for f in $FILES; do
    # Binder is chain B. Rg from CA atoms only; N from the same atoms, so rho is
    # length-free and a rod lands near 2.2 whatever its length.
    awk -v F="$f" '
      substr($0,1,4)=="ATOM" && substr($0,13,4)==" CA " && substr($0,22,1)=="B" {
        x=substr($0,31,8)+0; y=substr($0,39,8)+0; z=substr($0,47,8)+0
        sx+=x; sy+=y; sz+=z; q+=x*x+y*y+z*z; n++
      }
      END {
        if (n < 10) { printf "%s\t%d\tNA\tNA\tSHORT\n", F, n+0; exit }
        rg = sqrt(q/n - (sx/n)^2 - (sy/n)^2 - (sz/n)^2)
        printf "%s\t%d\t%.2f\t%.3f\t\n", F, n, rg, rg/(2.2*n^0.38)
      }' "$f" >> "$OUT"
  done

  # Survivors. Kept separate from the summary so the .pass file is written by
  # exactly one thing.
  awk -F'\t' -v C="$CUT" '$4 != "NA" && $4+0 < C { print $1 }' "$OUT" > "$PASS"

  NTOT=$(awk -F'\t' '$4 != "NA"' "$OUT" | wc -l | tr -d ' ')
  NPASS=$(wc -l < "$PASS" | tr -d ' ')
  if [ "$NTOT" -eq 0 ]; then
    printf "%-8s %6s %6s %7s %7s %s\n" "$NAME" 0 0 - - "no measurable designs"
    continue
  fi
  RATE=$(awk -v k="$NPASS" -v t="$NTOT" 'BEGIN{printf "%.0f%%", 100*k/t}')
  # Median without gawk's asort, which is not portable.
  MED=$(awk -F'\t' '$4 != "NA" {print $4}' "$OUT" | sort -n | awk '
    {v[NR]=$1} END{ if(NR==0){print "-"} else if(NR%2){printf "%.2f", v[(NR+1)/2]}
                    else {printf "%.2f", (v[NR/2]+v[NR/2+1])/2} }')
  NOTE=""
  # A median near 2.2 means the whole target came out as rods, which points at the
  # run (potential missing, hotspots scattered) rather than at luck.
  case "$MED" in
    2.0*|2.1*|2.2*|2.3*) NOTE="median is rod-like - check the ROG potential and hotspots" ;;
  esac
  [ "${NPASS:-0}" = "0" ] && NOTE="NOTHING PASSED - do not run stage 2 on this target"
  printf "%-8s %6s %6s %7s %7s %s\n" "$NAME" "$NTOT" "$NPASS" "$RATE" "$MED" "$NOTE"
  rm -f "$ROOT/gate1/.sum_$NAME"
done

echo
echo "survivor lists: $ROOT/gate1/<target>.pass    full tables: $ROOT/gate1/<target>.tsv"
echo "cut-off rho < $CUT   (expect roughly 30-60% to pass; 9% means the potential is off)"
