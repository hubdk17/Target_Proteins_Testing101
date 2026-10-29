#!/bin/bash
# GATE 2, part 3: read the self-consistency scores and apply the 2 A cut.
# CPU only, login node, instant. Run after 40_gate2_selfconsistency.sbatch.
#
#   cd ~/rfd && bash 41_gate2_report.sh            # all targets
#   cd ~/rfd && bash 41_gate2_report.sh pd1
#
# Writes gate2/<name>.pass listing the sequences that passed, and prints the
# per-backbone table. The deciding number per backbone is its BEST (lowest)
# binder_aligned_rmsd across its sequences.
#
# pLDDT is printed for information only. DO NOT rank by it: on IL-6 the two
# sequences that folded into rods scored 93-95 while the six correct folds scored
# 75-79, so confidence is anti-correlated with correctness here.

set -u
ROOT="$HOME/rfd"
CUT="${CUT:-2.0}"
TABLE="$ROOT/inputs/targets.tsv"
mkdir -p "$ROOT/gate2"

if [ $# -gt 0 ]; then NAMES="$*"; elif [ -s "$TABLE" ]; then NAMES=$(cut -f1 "$TABLE");
else echo "no targets.tsv and no target named"; exit 1; fi

for NAME in $NAMES; do
  # Newest score file for this target; a run writes a fresh timestamped name.
  SC=$(ls -t "$ROOT"/sc_${NAME}_*.sc 2>/dev/null | head -1)
  echo "============================================================"
  if [ -z "$SC" ]; then
    echo "$NAME: no score file (sc_${NAME}_*.sc) - has gate 2 run?"
    continue
  fi
  echo "$NAME   $(basename "$SC")"

  # Columns in a dl_binder_design monomer score line:
  #   $2 binder_aligned_rmsd ... $11 description
  # Taken from the IL-6 run's header; if the header differs, fix the indices here.
  HDR=$(grep -m1 '^SCORE:' "$SC")
  case "$HDR" in
    *binder_aligned_rmsd*) : ;;
    *) echo "  unexpected header, check column order:"; echo "  $HDR" ;;
  esac

  grep '^SCORE:' "$SC" | awk -v C="$CUT" -v P="$ROOT/gate2/$NAME.pass" '
    $2 == "binder_aligned_rmsd" { next }           # header line
    $2 + 0 != $2 { next }                          # not a number
    {
      desc = $NF
      bb = desc; sub(/_s[0-9]+_af2pred$/, "", bb)  # strip sequence + af2 suffix
      r = $2 + 0
      n[bb]++
      if (!(bb in best) || r < best[bb]) { best[bb] = r; bestseq[bb] = desc }
      if (r < C) { pass[bb]++; print desc > P }
      tot++
      if (r < C) ptot++
    }
    END {
      if (tot == 0) { print "  no numeric score lines"; exit }
      printf "  %-30s %6s %6s %7s\n", "backbone", "best", "pass", "of"
      # Selection sort by best RMSD. Deterministic order, and no reliance on
      # gawk-only asort or on an external sort that would reorder the summary.
      nk = 0
      for (k in best) keys[++nk] = k
      for (i = 1; i <= nk; i++) {
        m = i
        for (j = i + 1; j <= nk; j++)
          if (best[keys[j]] < best[keys[m]]) m = j
        t = keys[i]; keys[i] = keys[m]; keys[m] = t
      }
      for (i = 1; i <= nk; i++) {
        k = keys[i]
        printf "  %-30s %6.2f %6d %7d %s\n", k, best[k], pass[k] + 0, n[k],
               (pass[k] + 0 == 0 ? "<- undesignable, drop the backbone" : "")
      }
      printf "  --- %d of %d sequences under %.1f A\n", ptot + 0, tot, C
    }'

  if [ -s "$ROOT/gate2/$NAME.pass" ]; then
    echo "  survivors listed in gate2/$NAME.pass ($(wc -l < "$ROOT/gate2/$NAME.pass" | tr -d ' '))"
  fi
done

echo "============================================================"
echo "cut-off: binder_aligned_rmsd < $CUT A.  On IL-6, 31 of 40 passed."
echo "A backbone with zero passing sequences is not designable - drop it, do not"
echo "try to rescue it. Rank what survives by grip, never by confidence."
