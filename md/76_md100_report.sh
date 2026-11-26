#!/bin/bash
# Read a 100 ns run (all prod_??.nc chunks) and print a per-10 ns-window
# summary, so drift ONSET is visible instead of only a start/end number.
# Run INSIDE the nvhpc container, no GPU needed:
#   srun --cpus-per-task=4 --mem=16G --time=00:40:00 \
#     --container-image=/home/hpc_apps/nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh \
#     --container-mounts=/home/dakshcsed:/workspace \
#     bash /workspace/rfd/scripts_md/md/76_md100_report.sh p53_mdm2
#
# args: <run-name under runs100/> [binder_len=13] [anchor offsets, 1-based in
#       the binder, default "3 7 10" = p53 Phe19 / Trp23 / Leu26]
set -uo pipefail
S="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$S/md/_amber_env.sh"
cd "$S"

NAME="$1"; BLEN="${2:-13}"; ANCH="${3:-3 7 10}"
D="runs100/$NAME"
[ -d "$D" ] || { echo "no such run dir: $D"; exit 1; }

CH=$(ls "$D"/prod_??.nc 2>/dev/null | wc -l)
if [ "$CH" -eq 0 ]; then
  echo "NO prod chunks in $D -- build or GPU stage failed. Tails:"
  tail -6 "$D"/min.out "$D"/prod_01.out 2>/dev/null
  grep -i "D-RESIDUE\|NO CHAIN BREAK\|no kernel image\|Error" "$D"/*.out "$D"/tleap.log 2>/dev/null | head -5
  exit 1
fi
echo "=== $NAME : $CH chunk(s) = $((CH*10)) ns ==="

NP=$(awk 'substr($0,1,4)=="ATOM"{rn=substr($0,18,3);gsub(/ /,"",rn);
          if(rn!="WAT"&&rn!="HOH"&&rn!="Na+"&&rn!="Cl-"&&rn!="NA"&&rn!="CL"){seen[substr($0,23,4)]=1}}
          END{print length(seen)}' "$D/sys_built.pdb")
BST=$((NP-BLEN+1)); TEND=$((NP-BLEN))
echo "system: $NP protein residues -- target :1-$TEND, binder :$BST-$NP (len $BLEN)"

{
  echo "parm $D/sys.prmtop"
  for f in "$D"/prod_??.nc; do echo "trajin $f"; done
  echo "autoimage"
  echo "rms fitTarget :1-$TEND@CA first"
  echo "rms binderHeavy :$BST-$NP&!@H= nofit out $D/binder_rmsd.dat"
  echo "nativecontacts :$BST-$NP :1-$TEND distance 4.5 out $D/contacts.dat"
  for o in $ANCH; do
    r=$((BST+o-1))
    echo "nativecontacts :$r :1-$TEND distance 4.5 out $D/anchor_$r.dat"
  done
} > "$D/report.in"

cpptraj -i "$D/report.in" > "$D/report.log" 2>&1 || { echo "cpptraj FAILED"; tail -8 "$D/report.log"; exit 1; }

win () {  # per-10 ns-window mean of a 2-col cpptraj .dat
  awk -v nw="$2" 'NR>1 && $2!=""{n++; v[n]=$2}
    END{ if(n==0){print "  no data"; exit}
         w=n/nw; printf "   "; for(i=1;i<=nw;i++){s=0;k=0;
           for(j=int((i-1)*w)+1;j<=int(i*w);j++){s+=v[j];k++}
           printf "%7.1f", (k? s/k : 0)} printf "\n"}' "$1"
}
printf "   windows (ns):"; for i in $(seq 1 "$CH"); do printf "%7d" $((i*10)); done; echo
echo " binder heavy-atom RMSD vs start (A), per window:"; win "$D/binder_rmsd.dat" "$CH"
echo " interface contacts <4.5 A, per window:";            win "$D/contacts.dat" "$CH"
for o in $ANCH; do
  r=$((BST+o-1))
  rn=$(awk -v r="$r" 'substr($0,1,4)=="ATOM" && substr($0,23,4)+0==r {print substr($0,18,3); exit}' "$D/sys_built.pdb")
  echo " anchor $rn (binder pos $o, built :$r) contacts, per window:"; win "$D/anchor_$r.dat" "$CH"
done
echo "-----------------------------------------------------------------"
echo "Read: RMSD flat + contacts flat = stayed docked. A window where an"
echo "anchor's contacts fall to ~0 and stay there is the anchor letting go;"
echo "that is interface-local and is NOT visible in the total contact count."
echo "MD at 100 ns tests ASSEMBLY, not affinity: 'fell off' is strong"
echo "evidence, 'stayed on' is weak."
