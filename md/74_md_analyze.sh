#!/bin/bash
# Analyse the MD arms ON the cluster and print a compact summary (no raw
# trajectory leaves the cluster; scp off bcmnode1 does not work anyway).
# Run INSIDE the nvhpc container, no GPU needed:
#   srun --cpus-per-task=4 --mem=16G --time=00:20:00 \
#     --container-image=/home/hpc_apps/nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh \
#     --container-mounts=/home/dakshcsed:/workspace \
#     bash /workspace/rfd/scripts/md/74_md_analyze.sh
set -uo pipefail
source /workspace/rfd/scripts/md/_amber_env.sh
cd /workspace/rfd/scripts

summ () {  # summarise a cpptraj 2-col .dat: first, last, mean of last 20%
  awk 'NR>1 && $2!=""{n++; v[n]=$2}
       END{ if(n==0){print "no data"; exit}
            s=0; k=0; for(i=int(n*0.8)+1;i<=n;i++){s+=v[i];k++}
            printf "first=%.2f  last=%.2f  mean_last20%%=%.2f  (frames=%d)\n", v[1], v[n], s/k, n }' "$1"
}

# arm dirs: pass names as args (they map to runs/<name>); default to the
# original three. e.g.  bash md/74_md_analyze.sh armA armB ref design_01 L15_03
ARMS="${@:-armA armB ref}"
for arm in $ARMS; do
  D=runs/$arm
  echo "===================== $arm ====================="
  if [ ! -f "$D/prod.nc" ]; then
    echo "  NO prod.nc -- job aborted or crashed. tail of prod.out / tleap.log:"
    tail -6 "$D/prod.out" 2>/dev/null || echo "    (no prod.out)"
    grep -i "D-residue\|Error\|error while loading" "$D"/*.out "$D"/tleap.log 2>/dev/null | head -4
    continue
  fi
  # protein residue count (exclude water/ions); binder = last 13, target = rest
  NP=$(awk 'substr($0,1,4)=="ATOM"{rn=substr($0,18,3);gsub(/ /,"",rn);
            if(rn!="WAT"&&rn!="HOH"&&rn!="Na+"&&rn!="Cl-"){seen[substr($0,23,4)]=1}}
            END{print length(seen)}' "$D/sys_built.pdb")
  BST=$((NP-12)); TGT_END=$((NP-13))
  echo "  system: $NP protein residues  (target :1-$TGT_END, binder :$BST-$NP)"
  cat > "$D/analyze.in" <<CPP
parm $D/sys.prmtop
trajin $D/prod.nc
autoimage
rms fitTarget :1-$TGT_END@CA first
rms binderHeavy :$BST-$NP&!@H= nofit out $D/binder_rmsd.dat
nativecontacts :$BST-$NP :1-$TGT_END distance 4.5 out $D/contacts.dat
CPP
  cpptraj -i "$D/analyze.in" >"$D/cpptraj.log" 2>&1 || { echo "  cpptraj FAILED"; tail -5 "$D/cpptraj.log"; continue; }
  printf "  binder heavy-atom RMSD vs start (A): "; summ "$D/binder_rmsd.dat"
  printf "  interface contacts (<4.5 A)        : "; summ "$D/contacts.dat"
done
echo "================================================"
echo "Read: binder RMSD staying low + contacts holding = stays docked."
echo "RMSD climbing + contacts falling = binder drifted/peeled off."
