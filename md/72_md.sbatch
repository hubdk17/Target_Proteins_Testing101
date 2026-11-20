#!/bin/bash
#SBATCH --job-name=mdm2_md
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=32G
# NOTE: pass --gres=gpu:1 and --nodelist=dgxanode04 (A100) on the submit line.
# dgxanode03's enroot disk is intermittently broken; dgxhnode6 is H100 (use
# the A100 first, it is the safer prebuilt-binary target).
IMG=/home/hpc_apps/nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh
SYS="$1"      # directory holding sys.prmtop / sys.inpcrd
srun --container-image=$IMG --container-mounts=/home/dakshcsed:/workspace \
  bash -c '
    set -euo pipefail
    export AMBERHOME=/workspace/amber24
    source "$AMBERHOME/amber.sh"
    export LD_LIBRARY_PATH=/workspace/usr/lib:/workspace/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}
    cd '"$SYS"'
    P=/workspace/rfd/scripts/md
    pmemd.cuda -O -i $P/min.in  -p sys.prmtop -c sys.inpcrd -ref sys.inpcrd -o min.out  -r min.rst  -x min.nc
    pmemd.cuda -O -i $P/heat.in -p sys.prmtop -c min.rst   -ref min.rst   -o heat.out -r heat.rst -x heat.nc
    pmemd.cuda -O -i $P/equil.in -p sys.prmtop -c heat.rst -ref heat.rst  -o equil.out -r equil.rst -x equil.nc
    pmemd.cuda -O -i $P/prod.in -p sys.prmtop -c equil.rst               -o prod.out -r prod.rst -x prod.nc
  '
