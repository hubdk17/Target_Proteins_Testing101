#!/bin/bash
#SBATCH --job-name=mdm2_md
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=32G
# Submit one per arm, e.g.:
#   sbatch --gres=gpu:1 --nodelist=dgxanode04 --time=12:00:00 \
#          --output=logs/md_%x_%j.out --error=logs/md_%x_%j.err \
#          --job-name=md_armA 73_md_all.sbatch mdm2_blind/armA_design.pdb runs/armA
# Build (tleap, CPU, ~2 min) + chirality gate + min/heat/equil/prod, in one
# GPU job. The minimisation is itself the pmemd.cuda / GPU-load test: if the
# binary or CUDA is wrong we find out in the first minute, not after a day.
PDB="$1"      # complex PDB: MDM2 chain A + binder chain B
OUT="$2"      # output directory for this arm
IMG=/home/hpc_apps/nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh
srun --container-image=$IMG --container-mounts=/home/dakshcsed:/workspace \
  bash -c '
    set -euo pipefail
    source /workspace/rfd/scripts/md/_amber_env.sh
    cd /workspace/rfd/scripts
    bash md/70_md_build.sh "/workspace/rfd/scripts/'"$PDB"'" "/workspace/rfd/scripts/'"$OUT"'"
    cd "/workspace/rfd/scripts/'"$OUT"'"
    P=/workspace/rfd/scripts/md
    echo "=== minimisation (also the GPU-load test) ==="
    pmemd.cuda -O -i $P/min.in   -p sys.prmtop -c sys.inpcrd -ref sys.inpcrd -o min.out   -r min.rst   -x min.nc
    echo "=== heat ==="
    pmemd.cuda -O -i $P/heat.in  -p sys.prmtop -c min.rst    -ref min.rst    -o heat.out  -r heat.rst  -x heat.nc
    echo "=== equil ==="
    pmemd.cuda -O -i $P/equil.in -p sys.prmtop -c heat.rst   -ref heat.rst   -o equil.out -r equil.rst -x equil.nc
    echo "=== production ==="
    pmemd.cuda -O -i $P/prod.in  -p sys.prmtop -c equil.rst                   -o prod.out  -r prod.rst  -x prod.nc
    echo "=== DONE ==="
  '
