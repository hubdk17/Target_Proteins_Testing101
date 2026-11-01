#!/bin/bash
# STEP 1 of 4 -- the cheapest possible job. ~30 seconds of GPU time.
#
# Proves four things before you spend anything: you can submit to defq, the
# Enroot image loads, your home shows up at /workspace, and torch sees the GPU.
# No RFdiffusion, no downloads, nothing to go wrong.
#
#   module load slurm
#   mkdir -p ~/rfd/logs && cd ~/rfd
#   sbatch 10_hello_gpu.sbatch
#   cat logs/hello_*.out

#SBATCH --job-name=rfd-hello
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --gres=gpu:1
#SBATCH --time=00:05:00
#SBATCH --output=logs/hello_%j.out
#SBATCH --error=logs/hello_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh

echo "node: $(hostname)   job: $SLURM_JOB_ID"

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -c '
       echo "--- inside container ---"
       python --version
       nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
       python -c "import torch; print(\"torch\", torch.__version__, \"| cuda\", torch.version.cuda, \"| GPU visible:\", torch.cuda.is_available())"
       echo "--- home mounted at /workspace? ---"
       ls /workspace | head
       echo "--- can the compute node reach PyPI? ---"
       timeout 60 pip download --no-deps -q hydra-core -d /tmp/pipcheck \
         && echo "PyPI OK" || echo "PyPI UNREACHABLE from compute node -- tell Claude, it changes step 2"
     '
