#!/bin/bash
# Sample all three trained backbone checkpoints and score STRUCTURAL COHERENCE of
# each (not just C-N): bond geometry, steric clashes, secondary-structure content,
# compactness. One job -> a head-to-head scorecard for aux-loss vs capacity vs
# colored-prior. 98 auto-reads each checkpoint's prior. dgxanode03 excluded.
#
#   cd ~/rfd/scripts && git pull && sbatch model/101_sample_all.sbatch
#   (cancel the earlier interactive srun first if still queued: scancel 188637)
# Read the .out: per-model "coherence scorecard" block + COHERENT-ish verdict.

#SBATCH --job-name=bbsample
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh

echo "node: $(hostname)  job: $SLURM_JOB_ID"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -uo pipefail
cd /workspace/rfd/scripts
nvidia-smi --query-gpu=name --format=csv,noheader
for c in bbflow bbflow_d256 bbflow_colored; do
  CK=/workspace/rfd/corpus_v1/$c.pt
  [ -f "$CK" ] || { echo "=== $c: no checkpoint, skipping ==="; continue; }
  echo "================= $c ================="
  python model/98_sample_backbone.py --ckpt "$CK" --n 20 --length 80 --steps 200 \
    --out-dir /workspace/rfd/corpus_v1/samp_$c 2>&1 | tail -2
  python model/102_coherence.py /workspace/rfd/corpus_v1/samp_$c
done
echo "=== done ==="
'
