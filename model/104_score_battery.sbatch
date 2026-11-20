#!/bin/bash
# Sample + structural-coherence-score every battery variant (103) head-to-head
# against the baseline bbflow_colored_v2.pt. 98 auto-reads each checkpoint's prior
# AND rot_prior, so no per-model flags here. One .out = the full scorecard table.
#
#   cd ~/rfd/scripts && git checkout claude/project-thread-ho43or && git pull \
#     && sbatch model/104_score_battery.sbatch
# Read the .out: per-model "coherence scorecard" + COHERENT-ish verdict.

#SBATCH --job-name=bbbattscore
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
for c in bbflow_colored_v2 bbflow_rotsmooth_s05 bbflow_rotsmooth_s10 bbflow_alpha3; do
  CK=/workspace/rfd/corpus_v1/$c.pt
  [ -f "$CK" ] || { echo "=== $c: no checkpoint, skipping ==="; continue; }
  echo "================= $c ================="
  python model/98_sample_backbone.py --ckpt "$CK" --n 20 --length 80 --steps 200 \
    --out-dir /workspace/rfd/corpus_v1/samp_$c 2>&1 | tail -2
  python model/102_coherence.py /workspace/rfd/corpus_v1/samp_$c
done
echo "=== done ==="
'
