#!/bin/bash
# Honest re-measure: train the Riemannian flow on corpus_v1 restricted to UNIQUE
# sequences (exact-dedup keep-list from 90 --emit-keep), so recovery is not
# inflated by duplicate chains straddling train/test. Compare its balanced
# recovery to the redundant 0.267.
#
# First (login node): python3 model/90_redundancy_audit.py ~/rfd/corpus_v1/raw \
#                        --emit-keep ~/rfd/corpus_v1/keep_unique.txt
# Then:  cd ~/rfd/scripts && sbatch model/91_train_flow_dedup.sbatch
#
# dgxanode03 excluded (enroot flake). Judge by the recovery block, not Slurm State.

#SBATCH --job-name=flow-dedup
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
TABLE=/workspace/rfd/corpus_v1/interface_features.tsv
KEEP=/workspace/rfd/corpus_v1/keep_unique.txt
SAVE=/workspace/rfd/corpus_v1/flow_sphere_dedup.pt

test -s "$HOME/rfd/corpus_v1/keep_unique.txt" || {
  echo "keep_unique.txt missing - run 90 --emit-keep first"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
python /workspace/rfd/scripts/model/86_flow_model.py \
  "'"$TABLE"'" --geometry sphere --keep-list "'"$KEEP"'" \
  --epochs 120 --d-model 128 --layers 3 --save "'"$SAVE"'"
echo "=== done ==="
'
