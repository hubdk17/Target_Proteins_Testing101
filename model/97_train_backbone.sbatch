#!/bin/bash
# Step 2 run: train the SE(3) flow-matching BACKBONE generator (96) on the
# extracted backbone frames. Runs inside the enroot container (login python3
# has no torch); the GPU does the work. dgxanode03 excluded (enroot flake).
#
# The job SELF-CHECKS first: it runs 96 --smoke (seconds, catches wiring bugs)
# and aborts before training if that fails. Judge success by the loss curve in
# the .out and the saved checkpoint, never by Slurm State.
#
# First (login node) build the dataset if not done:
#   python3 model/95_extract_frames.py ~/rfd/corpus_v1/raw \
#       --keep-list ~/rfd/corpus_v1/keep_unique.txt \
#       --out ~/rfd/corpus_v1/backbones.tsv
# Then:  cd ~/rfd/scripts && sbatch model/97_train_backbone.sbatch

#SBATCH --job-name=bbflow
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
DATA=/workspace/rfd/corpus_v1/backbones.tsv
SAVE=/workspace/rfd/corpus_v1/bbflow.pt

test -s "$HOME/rfd/corpus_v1/backbones.tsv" || {
  echo "backbones.tsv missing - run 95_extract_frames.py first"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
echo "=== smoke test (abort train if it fails) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py --smoke
echo "=== training ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py \
  --data "'"$DATA"'" --epochs 300 --batch 16 \
  --d-model 128 --layers 4 --lr 3e-4 --log-every 10 \
  --w-dist 2.0 \
  --save "'"$SAVE"'"
echo "=== done ==="
test -f "'"$SAVE"'" && echo "checkpoint written: '"$SAVE"'"
'
