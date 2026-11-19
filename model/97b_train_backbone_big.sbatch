#!/bin/bash
# Step 2 variant: SAME as 97 (aux connectivity loss, --w-dist 2.0) but a BIGGER
# model (d256, L6 vs d128, L4). Differs from the running bbflow job (187213) by
# capacity ONLY, so if C-N / loss improve we can attribute it to model size and
# test the "406k params is underpowered" hypothesis. Distinct job-name + save so
# it NEVER clobbers the running bbflow.pt. dgxanode03 excluded (enroot flake).
#
# Run in parallel with 97:  cd ~/rfd/scripts && git pull && sbatch model/97b_train_backbone_big.sbatch
# After it finishes, sample + check C-N the same way:
#   srun ... 98_sample_backbone.py --ckpt /workspace/rfd/corpus_v1/bbflow_d256.pt ...

#SBATCH --job-name=bbflowBIG
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
SAVE=/workspace/rfd/corpus_v1/bbflow_d256.pt

test -s "$HOME/rfd/corpus_v1/backbones.tsv" || {
  echo "backbones.tsv missing - run 95_extract_frames.py first"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID  (d256 L6)"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
echo "=== smoke test (abort train if it fails) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py --smoke
echo "=== training (d256 L6) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py \
  --data "'"$DATA"'" --epochs 300 --batch 16 \
  --d-model 256 --layers 6 --lr 3e-4 --log-every 10 \
  --w-dist 2.0 \
  --save "'"$SAVE"'"
echo "=== done ==="
test -f "'"$SAVE"'" && echo "checkpoint written: '"$SAVE"'"
'
