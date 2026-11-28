#!/bin/bash
# Step 2 variant: SAME as 97 (d128 L4, aux connectivity loss) but draws the flow's
# base noise from the CHAIN-CORRELATED (colored) prior instead of iid N(0,I).
# Differs from the running bbflow job (187213) by the PRIOR ONLY, so an improvement
# is attributable to it. The prior starts the sampler from a connected chain at the
# right 3.8 A spacing (CPU-validated: overfit C-N 15.8 -> 9.5 A, training loss
# collapses). Distinct job-name + save so it NEVER clobbers the running bbflow.pt.
# Built on the COVALENT chain graph only -> no fold leakage. dgxanode03 excluded.
#
# Run in parallel:  cd ~/rfd/scripts && git pull && sbatch model/97c_train_backbone_colored.sbatch
# After it finishes, sample + check C-N (98 reads the prior from the ckpt):
#   srun ... 98_sample_backbone.py --ckpt /workspace/rfd/corpus_v1/bbflow_colored.pt ...

#SBATCH --job-name=bbflowCOL
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
SAVE=/workspace/rfd/corpus_v1/bbflow_colored.pt

test -s "$HOME/rfd/corpus_v1/backbones.tsv" || {
  echo "backbones.tsv missing - run 95_extract_frames.py first"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID  (colored prior)"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
echo "=== smoke test (abort train if it fails) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py --smoke
echo "=== training (colored prior, alpha=2 kappa=0.1) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py \
  --data "'"$DATA"'" --epochs 300 --batch 16 \
  --d-model 128 --layers 4 --lr 3e-4 --log-every 10 \
  --w-dist 2.0 --prior colored --alpha 2.0 --kappa 0.1 \
  --save "'"$SAVE"'"
echo "=== done ==="
test -f "'"$SAVE"'" && echo "checkpoint written: '"$SAVE"'"
'
