#!/bin/bash
# OUR FIRST NEURAL MODEL: train 84_train_mlp.py on a GPU inside the enroot
# container. This is the capacity test above the linear 82 baseline, on the
# same blind corpus and the same by-structure split.
#
# Submit (from the scripts repo root, after `module load slurm`):
#   cd ~/rfd/scripts && sbatch model/85_train_mlp.sbatch
#
# The container and mount pattern are copied from the AF2 job that actually ran
# (40_gate2_selfconsistency.sbatch, jobs 174676/174717), so nothing here needs
# building: the image ships torch + numpy. The login python3 has NO numpy, which
# is why this must run inside the container, never on the login node.
#
# dgxanode03's local enroot disk fills intermittently and kills a job in <1 min
# with an empty .out -- so exclude it. node04 (A100) and node06 (H100) are fine
# for this; the model is tiny.
#
# Judge success by the .out file: it must print a BALANCED acc line per arm.
# Never trust the Slurm State -- jobs here have reported FAILED on success.

#SBATCH --job-name=mlp-cond
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
TABLE=/workspace/rfd/corpus_v1/interface_features.tsv
SAVE=/workspace/rfd/corpus_v1/mlp_cond_model

# Fail before the GPU is allocated if the training table is missing.
test -s "$HOME/rfd/corpus_v1/interface_features.tsv" || {
  echo "interface_features.tsv missing - run 81_featurize_corpus.py first"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID"

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
python /workspace/rfd/scripts/model/84_train_mlp.py \
  "'"$TABLE"'" --compare --epochs 60 --hidden 128,64 \
  --save "'"$SAVE"'"
echo "=== done ==="
test -f "'"$SAVE"'.pt" && echo "model written: '"$SAVE"'.pt"
'
