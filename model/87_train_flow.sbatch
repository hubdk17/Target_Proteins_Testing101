#!/bin/bash
# OUR FLOW MODEL: train 86_flow_model.py on a GPU inside the enroot container.
# Generative flow matching over interface sequences, conditioned on target
# chemistry, on the blind corpus (MDM2/MDMX/p53 held out).
#
# Submit (from the scripts repo root, after `module load slurm`):
#   cd ~/rfd/scripts && sbatch model/87_train_flow.sbatch
#
# Runs on the EXISTING interface_features.tsv (base features) so it trains
# tonight with nothing to rebuild. After you re-featurize with the enriched 81,
# add the enriched columns by editing the python line below to add --enriched.
#
# Container/mount pattern copied from the AF2 job that actually ran
# (40_gate2_selfconsistency.sbatch) -- the image ships torch + numpy; the login
# python3 has neither, which is why this runs inside the container.
#
# dgxanode03's local enroot disk fills intermittently (empty .out, FAILED <1min)
# -> excluded. Judge success by the recovery block in the .out, never by the
# Slurm State: jobs here have reported FAILED on success.

#SBATCH --job-name=flow-seq
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
SAVE=/workspace/rfd/corpus_v1/flow_seq_model.pt

test -s "$HOME/rfd/corpus_v1/interface_features.tsv" || {
  echo "interface_features.tsv missing - run 81_featurize_corpus.py first"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID"

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
python /workspace/rfd/scripts/model/86_flow_model.py \
  "'"$TABLE"'" --epochs 120 --d-model 128 --layers 3 \
  --save "'"$SAVE"'"
echo "=== done ==="
test -f "'"$SAVE"'" && echo "model written: '"$SAVE"'"
'
