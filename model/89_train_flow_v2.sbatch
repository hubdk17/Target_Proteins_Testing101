#!/bin/bash
# Train the WINNING config -- Riemannian (sphere) flow with ENRICHED features --
# on corpus_v2 (the big ~20k-complex corpus). Concats the featurize shards first
# if interface_features.tsv isn't there yet, then trains on a GPU in the enroot
# container (login python3 has no torch/numpy).
#
# Submit from the scripts repo root, AFTER the 88 featurize array has finished:
#   cd ~/rfd/scripts && sbatch model/89_train_flow_v2.sbatch
#
# dgxanode03 excluded (enroot flake). Judge success by the recovery block in the
# .out, never the Slurm State.

#SBATCH --job-name=flow-v2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
CORP="$HOME/rfd/corpus_v2"
TSV="$CORP/interface_features.tsv"

# Concat the shards if the combined table isn't built yet (keep one header).
if [ ! -s "$TSV" ]; then
  shopt -s nullglob
  shards=("$CORP"/if_*.tsv)
  test "${#shards[@]}" -gt 0 || { echo "no if_*.tsv shards in $CORP - run 88 first"; exit 1; }
  echo "concatenating ${#shards[@]} shards -> $TSV"
  head -1 "${shards[0]}" > "$TSV"
  for f in "${shards[@]}"; do tail -n +2 "$f"; done >> "$TSV"
fi
echo "training table: $TSV ($(wc -l < "$TSV") rows)"

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
SAVE=/workspace/rfd/corpus_v2/flow_sphere_enriched.pt

echo "node: $(hostname)  job: $SLURM_JOB_ID"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
python /workspace/rfd/scripts/model/86_flow_model.py \
  /workspace/rfd/corpus_v2/interface_features.tsv \
  --geometry sphere --enriched --epochs 150 --d-model 128 --layers 3 \
  --save "'"$SAVE"'"
echo "=== done ==="
test -f "'"$SAVE"'" && echo "model written: '"$SAVE"'"
'
