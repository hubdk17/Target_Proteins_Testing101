#!/bin/bash
# A BATTERY of architecturally-distinct backbone-generator variants, each a
# different MATHEMATICAL modelling of the flow's base noise / geometry, trained
# in parallel as a Slurm array on the big corpus (corpus_v2). The point is to let
# the structural-coherence scorecard (104) pick a winner, since the single-chain
# overfit micro-test is confounded (overfitting one example breaks sampling-time
# generalisation, so it can't adjudicate a prior). All share the winning colored
# TRANSLATION prior + aux chain loss; each changes ONE lever:
#
#   task 0  rotsmooth_s05 : + smooth SO(3) ROTATION prior, sigma 0.5
#           (neighbouring start-frames vary smoothly -- the rotation analog of the
#            colored translation prior; the untouched channel where SS lives)
#   task 1  rotsmooth_s10 : + smooth ROTATION prior, sigma 1.0
#           (stronger per-residue turn, nearer a real helix's ~100 deg/residue)
#   task 2  alpha3        : colored translation with alpha 3.0 (vs 2.0 baseline)
#           (more COLLECTIVE chain modes -- a stiffer, longer-correlated prior)
# Baseline for comparison is the existing bbflow_colored_v2.pt (alpha 2, iid rot).
# This is a BLIND backbone stage (no target), so the levers are geometric/SO(3),
# not chemistry; target chemistry belongs to the conditioned model, not here.
#
#   cd ~/rfd/scripts && git checkout claude/project-thread-ho43or && git pull \
#     && sbatch --array=0-2 model/103_train_battery.sbatch
# Then, once all three finish (squeue empty):
#   sbatch model/104_score_battery.sbatch

#SBATCH --job-name=bbbatt
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=05:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
DATA=/workspace/rfd/corpus_v2/backbones.tsv

NAMES=(rotsmooth_s05 rotsmooth_s10 alpha3)
EXTRA=("--prior colored --alpha 2.0 --kappa 0.1 --rot-prior smooth --rot-sigma 0.5"
       "--prior colored --alpha 2.0 --kappa 0.1 --rot-prior smooth --rot-sigma 1.0"
       "--prior colored --alpha 3.0 --kappa 0.1 --rot-prior iid")
ID=${SLURM_ARRAY_TASK_ID:-0}
NAME=${NAMES[$ID]}
CFG=${EXTRA[$ID]}
SAVE=/workspace/rfd/corpus_v1/bbflow_$NAME.pt

test -s "$HOME/rfd/corpus_v2/backbones.tsv" || {
  echo "corpus_v2/backbones.tsv missing - run 95_extract_frames.py first"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID  array task: $ID -> $NAME"
echo "config: $CFG"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
echo "=== smoke test (abort train if it fails) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py --smoke
echo "=== training '"$NAME"' (corpus_v2, 40 epochs) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py \
  --data "'"$DATA"'" --epochs 40 --batch 16 \
  --d-model 128 --layers 4 --lr 3e-4 --log-every 2 \
  --w-dist 2.0 --w-dgram 0.0 '"$CFG"' \
  --save "'"$SAVE"'"
echo "=== done ==="
test -f "'"$SAVE"'" && echo "checkpoint written: '"$SAVE"'"
'
