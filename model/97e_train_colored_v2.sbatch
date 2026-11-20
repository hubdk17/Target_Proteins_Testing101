#!/bin/bash
# Retrain the WINNING config (colored prior + aux loss) on the 15x-bigger backbone
# corpus (corpus_v2: 18,675 monomers / 2.3M frames vs corpus_v1's 1,248 / 154k).
# Data is the one lever changed vs 97c, to test whether the missing secondary
# structure (all three v1 models were ~75% loop) is a data-starvation problem.
# Fewer epochs because each epoch now sees 15x more chains. dgram stays OFF
# (--w-dgram 0; CPU-validated as unhelpful). dgxanode03 excluded.
#
#   cd ~/rfd/scripts && git checkout claude/project-thread-ho43or && git pull && sbatch model/97e_train_colored_v2.sbatch
# Then sample + coherence-score it (98 reads the prior from the ckpt):
#   ...98_sample_backbone.py --ckpt /workspace/rfd/corpus_v1/bbflow_colored_v2.pt ...
#   ...102_coherence.py <out-dir>

#SBATCH --job-name=bbflowCV2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=05:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
DATA=/workspace/rfd/corpus_v2/backbones.tsv
SAVE=/workspace/rfd/corpus_v1/bbflow_colored_v2.pt

test -s "$HOME/rfd/corpus_v2/backbones.tsv" || {
  echo "corpus_v2/backbones.tsv missing - run 95_extract_frames.py first"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID  (colored prior, corpus_v2)"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -euo pipefail
nvidia-smi --query-gpu=name --format=csv,noheader
echo "=== smoke test (abort train if it fails) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py --smoke
echo "=== training (colored prior, corpus_v2, 40 epochs) ==="
python /workspace/rfd/scripts/model/96_flow_backbone.py \
  --data "'"$DATA"'" --epochs 40 --batch 16 \
  --d-model 128 --layers 4 --lr 3e-4 --log-every 2 \
  --w-dist 2.0 --prior colored --alpha 2.0 --kappa 0.1 --w-dgram 0.0 \
  --save "'"$SAVE"'"
echo "=== done ==="
test -f "'"$SAVE"'" && echo "checkpoint written: '"$SAVE"'"
'
