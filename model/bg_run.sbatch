#!/bin/bash
# Run a BoltzGen binder-design job from a spec YAML. Runs BARE in the uv venv
# (GPU works bare on this cluster via cu12.8 torch — no enroot container needed).
# Skips the ESMFold2 analysis step (needs a CUDA-13 driver we don't have); designs
# are ranked on Boltz confidence + geometry. Model weights are pre-cached in ~/.cache
# from the login node (compute nodes can't reach HuggingFace).
#
#   cd ~/rfd/scripts && git checkout claude/project-thread-ho43or && git pull \
#     && sbatch model/bg_run.sbatch
#   # knobs: SPEC=mdm2_blind/mdm2_bindany.yaml N=20 PROTO=protein-anything sbatch model/bg_run.sbatch
#
# Judge by the output tree + final_ranked_designs, never Slurm State.

#SBATCH --job-name=bgrun
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -uo pipefail
VENV=$HOME/rfd/venv_boltzgen
SPEC="${SPEC:-mdm2_blind/mdm2_bindany.yaml}"
N="${N:-20}"
PROTO="${PROTO:-protein-anything}"
# Steps override. Default = full pipeline minus ESMFold2 analysis. For a first run
# with SolubleMPNN not yet cached, use STEPS="design" (design step needs only the
# cached design checkpoints + moldir — no SolubleMPNN, no HF access at runtime).
STEPS="${STEPS:-design inverse_folding folding filtering}"
STAMP=$(date +%m%d_%H%M)
OUT="${OUT:-$HOME/rfd/outputs/bg_$(basename ${SPEC%.yaml})_${STAMP}}"

test -f "$VENV/bin/activate" || { echo "venv_boltzgen missing"; exit 1; }
test -f "$HOME/rfd/scripts/$SPEC" || { echo "spec $SPEC missing (git pull?)"; exit 1; }
source "$VENV/bin/activate"
cd "$HOME/rfd/scripts/$(dirname "$SPEC")"
SPECFILE="$(basename "$SPEC")"

echo "node: $(hostname)  job: $SLURM_JOB_ID  spec=$SPECFILE N=$N proto=$PROTO"
echo "out=$OUT"
nvidia-smi --query-gpu=name --format=csv,noheader
python -c "import torch;print('torch',torch.__version__,'cuda avail',torch.cuda.is_available())"

echo "=== 1/2 check spec ==="
boltzgen check "$SPECFILE"

echo "=== 2/2 run (steps: $STEPS) ==="
boltzgen run "$SPECFILE" --output "$OUT" --protocol "$PROTO" \
  --num_designs "$N" --steps $STEPS --reuse

echo "=== output tree ==="
ls -R "$OUT" 2>/dev/null | head -80
echo "=== done ==="
