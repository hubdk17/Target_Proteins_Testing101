#!/bin/bash
# Test the SAMPLER FIX on the EXISTING trained checkpoints -- no retraining.
# CPU diagnosis (diag_sampler/diag_fix) showed the sampled chain inflates because
# the OT velocity (px1-x)/(1-t) blows up near t=1 and compounds off-manifold drift:
# the model predicts a well-connected chain on-path (CA-CA ~3.8) but free
# integration walks off-manifold, and MORE steps make it WORSE (not a discretisation
# error). Two inference-only knobs tame it: --vclamp (cap the velocity) and --t-stop
# (end before the unstable late-t region). This sweeps them on the two best models
# (alpha3 = best balance, rotsmooth_s05 = best raw C-N) and coherence-scores each,
# so the real-model scorecard decides whether the fix lifts connectivity.
#
#   cd ~/rfd/scripts && git checkout claude/project-thread-ho43or && git pull \
#     && sbatch model/105_sampler_fix_sweep.sbatch
# Read the .out: per (model,config) "coherence scorecard"; compare C-N / CA-CA /
# clashes against each model's no-clamp baseline row.

#SBATCH --job-name=bbfix
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh

echo "node: $(hostname)  job: $SLURM_JOB_ID"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -uo pipefail
cd /workspace/rfd/scripts
nvidia-smi --query-gpu=name --format=csv,noheader
# config tag -> extra sampler flags
declare -A CFG=(
  [baseline]=""
  [vclamp1]="--vclamp 1.0"
  [vclamp05]="--vclamp 0.5"
  [tstop07]="--t-stop 0.7"
)
for m in bbflow_alpha3 bbflow_rotsmooth_s05; do
  CK=/workspace/rfd/corpus_v1/$m.pt
  [ -f "$CK" ] || { echo "=== $m: no checkpoint, skipping ==="; continue; }
  for tag in baseline vclamp1 vclamp05 tstop07; do
    OUT=/workspace/rfd/corpus_v1/fix_${m}_${tag}
    echo "================= $m | $tag ${CFG[$tag]} ================="
    python model/98_sample_backbone.py --ckpt "$CK" --n 20 --length 80 --steps 200 \
      ${CFG[$tag]} --out-dir "$OUT" 2>&1 | tail -1
    python model/102_coherence.py "$OUT"
  done
done
echo "=== done ==="
'
