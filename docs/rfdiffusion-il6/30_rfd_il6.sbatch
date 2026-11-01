#!/bin/bash
# IL-6 binder, no hotspots, no mutation. Deliberately tiny: 1 design, T=20.
# Submit from ~/rfd  (the log paths below are relative to where you sbatch).
#
#   cd ~/rfd && sbatch 30_rfd_il6.sbatch
#
# Scale up later with env vars, and raise --time to match:
#   NUM_DESIGNS=10 DIFFUSER_T=50 sbatch 30_rfd_il6.sbatch

#SBATCH --job-name=rfd-il6
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --output=logs/rfd_il6_%j.out
#SBATCH --error=logs/rfd_il6_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh

NUM_DESIGNS="${NUM_DESIGNS:-1}"
DIFFUSER_T="${DIFFUSER_T:-20}"
BINDER_LEN="${BINDER_LEN:-70-100}"
# Confirmed from 1ALU on 2026-10-21: chain A, residues 19-184, with 52-60
# genuinely absent from the coordinates. Hence two segments, not one range.
TARGET_CONTIG="${TARGET_CONTIG:-A19-51/A61-184}"

echo "node: $(hostname)  job: $SLURM_JOB_ID"
echo "contig: [${TARGET_CONTIG}/0 ${BINDER_LEN}]"
test -f "$HOME/rfd/venv/.complete"        || { echo "venv missing - run 20_buildenv.sbatch"; exit 1; }
test -f "$HOME/rfd/inputs/il6_clean.pdb"  || { echo "il6_clean.pdb missing - redo step 1c"; exit 1; }

cat > "$HOME/rfd/_run_il6.sh" <<'INNER'
set -euo pipefail
ROOT=/workspace/rfd
source "$ROOT/venv/bin/activate"
export DGLBACKEND=pytorch
cd "$ROOT/RFdiffusion"
mkdir -p "$ROOT/outputs"
nvidia-smi --query-gpu=name --format=csv,noheader

python scripts/run_inference.py \
  inference.output_prefix="$ROOT/outputs/il6_binder_nohotspot" \
  inference.model_directory_path="$ROOT/models" \
  inference.input_pdb="$ROOT/inputs/il6_clean.pdb" \
  inference.ckpt_override_path="$ROOT/models/Complex_base_ckpt.pt" \
  inference.num_designs="$NUM_DESIGNS" \
  diffuser.T="$DIFFUSER_T" \
  "contigmap.contigs=[${TARGET_CONTIG}/0 ${BINDER_LEN}]"

echo "=== done ==="
ls -lh "$ROOT/outputs" | tail -5
INNER

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL,NUM_DESIGNS="$NUM_DESIGNS",DIFFUSER_T="$DIFFUSER_T",BINDER_LEN="$BINDER_LEN",TARGET_CONTIG="$TARGET_CONTIG" \
     bash /workspace/rfd/_run_il6.sh
