#!/bin/bash
# STEP 4 -- generate a backbone. Deliberately small: 1 design, 50 residues,
# 20 diffusion timesteps (RFdiffusion's own default is 50), 20-minute cap.
#
#   sbatch 30_rfd.sbatch                       # unconditional smoke test
#   MODE=binder sbatch 30_rfd.sbatch           # trypsin target, no hotspots
#
# Scale up once it works, and raise --time to match. Cost goes as
# designs x timesteps x residues^2:
#   NUM_DESIGNS=10 LENGTH=100 DIFFUSER_T=50 sbatch 30_rfd.sbatch

#SBATCH --job-name=rfd-run
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --output=logs/rfd_%j.out
#SBATCH --error=logs/rfd_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh

MODE="${MODE:-unconditional}"
NUM_DESIGNS="${NUM_DESIGNS:-1}"
LENGTH="${LENGTH:-50}"
DIFFUSER_T="${DIFFUSER_T:-20}"
BINDER_LEN="${BINDER_LEN:-70-100}"
# Confirmed from the structure: 2PTN chain A is numbered 16-245 (223 residues).
TARGET_CONTIG="${TARGET_CONTIG:-A16-245}"

echo "node: $(hostname)  job: $SLURM_JOB_ID  mode: $MODE"
test -f "$HOME/rfd/venv/.complete" || { echo "run 20_buildenv.sbatch first"; exit 1; }

cat > "$HOME/rfd/_run_inner.sh" <<'INNER'
set -euo pipefail
ROOT=/workspace/rfd
source "$ROOT/venv/bin/activate"
export DGLBACKEND=pytorch     # dgl warns and guesses otherwise
cd "$ROOT/RFdiffusion"
mkdir -p "$ROOT/outputs"

nvidia-smi --query-gpu=name --format=csv,noheader

case "$MODE" in
  unconditional)
    echo "=== unconditional: ${LENGTH} aa, ${NUM_DESIGNS} design(s), T=${DIFFUSER_T} ==="
    # No target at all. This proves the pipeline runs; the output has nothing
    # to do with trypsin.
    python scripts/run_inference.py \
      inference.output_prefix="$ROOT/outputs/uncond_${LENGTH}aa" \
      inference.model_directory_path="$ROOT/models" \
      inference.num_designs="$NUM_DESIGNS" \
      diffuser.T="$DIFFUSER_T" \
      "contigmap.contigs=[${LENGTH}-${LENGTH}]"
    ;;
  binder)
    echo "=== binder vs trypsin, no hotspots: ${NUM_DESIGNS} design(s), T=${DIFFUSER_T} ==="
    # '/0 ' with the trailing space is the chain break between target and the
    # designed binder. Complex_base_ckpt.pt is the protein-protein model.
    python scripts/run_inference.py \
      inference.output_prefix="$ROOT/outputs/trypsin_binder_nohotspot" \
      inference.model_directory_path="$ROOT/models" \
      inference.input_pdb="$ROOT/inputs/trypsin_clean.pdb" \
      inference.ckpt_override_path="$ROOT/models/Complex_base_ckpt.pt" \
      inference.num_designs="$NUM_DESIGNS" \
      diffuser.T="$DIFFUSER_T" \
      "contigmap.contigs=[${TARGET_CONTIG}/0 ${BINDER_LEN}]"
    ;;
  *) echo "MODE must be unconditional or binder"; exit 2 ;;
esac

echo "=== done ==="
ls -lh "$ROOT/outputs"
INNER

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL,MODE="$MODE",NUM_DESIGNS="$NUM_DESIGNS",LENGTH="$LENGTH",DIFFUSER_T="$DIFFUSER_T",BINDER_LEN="$BINDER_LEN",TARGET_CONTIG="$TARGET_CONTIG" \
     bash /workspace/rfd/_run_inner.sh
