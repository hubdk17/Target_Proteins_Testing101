#!/bin/bash
# GATE 2, part 2: predict each threaded sequence ALONE and measure how far it
# lands from the backbone it was designed for. One array task per target.
#
# Method: /mnt/project-files/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md
# Keep sequences under 2 A. Never rank by pLDDT - on IL-6 the two sequences that
# folded into rods had the HIGHEST confidence (93-95) and the six correct folds
# sat at 75-79.
#
#   cd ~/rfd && sbatch --array=1 40_gate2_selfconsistency.sbatch      # pd1 only
#   cd ~/rfd && sbatch --array=1-6 40_gate2_selfconsistency.sbatch    # batch 1
#
# The environment lines are copied from the jobs that actually ran today
# (174676, 174717), so nothing here needs building: venv_af2 and the two one-time
# patches are already in place on this account. Do NOT rebuild them.
#
# If a run ever dies on a `Bio.Data.SCOPData` import or a `take_current_trace`
# AttributeError, the tree was rebuilt and lost those two patches:
#   - mmcif_parsing.py must use PDBData.protein_letters_3to1_extended
#   - venv_af2 must hold dm-haiku==0.0.13 (0.0.14+ breaks jax 0.4.20)

#SBATCH --job-name=sc-b1
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
TABLE="${TABLE:-$HOME/rfd/inputs/targets.tsv}"
TAG="${TAG:-b1}"

test -s "$TABLE" || { echo "$TABLE missing - run 02_prep_targets.py"; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:-1}"
NAME=$(sed -n "${IDX}p" "$TABLE" | cut -f1)
test -n "$NAME" || { echo "no line $IDX in $TABLE"; exit 1; }

INDIR="$HOME/rfd/sc_in/$NAME"
NIN=$(ls "$INDIR"/*.pdb 2>/dev/null | wc -l)
test "$NIN" -gt 0 || { echo "no threaded inputs in $INDIR - run thread_multi.py"; exit 1; }

# Fail here rather than after the GPU has been allocated and the container started.
test -f "$HOME/rfd/venv_af2/bin/activate" || { echo "venv_af2 missing - do not rebuild it blindly, check why"; exit 1; }
test -f "$HOME/rfd/dl_binder_design/af2_initial_guess/predict.py" || { echo "predict.py missing"; exit 1; }

# Fresh score filename every run: the score file is append-only, and a reused
# name mixes two runs together with no warning.
STAMP=$(date +%m%d_%H%M)
SCORE="sc_${NAME}_${TAG}_${STAMP}.sc"

echo "node: $(hostname)  job: $SLURM_JOB_ID  task: $IDX"
echo "target: $NAME   inputs: $NIN   score file: $SCORE"
echo "expect ~25 s per structure (measured: job 174717 did 40 in 15:58), so roughly $(( NIN * 25 / 60 )) min"

INNER="$HOME/rfd/_run_sc_${IDX}.sh"

cat > "$INNER" <<'EOS'
set -euo pipefail
ROOT=/workspace/rfd

# --- AF2 environment, copied from the jobs that ran (174676, 174717) ---
# venv_af2, NOT the RFdiffusion venv. predict.py resolves its imports from this
# directory, so the cd is load-bearing and no PYTHONPATH is needed.
source "$ROOT/venv_af2/bin/activate"
cd "$ROOT/dl_binder_design/af2_initial_guess"
# ----------------------------------------------------------------------

mkdir -p "$ROOT/sc_out/$NAME"
nvidia-smi --query-gpu=name --format=csv,noheader

# -force_monomer: fold the binder on its own, which is the whole question here.
# No -debug: one bad input must not kill the rest of the batch.
python predict.py \
  -pdbdir "$ROOT/sc_in/$NAME" \
  -outpdbdir "$ROOT/sc_out/$NAME" \
  -scorefilename "$ROOT/$SCORE" \
  -maintain_res_numbering \
  -force_monomer

echo "=== done: $NAME ==="
# Judge success by output count and elapsed time, never by exit code: a job here
# has reported FAILED on success and COMPLETED having produced nothing.
echo "predicted structures: $(ls "$ROOT/sc_out/$NAME" | wc -l) of $NIN inputs"
echo "score lines: $(grep -c '^SCORE:' "$ROOT/$SCORE" || true)"
EOS

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL,NAME="$NAME",SCORE="$SCORE",NIN="$NIN" \
     bash "/workspace/rfd/_run_sc_${IDX}.sh"
