#!/bin/bash
# STAGE 3 for the arm groups: AF2 initial guess on the assembled COMPLEX, one
# array task per arm.
#
# Method: /mnt/project-files/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md
# Sibling of 40_gate2_selfconsistency.sbatch, which folds the binder ALONE with
# -force_monomer. This one scores the complex and so does NOT pass that flag;
# that one flag is the whole difference in what the two jobs measure.
#
#   cd ~/rfd && TABLE=$HOME/rfd/inputs/arms_m1.tsv sbatch --array=1-6%2 \
#                 --nodelist=dgxanode04 --exclude=dgxanode03 62_af2_arms.sbatch
#
# --nodelist is not optional for a comparison. AF2 here is not deterministic even
# on one node and the spread is heavy-tailed, so letting the arms land on
# different hardware adds a between-node term to the very difference being read.
# Pin all arms of a group to one node, and pin the same node for a re-run.
#
# --exclude=dgxanode03 as well: its container disk fails intermittently while
# unpacking the squashfs (FAILED inside 2 min, zero-byte log). One passing short
# test there does not clear it.
#
# Environment lines are copied from the jobs that ran (174676, 174717). venv_af2
# and its two one-time patches already exist on this account - do NOT rebuild
# them. A `Bio.Data.SCOPData` import error or a `take_current_trace`
# AttributeError means the tree was rebuilt and lost them:
#   - mmcif_parsing.py must use PDBData.protein_letters_3to1_extended
#   - venv_af2 must hold dm-haiku==0.0.13 (0.0.14+ breaks jax 0.4.20)

#SBATCH --job-name=af2-arms
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --exclude=dgxanode03
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
TABLE="${TABLE:-$HOME/rfd/inputs/arms_m1.tsv}"
TAG="${TAG:-m1}"

test -s "$TABLE" || { echo "$TABLE missing - run 60_equal_n_arms.py"; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:-1}"
NAME=$(sed -n "${IDX}p" "$TABLE" | cut -f1)
test -n "$NAME" || { echo "no line $IDX in $TABLE"; exit 1; }

INDIR="$HOME/rfd/cplx_in/$NAME"
NIN=$(ls "$INDIR"/*.pdb 2>/dev/null | wc -l)
test "$NIN" -gt 0 || { echo "no complexes in $INDIR - run 43_build_complex.py $NAME"; exit 1; }

# Fail before the GPU is allocated and the container started, not after.
test -f "$HOME/rfd/venv_af2/bin/activate" || { echo "venv_af2 missing - do not rebuild it blindly, check why"; exit 1; }
test -f "$HOME/rfd/dl_binder_design/af2_initial_guess/predict.py" || { echo "predict.py missing"; exit 1; }

# Equal n across the arms of a group is the whole point of 60_equal_n_arms.py and
# it is cheap to confirm here, where the number that will be compared is in hand.
GROUP="${NAME%_*}"
for A in wt mut ctl; do
  O=$(ls "$HOME/rfd/cplx_in/${GROUP}_${A}"/*.pdb 2>/dev/null | wc -l)
  echo "  arm ${GROUP}_${A}: $O complex input(s)"
done

# Fresh score filename every run. The score file is APPEND-ONLY: a reused name
# interleaves two runs with no warning and no way to tell them apart afterwards.
STAMP=$(date +%m%d_%H%M%S)
SCORE="af2_${NAME}_${TAG}_${STAMP}.sc"

echo "node: $(hostname)  job: $SLURM_JOB_ID  task: $IDX"
echo "arm: $NAME   inputs: $NIN   score file: $SCORE"

INNER="$HOME/rfd/_run_af2_arm_${IDX}.sh"

cat > "$INNER" <<'EOS'
set -euo pipefail
ROOT=/workspace/rfd

# --- AF2 environment, exactly as in the jobs that ran (174676, 174717) ---
# venv_af2, NOT the RFdiffusion venv. predict.py resolves its imports from this
# directory, so the cd is load-bearing and no PYTHONPATH is needed; a parent on
# PYTHONPATH risks shadowing.
source "$ROOT/venv_af2/bin/activate"
cd "$ROOT/dl_binder_design/af2_initial_guess"
# ------------------------------------------------------------------------

mkdir -p "$ROOT/af2_out/$NAME"
nvidia-smi --query-gpu=name --format=csv,noheader

# No -force_monomer: the complex is what is being scored here.
# No -debug: one bad input must not kill the rest of the arm.
python predict.py \
  -pdbdir "$ROOT/cplx_in/$NAME" \
  -outpdbdir "$ROOT/af2_out/$NAME" \
  -scorefilename "$ROOT/$SCORE" \
  -maintain_res_numbering

echo "=== done: $NAME ==="
# Success is output files and elapsed time, never the exit code: a job here has
# reported FAILED on success and COMPLETED having produced nothing.
echo "predicted structures: $(ls "$ROOT/af2_out/$NAME" | wc -l) of $NIN inputs"
# One of these SCORE: lines is the header, so this count is one more than the
# number of scored designs.
echo "score lines incl. header: $(grep -c '^SCORE:' "$ROOT/$SCORE" || true)"
EOS

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL,NAME="$NAME",SCORE="$SCORE",NIN="$NIN" \
     bash "/workspace/rfd/_run_af2_arm_${IDX}.sh"
