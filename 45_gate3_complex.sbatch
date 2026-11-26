#!/bin/bash
# GATE 3: score the designed COMPLEX with AF2 initial guess and read
# pae_interaction. This is the real discriminator - gate 2 only asks whether a
# binder folds back to its backbone; gate 3 asks whether the binder+target
# complex is predicted to be a confident interface.
#
# Method: /mnt/project-files/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md
# Copied from 40_gate2_selfconsistency.sbatch; the ONLY scientific differences
# are the input dir (cplx_in, built by 43_build_complex.py) and that there is NO
# -force_monomer, so AF2 folds the two-chain complex and reports pae_interaction.
#
#   cd ~/rfd && bash 41_gate2_report.sh mdm2            # -> gate2/mdm2.pass
#   cd ~/rfd && python3 43_build_complex.py mdm2        # -> cplx_in/mdm2/*.pdb
#   cd ~/rfd && TABLE=$HOME/rfd/scripts/mdm2_blind/mdm2_rfd.tsv TAG=mdm2b1 \
#     sbatch --exclude=dgxanode03 --array=1 45_gate3_complex.sbatch
#
# GATE: keep designs with pae_interaction below ~10. The metric saturates near
# ~31.75 (off-site designs pin there), so it is a good gate but a poor ranker,
# and noisiest in the 9-20 band that matters. A design within ~1 of the gate
# needs replicate scoring before it is called; pin --nodelist= when comparing.
# Never rank by pLDDT.

#SBATCH --job-name=g3-b1
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
# Overridable so a single-target table (e.g. mdm2_blind/mdm2_rfd.tsv) can drive it.
TABLE="${TABLE:-$HOME/rfd/inputs/targets.tsv}"
TAG="${TAG:-b1}"

test -s "$TABLE" || { echo "$TABLE missing - run 02_prep_targets.py"; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:-1}"
NAME=$(sed -n "${IDX}p" "$TABLE" | cut -f1)
test -n "$NAME" || { echo "no line $IDX in $TABLE"; exit 1; }

INDIR="$HOME/rfd/cplx_in/$NAME"
NIN=$(ls "$INDIR"/*.pdb 2>/dev/null | wc -l)
test "$NIN" -gt 0 || { echo "no complexes in $INDIR - run 43_build_complex.py"; exit 1; }

# Fail here rather than after the GPU has been allocated and the container started.
test -f "$HOME/rfd/venv_af2/bin/activate" || { echo "venv_af2 missing - do not rebuild it blindly, check why"; exit 1; }
test -f "$HOME/rfd/dl_binder_design/af2_initial_guess/predict.py" || { echo "predict.py missing"; exit 1; }

# Fresh score filename every run: the score file is append-only, and a reused
# name mixes two runs together with no warning.
STAMP=$(date +%m%d_%H%M)
SCORE="cp_${NAME}_${TAG}_${STAMP}.sc"

echo "node: $(hostname)  job: $SLURM_JOB_ID  task: $IDX"
echo "target: $NAME   complexes: $NIN   score file: $SCORE"
echo "expect ~25 s per structure, so roughly $(( NIN * 25 / 60 )) min"

INNER="$HOME/rfd/_run_g3_${IDX}.sh"

cat > "$INNER" <<'EOS'
set -euo pipefail
ROOT=/workspace/rfd

# --- AF2 environment, same as gate 2 (jobs 174676, 174717) ---
source "$ROOT/venv_af2/bin/activate"
cd "$ROOT/dl_binder_design/af2_initial_guess"
# ------------------------------------------------------------

mkdir -p "$ROOT/cplx_out/$NAME"
nvidia-smi --query-gpu=name --format=csv,noheader

# No -force_monomer: fold the two-chain complex so pae_interaction is reported.
# No -debug: one bad input must not kill the rest of the batch.
python predict.py \
  -pdbdir "$ROOT/cplx_in/$NAME" \
  -outpdbdir "$ROOT/cplx_out/$NAME" \
  -scorefilename "$ROOT/$SCORE" \
  -maintain_res_numbering

echo "=== done: $NAME ==="
# Judge success by output count and elapsed time, never by exit code.
echo "predicted structures: $(ls "$ROOT/cplx_out/$NAME" | wc -l) of $NIN inputs"
echo "score lines: $(grep -c '^SCORE:' "$ROOT/$SCORE" || true)"
EOS

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL,NAME="$NAME",SCORE="$SCORE",NIN="$NIN" \
     bash "/workspace/rfd/_run_g3_${IDX}.sh"
