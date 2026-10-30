#!/bin/bash
# Batch 1 RFdiffusion backbones: one array task per target.
# Generalised from 30_rfd_il6.sbatch. Everything target-specific is read from
# inputs/targets.tsv, which 02_prep_targets.py writes after measuring each PDB.
#
# Submit from ~/rfd, after running 02_prep_targets.py and READING its output:
#
#   cd ~/rfd && sbatch --array=1-6 31_rfd_multi.sbatch
#
# One target only (line 1 of the table is pd1, the diagnostic control):
#   cd ~/rfd && sbatch --array=1 31_rfd_multi.sbatch
#
# Smaller/cheaper first, which is the house preference:
#   NUM_DESIGNS=5 sbatch --array=1 31_rfd_multi.sbatch

#SBATCH --job-name=rfd-b1
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=logs/rfd_b1_%A_%a.out
#SBATCH --error=logs/rfd_b1_%A_%a.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
TABLE="$HOME/rfd/inputs/targets.tsv"

NUM_DESIGNS="${NUM_DESIGNS:-50}"
DIFFUSER_T="${DIFFUSER_T:-50}"
# Range, not a fixed length: a model free to choose its own length does not have
# to stretch a fixed residue count into a rod. 60% pass vs 31% on IL-6.
BINDER_LEN="${BINDER_LEN:-40-55}"
# Run tag keeps output names distinct; a duplicate name has overwritten a batch
# here before.
TAG="${TAG:-b1}"

test -f "$HOME/rfd/venv/.complete" || { echo "venv missing - run 20_buildenv.sbatch"; exit 1; }
test -s "$TABLE"                   || { echo "$TABLE missing - run 02_prep_targets.py"; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:-1}"
ROW=$(sed -n "${IDX}p" "$TABLE")
test -n "$ROW" || { echo "no line $IDX in $TABLE"; exit 1; }

NAME=$(echo "$ROW"   | cut -f1)
PDBID=$(echo "$ROW"  | cut -f2)
CONTIG=$(echo "$ROW" | cut -f5)
HOTSPOTS=$(echo "$ROW" | cut -f6)
INPDB="$HOME/rfd/inputs/${NAME}_clean.pdb"

test -f "$INPDB"      || { echo "$INPDB missing - rerun 02_prep_targets.py"; exit 1; }
test -n "$CONTIG"     || { echo "empty contig for $NAME"; exit 1; }
test -n "$HOTSPOTS"   || { echo "empty hotspots for $NAME"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID  task: $IDX"
echo "target: $NAME ($PDBID)"
echo "contig: [${CONTIG}/0 ${BINDER_LEN}]"
NHOT=$(printf '%s' "$HOTSPOTS" | tr ',' '\n' | grep -c .)
echo "hotspots: ${HOTSPOTS}  (${NHOT})"
echo "designs: $NUM_DESIGNS   T: $DIFFUSER_T"

# Per-task inner script. A shared filename would have array tasks overwrite
# each other's script while it is being read.
INNER="$HOME/rfd/_run_multi_${IDX}.sh"

cat > "$INNER" <<'EOS'
set -euo pipefail
ROOT=/workspace/rfd
source "$ROOT/venv/bin/activate"
export DGLBACKEND=pytorch
cd "$ROOT/RFdiffusion"
mkdir -p "$ROOT/outputs"
nvidia-smi --query-gpu=name --format=csv,noheader

# The radius-of-gyration potential is a switch so it can be A/B tested. Whether
# it fires at all is an open question as of 2026-10-30: the IL-6 arms came back
# at 23-24% pass with it supposedly on, no better than a run without it.
# POTENTIALS=0 turns it off and changes nothing else.
POT_ARGS=()
if [ "${POTENTIALS:-1}" = "1" ]; then
  POT_ARGS=(
    'potentials.guiding_potentials=["type:binder_ROG,weight:3"]'
    potentials.guide_scale=2
    potentials.guide_decay=quadratic
  )
fi

# Print the exact arguments, so a later log says what was actually passed
# instead of leaving us to infer it. Absence of a log line is not evidence.
echo "POTENTIALS=${POTENTIALS:-1}"
echo "potential args: ${POT_ARGS[*]-(none)}"

# Compare what arrived with what was sent. The hotspot list crossed a process
# boundary that eats commas once already; this is the check that would have
# caught it on design 1 instead of after 300.
GOT=$(printf '%s' "$HOTSPOTS" | tr ',' '\n' | grep -c .)
echo "hotspots received: $HOTSPOTS ($GOT of ${NHOT:-?} expected)"
if [ -n "${NHOT:-}" ] && [ "$GOT" != "$NHOT" ]; then
  echo "*** HOTSPOTS TRUNCATED in transit - refusing to waste a GPU on a mis-aimed run"
  exit 1
fi

python scripts/run_inference.py \
  inference.output_prefix="$ROOT/outputs/${NAME}_${TAG}" \
  inference.model_directory_path="$ROOT/models" \
  inference.input_pdb="/workspace/rfd/inputs/${NAME}_clean.pdb" \
  inference.ckpt_override_path="$ROOT/models/Complex_base_ckpt.pt" \
  inference.num_designs="$NUM_DESIGNS" \
  diffuser.T="$DIFFUSER_T" \
  "contigmap.contigs=[${CONTIG}/0 ${BINDER_LEN}]" \
  "ppi.hotspot_res=[${HOTSPOTS}]" \
  ${POT_ARGS[@]+"${POT_ARGS[@]}"}

echo "=== done: ${NAME} ==="
ls -lh "$ROOT/outputs" | grep "${NAME}_${TAG}" | tail -3
EOS

# Export into the job environment and propagate the whole thing with
# --export=ALL. Do NOT list NAME=VALUE pairs in --export: that list is itself
# comma-separated, so Slurm truncated HOTSPOTS="A107,A104,A108,A118,A146" at the
# first comma and RFdiffusion received a single hotspot residue. That is what
# happened to all 300 backbones on 2026-10-29, silently, because nothing
# compared what was sent with what arrived.
export POTENTIALS="${POTENTIALS:-1}"
export NUM_DESIGNS DIFFUSER_T BINDER_LEN CONTIG HOTSPOTS NAME TAG NHOT

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL \
     bash "/workspace/rfd/_run_multi_${IDX}.sh"
