#!/bin/bash
# STAGE 2 across targets: ProteinMPNN sequences for the backbones that passed
# gate 1. One array task per target. This is the missing link between
# 32_gate1_compactness.sh and thread_multi.py.
#
# Method: /mnt/project-files/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md
#
#   cd ~/rfd && sbatch --array=1 33_mpnn_multi.sbatch        # pd1 only
#   cd ~/rfd && sbatch --array=1-6 33_mpnn_multi.sbatch      # batch 1
#   cd ~/rfd && NUM_SEQ=48 sbatch --array=1 33_mpnn_multi.sbatch
#
# Only backbones listed in gate1/<name>.pass are used. Writing sequences for a
# backbone that already failed compactness is work the gate order exists to save.
#
# A GPU is requested deliberately. The single-target run_mpnn.sh runs on CPU
# because a handful of backbones is fine there, but CPU is ~11 min per structure
# against seconds on a GPU, and a batch is tens of structures per target.

#SBATCH --job-name=mpnn-b1
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=00:40:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
TABLE="$HOME/rfd/inputs/targets.tsv"

# 8 exposed the problem on IL-6 but is on the low side; gate 2 picks among these,
# so give it a real choice.
NUM_SEQ="${NUM_SEQ:-16}"
SAMPLING_TEMP="${SAMPLING_TEMP:-0.1}"
SEED="${SEED:-37}"

test -s "$TABLE" || { echo "$TABLE missing - run 02_prep_targets.py"; exit 1; }

IDX="${SLURM_ARRAY_TASK_ID:-1}"
NAME=$(sed -n "${IDX}p" "$TABLE" | cut -f1)
test -n "$NAME" || { echo "no line $IDX in $TABLE"; exit 1; }

PASS="$HOME/rfd/gate1/${NAME}.pass"
test -s "$PASS" || { echo "no $PASS - run 32_gate1_compactness.sh first"; exit 1; }

# Stage the survivors into a per-target input directory. Copies, not symlinks:
# the container sees a different root and a symlink to a host path would dangle.
INDIR="$HOME/rfd/mpnn_in/$NAME"
rm -rf "$INDIR"; mkdir -p "$INDIR"
N=0
while read -r f; do
  [ -f "$f" ] && { cp "$f" "$INDIR/"; N=$((N+1)); }
done < "$PASS"
test "$N" -gt 0 || { echo "none of the paths in $PASS exist"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID  task: $IDX"
echo "target: $NAME   backbones: $N   sequences each: $NUM_SEQ   total: $((N*NUM_SEQ))"

INNER="$HOME/rfd/_run_mpnn_${IDX}.sh"

cat > "$INNER" <<'EOS'
set -euo pipefail
ROOT=/workspace/rfd
source "$ROOT/venv/bin/activate"
cd "$ROOT/ProteinMPNN"
OUT="$ROOT/mpnn_out_$NAME"
mkdir -p "$OUT"

echo "=== 1/3 parsing structures ==="
python helper_scripts/parse_multiple_chains.py \
  --input_path="$ROOT/mpnn_in/$NAME" \
  --output_path="$OUT/parsed.jsonl"

echo "=== 2/3 marking chain B as the one to design ==="
# Chain A is the target and must keep its real sequence. Only chain B, the
# binder, gets new amino acids.
python helper_scripts/assign_fixed_chains.py \
  --input_path="$OUT/parsed.jsonl" \
  --output_path="$OUT/assigned.jsonl" \
  --chain_list "B"

echo "=== 3/3 designing sequences ==="
python protein_mpnn_run.py \
  --jsonl_path "$OUT/parsed.jsonl" \
  --chain_id_jsonl "$OUT/assigned.jsonl" \
  --out_folder "$OUT" \
  --num_seq_per_target "$NUM_SEQ" \
  --sampling_temp "$SAMPLING_TEMP" \
  --seed "$SEED" \
  --batch_size 1

echo "=== done: $NAME ==="
# Success is output files, never the exit code. One .fa per backbone, and each
# should hold NUM_SEQ sequences.
echo "fasta files: $(ls "$OUT/seqs" | wc -l) of $N backbones"
echo "sequences in the first file: $(grep -c '^>T=' "$(ls "$OUT"/seqs/*.fa | head -1)" || true)"
EOS

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL,NAME="$NAME",NUM_SEQ="$NUM_SEQ",SAMPLING_TEMP="$SAMPLING_TEMP",SEED="$SEED",N="$N" \
     bash "/workspace/rfd/_run_mpnn_${IDX}.sh"
