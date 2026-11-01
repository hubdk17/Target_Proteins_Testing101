#!/usr/bin/env bash
# Stage 2 -- ProteinMPNN: write sequences for the designed binder backbones.
#
# Runs INSIDE the Enroot container (paths are container paths). Launch with:
#
#   srun --partition=defq --cpus-per-task=4 --mem=16G --time=00:15:00 \
#        --container-image=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh \
#        --container-mounts=/home/$USER:/workspace \
#        --container-workdir=/workspace \
#        bash /workspace/rfd/run_mpnn.sh
#
# No GPU requested: ProteinMPNN is fast enough on CPU for a handful of
# backbones, and asking for a GPU only makes the job queue longer.
#
# Input:  ~/rfd/mpnn_in/*.pdb   (the selected RFdiffusion outputs)
# Output: ~/rfd/mpnn_out/seqs/*.fa

set -euo pipefail
ROOT=/workspace/rfd
source "$ROOT/venv/bin/activate"
cd "$ROOT/ProteinMPNN"

echo "=== 1/3 parsing the input structures ==="
python helper_scripts/parse_multiple_chains.py \
  --input_path="$ROOT/mpnn_in" \
  --output_path="$ROOT/mpnn_out/parsed.jsonl"

echo "=== 2/3 marking chain B as the one to design ==="
# Chain A is IL-6 and must keep its real sequence. Only chain B, the designed
# binder, gets new amino acids. Without this step ProteinMPNN would happily
# rewrite the target too.
python helper_scripts/assign_fixed_chains.py \
  --input_path="$ROOT/mpnn_out/parsed.jsonl" \
  --output_path="$ROOT/mpnn_out/assigned.jsonl" \
  --chain_list "B"

echo "=== 3/3 designing sequences ==="
# sampling_temp 0.1 is the usual low-temperature setting for binder design:
# it keeps sequences close to what the model thinks is best, rather than
# exploring. 8 sequences per backbone.
python protein_mpnn_run.py \
  --jsonl_path "$ROOT/mpnn_out/parsed.jsonl" \
  --chain_id_jsonl "$ROOT/mpnn_out/assigned.jsonl" \
  --out_folder "$ROOT/mpnn_out" \
  --num_seq_per_target 8 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1

echo "=== done ==="
ls -lh "$ROOT/mpnn_out/seqs/"
