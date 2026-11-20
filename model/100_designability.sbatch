#!/bin/bash
# Step 4 run: designability / self-consistency of our generated backbones.
#   bb_samples (from 98) --ProteinMPNN--> sequences --thread--> AF2 monomer fold
#   --> scRMSD(prediction, generated backbone).  Designable if scRMSD < 2 A.
#
# This is the SAME self-consistency AF2 (initial-guess, -force_monomer) the lab
# uses for binders, so it is a calibrated, discriminating test (it failed ~15% of
# real binder designs, so it is not a check-that-cannot-fail). Judge by the
# designability block at the end, never by Slurm State.
#
# First generate the backbones (quick, separate):
#   srun --gres=gpu:1 --container-image=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh \
#     --container-mounts=/home/$USER:/workspace --container-workdir=/workspace \
#     python /workspace/rfd/scripts/model/98_sample_backbone.py \
#       --ckpt /workspace/rfd/corpus_v1/bbflow.pt --n 20 --length 80 --steps 100 \
#       --out-dir /workspace/rfd/corpus_v1/bb_samples
# Then:
#   cd ~/rfd/scripts && sbatch model/100_designability.sbatch
#   # NUM_SEQ=8 by default; a bigger choice costs AF2 time: NUM_SEQ=16 sbatch ...

#SBATCH --job-name=bbdesign
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --gres=gpu:1
#SBATCH --exclude=dgxanode03
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
NUM_SEQ="${NUM_SEQ:-8}"
STAMP=$(date +%m%d_%H%M)

test -n "$(ls $HOME/rfd/corpus_v1/bb_samples/*.pdb 2>/dev/null)" || {
  echo "no backbones in ~/rfd/corpus_v1/bb_samples - run 98_sample_backbone.py first"; exit 1; }
test -f "$HOME/rfd/venv/bin/activate" || { echo "ProteinMPNN venv missing"; exit 1; }
test -f "$HOME/rfd/venv_af2/bin/activate" || { echo "venv_af2 missing"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID  NUM_SEQ=$NUM_SEQ"
INNER="$HOME/rfd/_run_bbdesign.sh"

cat > "$INNER" <<'EOS'
set -euo pipefail
ROOT=/workspace/rfd
BB=$ROOT/corpus_v1/bb_samples
MP=$ROOT/corpus_v1/mpnn_bb
TH=$ROOT/corpus_v1/bb_threaded
AF=$ROOT/corpus_v1/bb_af2
SC=$ROOT/corpus_v1/bb_designability_${STAMP}.sc
rm -rf "$MP" "$TH" "$AF"; mkdir -p "$MP" "$TH" "$AF"
nvidia-smi --query-gpu=name --format=csv,noheader

echo "=== 1/4 ProteinMPNN: design the whole chain ==="
source "$ROOT/venv/bin/activate"
cd "$ROOT/ProteinMPNN"
python helper_scripts/parse_multiple_chains.py --input_path="$BB" --output_path="$MP/parsed.jsonl"
# No assign_fixed_chains: nothing fixed -> MPNN designs every residue (monomer).
python protein_mpnn_run.py --jsonl_path "$MP/parsed.jsonl" --out_folder "$MP" \
  --num_seq_per_target "$NUM_SEQ" --sampling_temp 0.1 --seed 37 --batch_size 1
deactivate
echo "fasta files: $(ls "$MP/seqs" | wc -l)"

echo "=== 2/4 thread sequences onto backbones ==="
python "$ROOT/scripts/model/thread_monomer.py" --backbones "$BB" --seqs "$MP/seqs" --out "$TH"

echo "=== 3/4 AF2 monomer fold (self-consistency) ==="
source "$ROOT/venv_af2/bin/activate"
cd "$ROOT/dl_binder_design/af2_initial_guess"
python predict.py -pdbdir "$TH" -outpdbdir "$AF" -scorefilename "$SC" \
  -maintain_res_numbering -force_monomer

echo "=== 4/4 designability (scRMSD < 2 A) ==="
python "$ROOT/scripts/model/99_designability.py" --ref-dir "$BB" --preds-flat "$AF"
echo "threaded: $(ls "$TH"/*.pdb 2>/dev/null | wc -l)   af2 preds: $(ls "$AF"/*.pdb 2>/dev/null | wc -l)"
echo "=== done ==="
EOS

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL,NUM_SEQ="$NUM_SEQ",STAMP="$STAMP" \
     bash "/workspace/rfd/_run_bbdesign.sh"
