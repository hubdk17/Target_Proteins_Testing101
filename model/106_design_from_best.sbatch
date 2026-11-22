#!/bin/bash
# FULL design generation in ONE job: sample backbones from our best model, then run
# the calibrated self-consistency pipeline end to end --
#   98 sample (alpha3 + mild velocity-clamp) --ProteinMPNN--> sequences --thread-->
#   AF2 monomer fold --> scRMSD(prediction, backbone).  Designable if scRMSD < 2 A.
# Produces actual DESIGNS (sequences + AF2-predicted structures) plus a real
# designability rate. Judge by the designability block at the end, never Slurm State.
#
# HONEST NOTE: our backbones are not yet clean (sampler analysis shows a
# connectivity<->clash tension: the best ops point still has ~19 clashes/chain), so
# expect a LOW designable rate on this first pass. This gives us designs to look at
# and a baseline number, not a finished generator.
#
#   cd ~/rfd/scripts && git checkout claude/project-thread-ho43or && git pull \
#     && sbatch model/106_design_from_best.sbatch
#   # knobs: CKPT=bbflow_rotsmooth_s05 N=60 VCLAMP=0.5 NUM_SEQ=16 sbatch model/106_...

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
CKPT="${CKPT:-bbflow_alpha3}"      # best-balanced model
N="${N:-40}"                       # backbones to generate
LEN="${LEN:-80}"
VCLAMP="${VCLAMP:-1.0}"            # mild clamp: tightens connectivity w/o full collapse
NUM_SEQ="${NUM_SEQ:-8}"
STAMP=$(date +%m%d_%H%M)

test -f "$HOME/rfd/corpus_v1/$CKPT.pt" || { echo "checkpoint $CKPT.pt missing"; exit 1; }
test -f "$HOME/rfd/venv/bin/activate" || { echo "ProteinMPNN venv missing"; exit 1; }
test -f "$HOME/rfd/venv_af2/bin/activate" || { echo "venv_af2 missing"; exit 1; }

echo "node: $(hostname)  job: $SLURM_JOB_ID  CKPT=$CKPT N=$N VCLAMP=$VCLAMP NUM_SEQ=$NUM_SEQ"
INNER="$HOME/rfd/_run_bbdesign_best.sh"

cat > "$INNER" <<'EOS'
set -euo pipefail
ROOT=/workspace/rfd
BB=$ROOT/corpus_v1/design_bb_${STAMP}
MP=$ROOT/corpus_v1/design_mpnn_${STAMP}
TH=$ROOT/corpus_v1/design_threaded_${STAMP}
AF=$ROOT/corpus_v1/design_af2_${STAMP}
SC=$ROOT/corpus_v1/design_${STAMP}.sc
mkdir -p "$BB" "$MP" "$TH" "$AF"
nvidia-smi --query-gpu=name --format=csv,noheader

echo "=== 0/4 sample $N backbones from $CKPT (vclamp $VCLAMP) ==="
python "$ROOT/scripts/model/98_sample_backbone.py" \
  --ckpt "$ROOT/corpus_v1/$CKPT.pt" --n "$N" --length "$LEN" --steps 200 \
  --vclamp "$VCLAMP" --out-dir "$BB"
echo "backbones: $(ls "$BB"/*.pdb 2>/dev/null | wc -l)"

echo "=== 1/4 ProteinMPNN: design the whole chain ==="
source "$ROOT/venv/bin/activate"
cd "$ROOT/ProteinMPNN"
python helper_scripts/parse_multiple_chains.py --input_path="$BB" --output_path="$MP/parsed.jsonl"
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
echo "backbones: $(ls "$BB"/*.pdb|wc -l)  threaded: $(ls "$TH"/*.pdb|wc -l)  af2: $(ls "$AF"/*.pdb|wc -l)"
echo "designs (AF2 structures) in: $AF    scores: $SC"
echo "=== done ==="
EOS

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL,CKPT="$CKPT",N="$N",LEN="$LEN",VCLAMP="$VCLAMP",NUM_SEQ="$NUM_SEQ",STAMP="$STAMP" \
     bash "/workspace/rfd/_run_bbdesign_best.sh"
