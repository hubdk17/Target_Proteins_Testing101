#!/bin/bash
# Featurize corpus_v2 in parallel shards. 81_featurize_corpus.py is pure Python
# (stdlib only -- no numpy), so this runs on the bare compute node, NO container
# and NO GPU. Each array task featurizes 1/N of corpus_v2/raw with the enriched
# columns and the self-verifying holdout leak-check.
#
# Submit from the scripts repo root AFTER the corpus_v2 fetch has finished:
#   cd ~/rfd/scripts && sbatch --array=1-20 model/88_featurize_v2.sbatch
# (N defaults to 20; override with SHARDS=40 and --array=1-40 to match.)
#
# After the array finishes, concat the shards (89 does this for you, or by hand):
#   cd ~/rfd/corpus_v2 && head -1 if_1.tsv > interface_features.tsv \
#     && for f in if_*.tsv; do tail -n +2 "$f"; done >> interface_features.tsv
#
# Judge success by line counts in the shard .tsv files, not the Slurm State.

#SBATCH --job-name=feat-v2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
N="${SHARDS:-20}"
I="${SLURM_ARRAY_TASK_ID:-1}"
RAW="$HOME/rfd/corpus_v2/raw"
OUT="$HOME/rfd/corpus_v2/if_${I}.tsv"
SCRIPT="$HOME/rfd/scripts/model/81_featurize_corpus.py"

test -d "$RAW" || { echo "$RAW missing - fetch corpus_v2 first"; exit 1; }

echo "node: $(hostname)  shard $I of $N"
# 81 refuses to overwrite its --out (no silent append), so clear a stale shard
# from a previous attempt before re-running this task.
rm -f "$OUT"
python3 "$SCRIPT" "$RAW" --out "$OUT" --shard "${I}/${N}"
echo "=== shard $I/$N done: $(wc -l < "$OUT") rows ==="
