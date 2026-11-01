#!/bin/bash
# BindCraft step 4: run it on IL-6 site III. Do NOT submit until
# 62_bindcraft_imports.sbatch has reported a clean sweep.
#
#   cd ~/rfd && sbatch scripts/63_bindcraft_run.sbatch              # 1 trajectory, a smoke test
#   cd ~/rfd && MAXTRAJ=20 NFINAL=5 sbatch scripts/63_bindcraft_run.sbatch
#
# Default is ONE trajectory on purpose: it proves the install end to end and
# measures the per-trajectory cost on this hardware before we commit an
# overnight array to a rate we guessed. The docs quote ~5 min for a
# 250-residue trajectory on an H100 and an A100 roughly 4x slower; IL-6 (157)
# plus a 40-55 binder is ~200-212 residues, so a single trajectory is the
# cheapest possible way to replace that quote with our own number.

#SBATCH --job-name=bc-il6
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=42G
#SBATCH --gres=gpu:1
#SBATCH --time=08:00:00
#SBATCH --exclude=dgxanode03
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh

# Overridable so the same tested script drives the pilot and the real batch.
MAXTRAJ="${MAXTRAJ:-1}"
NFINAL="${NFINAL:-1}"
TAG="${TAG:-p1}"
# Our own hotspots and binder band, so the comparison is against the same site
# our own arms were aimed at: IL-6 site III, Trp157 dominant.
HOTSPOTS="${HOTSPOTS:-A47,A155,A157,A158}"
LEN_MIN="${LEN_MIN:-40}"
LEN_MAX="${LEN_MAX:-55}"

mkdir -p "$HOME/rfd/logs"
test -f "$HOME/rfd/inputs/il6_clean.pdb" || { echo "il6_clean.pdb missing"; exit 1; }
test -d "$HOME/rfd/BindCraft"            || { echo "BindCraft missing - run 61"; exit 1; }

INNER="$HOME/rfd/_bc_run_${TAG}.sh"

cat > "$INNER" <<'EOS'
set -uo pipefail
ROOT=/workspace/rfd
BC="$ROOT/BindCraft"
PY="$ROOT/venv_bc/bin/python"
OUT="$ROOT/outputs/bc_il6_${TAG}"

echo "node: $(hostname)   tag: $TAG"
nvidia-smi --query-gpu=name --format=csv,noheader
date +'start: %F %T'

# Same guard as the import job: call python by path and refuse to proceed if
# that is not venv_bc.
PREFIX=$("$PY" -c 'import sys; print(sys.prefix)')
case "$PREFIX" in
  *venv_bc) echo "env: $PREFIX" ;;
  *) echo "*** WRONG ENV ($PREFIX) - stopping"; exit 1 ;;
esac

mkdir -p "$OUT"

# --- target settings -------------------------------------------------------
# Written here rather than committed, so the paths are always this cluster's
# and a stale copy cannot be picked up.
"$PY" - "$OUT" "$HOTSPOTS" "$LEN_MIN" "$LEN_MAX" "$NFINAL" <<'PYEOF'
import json, sys
out, hot, lmin, lmax, nfinal = sys.argv[1:6]
cfg = {
    "design_path": out + "/",
    "binder_name": "il6bc",
    "starting_pdb": "/workspace/rfd/inputs/il6_clean.pdb",
    "chains": "A",
    "target_hotspot_residues": hot,
    "lengths": [int(lmin), int(lmax)],
    "number_of_final_designs": int(nfinal),
}
with open("/workspace/rfd/inputs/bc_il6_target.json", "w") as fh:
    json.dump(cfg, fh, indent=4)
print("target settings:")
print(json.dumps(cfg, indent=2))
PYEOF

# --- advanced settings: theirs, with our overrides -------------------------
# Built by reading BindCraft's own default and changing named keys, so we
# inherit their defaults rather than retyping 60 of them, and every difference
# from stock is printed.
"$PY" - "$BC" "$MAXTRAJ" <<'PYEOF'
import json, sys
bc, maxtraj = sys.argv[1], int(sys.argv[2])
src = f"{bc}/settings_advanced/default_4stage_multimer.json"
cfg = json.load(open(src))
overrides = {
    # Bound the run. Stock is False, which means "keep going until you have
    # number_of_final_designs", i.e. open-ended.
    "max_trajectories": maxtraj,
    # Stock aborts the campaign if the acceptance rate falls under 1%. Over a
    # handful of trajectories that threshold is noise, and an early abort would
    # look like a failed install.
    "enable_rejection_check": False,
    # KEEP THE INTERMEDIATES. These three are True by default, and each one
    # deletes something our comparison needs: the binder predicted alone is the
    # fold-back measurement, and the unrelaxed complex is what our own scoring
    # path takes as input. Letting the tool clean up after itself would leave us
    # with only its own verdict and nothing to re-score.
    "remove_binder_monomer": False,
    "remove_unrelaxed_complex": False,
    "remove_unrelaxed_trajectory": False,
    # Plots and animations are per-trajectory disk and time we have no use for.
    "save_design_animations": False,
    "save_design_trajectory_plots": False,
    "zip_animations": False,
    "zip_plots": False,
}
print("overrides against stock default_4stage_multimer.json:")
for k, v in overrides.items():
    print(f"  {k}: {cfg.get(k, '<absent>')} -> {v}")
cfg.update(overrides)
dst = "/workspace/rfd/inputs/bc_il6_advanced.json"
json.dump(cfg, open(dst, "w"), indent=4)
print(f"wrote {dst}")
PYEOF

# --- run -------------------------------------------------------------------
cd "$BC"
echo "--- bindcraft.py ---"
"$PY" -u bindcraft.py \
  --settings /workspace/rfd/inputs/bc_il6_target.json \
  --filters  "$BC/settings_filters/default_filters.json" \
  --advanced /workspace/rfd/inputs/bc_il6_advanced.json

echo "--- what arrived ---"
date +'end: %F %T'
# A Slurm state has lied in both directions on this cluster, so read success off
# the output files and the elapsed time instead.
for f in trajectory_stats.csv mpnn_design_stats.csv final_design_stats.csv failure_csv.csv; do
  if [ -f "$OUT/$f" ]; then
    printf '  %-26s %s data rows\n' "$f" "$(( $(wc -l < "$OUT/$f") - 1 ))"
  else
    printf '  %-26s ABSENT\n' "$f"
  fi
done
echo "  pdb files: $(find "$OUT" -name '*.pdb' 2>/dev/null | wc -l)"
find "$OUT" -maxdepth 1 -type d | sed 's|^|  dir: |'
exit 0
EOS

export TAG MAXTRAJ NFINAL HOTSPOTS LEN_MIN LEN_MAX
# --export=ALL and never a NAME=VALUE list: that list is comma-separated, and it
# truncated a hotspot list at the first comma on this project once already.
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL \
     bash "/workspace/rfd/_bc_run_${TAG}.sh"
