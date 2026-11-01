#!/bin/bash
# BindCraft step 5: install the five packages the import sweep named, then
# re-run the sweep. Short job, one GPU.
#
#   cd ~/rfd && sbatch scripts/64_bindcraft_deps.sbatch
#
# What 62 established (job 176356, dgxanode04): colabdesign 1.1.3 built and
# installed, jax 0.4.20 sees cuda(id=0) from venv_bc, nothing is broken by a
# version conflict, and the gap is optax, chex, immutabledict, py3Dmol, seaborn.
#
# Two corrections to 62 are folded in here.
#
# 1. The sweep MISATTRIBUTED a failure. It listed colabdesign as missing when
#    pip had just installed it successfully: importing colabdesign raises
#    ModuleNotFoundError for *optax*, and the sweep blamed the module it was
#    importing rather than the name that was actually absent. The entry-point
#    check gave the true reason. Fixed below by comparing the error's .name
#    against the module asked for, so a dependency's absence is reported as
#    such.
#
# 2. VERSIONS ARE PINNED, not latest. Current optax requires jax >= 0.4.27 and
#    current chex is built against the same era. Ours is deliberately 0.4.20,
#    because the only GPU wheel we have is jaxlib 0.4.20+cuda12.cudnn89 and the
#    cluster has cuDNN 8.9. So pick the contemporaries of jax 0.4.20:
#    optax 0.1.7 and chex 0.1.85, both from late 2023, both declaring a jax
#    floor at or below 0.4.20. Installing latest with --no-deps would succeed
#    and then fail at import on an attribute 0.4.20 does not have, which is
#    exactly how dm-haiku 0.0.14 cost us a job already.

#SBATCH --job-name=bc-deps
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --exclude=dgxanode03
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
mkdir -p "$HOME/rfd/logs"
INNER="$HOME/rfd/_bc_deps.sh"

cat > "$INNER" <<'EOS'
set -uo pipefail
ROOT=/workspace/rfd
PY="$ROOT/venv_bc/bin/python"

echo "node: $(hostname)"
PREFIX=$("$PY" -c 'import sys; print(sys.prefix)')
case "$PREFIX" in
  *venv_bc) echo "env: $PREFIX" ;;
  *) echo "*** WRONG ENV ($PREFIX) - stopping"; exit 1 ;;
esac

# The NGC image sets PIP_CONSTRAINT to a file pinning everything it ships, and
# it makes unrelated installs fail with ResolutionImpossible. Replace it with a
# file holding only the container's torch pin: nothing is left to collide, and
# pip still will not try to replace the pre-release torch with a 2.5 GB
# download.
"$PY" -c "import torch;print('torch=='+torch.__version__)" > /workspace/rfd/pipc_bc.txt
export PIP_CONSTRAINT=/workspace/rfd/pipc_bc.txt
echo "constraint file: $(cat /workspace/rfd/pipc_bc.txt)"

# absl-py and toolz are chex's own dependencies; named explicitly because
# --no-deps will not pull them, and absl-py is probably already here via
# tensorflow-cpu. Naming a satisfied package costs nothing.
PKGS="optax==0.1.7 chex==0.1.85 absl-py toolz immutabledict py3Dmol seaborn"
echo "--- dry run first ---"
# A dry run catches a version number pip would pick that we do not want. It
# cannot catch an API mismatch that only appears at import, which is why the
# sweep below is the real check and not a formality.
"$PY" -m pip install --no-deps --dry-run $PKGS 2>&1 | tail -12

echo "--- installing ---"
"$PY" -m pip install --no-deps $PKGS 2>&1 | tail -8

echo "--- import sweep (corrected attribution) ---"
"$PY" - <<'PYEOF'
import importlib
mods = [
    "jax", "jaxlib", "haiku", "pyrosetta",
    "colabdesign", "optax", "chex", "ml_collections", "tree", "immutabledict",
    "joblib", "py3Dmol", "absl", "toolz",
    "numpy", "scipy", "pandas", "matplotlib", "seaborn", "Bio", "tqdm",
]
ok, missing, viadep, broken = [], [], [], []
for m in mods:
    try:
        mod = importlib.import_module(m)
        ok.append((m, getattr(mod, "__version__", "?")))
    except ModuleNotFoundError as e:
        # The fix: blame the name that is actually absent, not the one we asked
        # for. Anything else reads as "the package we just installed is
        # missing", which is what 62 reported.
        if getattr(e, "name", None) in (None, m) or m.startswith(str(e.name)):
            missing.append(m)
        else:
            viadep.append((m, e.name))
    except Exception as e:
        broken.append((m, f"{type(e).__name__}: {e}"))

print("\nPRESENT:")
for m, v in ok:
    print(f"  {m:16s} {v}")
print("\nMISSING (the package itself):")
print("  none" if not missing else "\n".join(f"  {m}" for m in missing))
print("\nPRESENT BUT BLOCKED (installed; a dependency is absent):")
print("  none" if not viadep else "\n".join(f"  {m:16s} needs {d}" for m, d in viadep))
print("\nBROKEN (imports but errors - version conflict):")
print("  none" if not broken else "\n".join(f"  {m}\n      {e}" for m, e in broken))

print("\n--- the entry point BindCraft calls ---")
try:
    from colabdesign import mk_afdesign_model
    print("  mk_afdesign_model OK")
except Exception as e:
    print(f"  FAILED  {type(e).__name__}: {e}")

print("\n--- and the one BindCraft itself imports ---")
# If this passes, the environment is done and the next job is the real run.
try:
    from colabdesign import clear_mem, mk_afdesign_model
    from colabdesign.mpnn import mk_mpnn_model
    from colabdesign.shared.utils import copy_dict
    print("  colabdesign AF + MPNN + utils OK")
except Exception as e:
    print(f"  FAILED  {type(e).__name__}: {e}")

print("\n--- jax device ---")
import jax
print("  jax", jax.__version__, jax.devices())
PYEOF

# BindCraft's own top-level imports, exercised through its own module rather
# than our list of guesses. A module we never thought to check is the failure
# this catches.
echo "--- BindCraft's own imports ---"
cd "$ROOT/BindCraft"
"$PY" -c "
import sys, traceback
sys.path.insert(0, '.')
try:
    from functions import *
    print('  functions package OK')
except Exception as e:
    print(f'  FAILED  {type(e).__name__}: {e}')
    traceback.print_exc(limit=3)
" 2>&1 | tail -20

echo "=== deps job done ==="
exit 0
EOS

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL \
     bash /workspace/rfd/_bc_deps.sh
