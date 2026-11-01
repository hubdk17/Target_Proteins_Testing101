#!/bin/bash
# BindCraft step 3: install ColabDesign into venv_bc and find out what is
# actually missing. Short job, one GPU, nothing downloaded.
#
#   cd ~/rfd && sbatch scripts/62_bindcraft_imports.sbatch
#   # then, when it finishes:
#   cat logs/bc-imports_<jobid>.out
#
# Why a GPU for an import check: it is a few minutes, and it lets the same job
# answer "does jax see the device from venv_bc" instead of leaving that for a
# later job to discover.
#
# Two things this job does NOT do: touch the network (compute nodes are
# host-filtered, so ColabDesign installs from the clone 61 made), and touch
# venv_af2.

#SBATCH --job-name=bc-imports
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
# dgxanode03's container disk is intermittently full: enroot dies unpacking the
# squashfs, FAILED in under 2 minutes with a zero-byte log. 23 failures to 2
# successes there, none elsewhere.
#SBATCH --exclude=dgxanode03
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# The container venv_af2 was built inside, so its pins match what is here.
IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh

mkdir -p "$HOME/rfd/logs"
INNER="$HOME/rfd/_bc_imports.sh"

cat > "$INNER" <<'EOS'
set -uo pipefail          # not -e: a failed import is the RESULT of this job
ROOT=/workspace/rfd
BC="$ROOT/BindCraft"
PY="$ROOT/venv_bc/bin/python"

echo "node: $(hostname)"
nvidia-smi --query-gpu=name --format=csv,noheader

# --- the path trap, fixed before anything else -----------------------------
# A venv records its own absolute path in bin/activate and in every console
# script shebang. venv_af2 was built INSIDE this container, so what it records
# is /workspace/rfd/venv_af2 - not the /home/... path the login-node copy
# rewrote. So the copy's scripts can still point at venv_af2, which is the one
# env we are trying not to disturb.
#
# Two defences. First, rewrite the container-side path too:
if grep -rqs --binary-files=without-match "$ROOT/venv_af2" "$ROOT/venv_bc/bin"; then
  echo "fixing container-side paths in venv_bc/bin"
  grep -rls --binary-files=without-match "$ROOT/venv_af2" "$ROOT/venv_bc/bin" \
    | while read -r f; do sed -i "s|$ROOT/venv_af2|$ROOT/venv_bc|g" "$f"; done
fi
echo "venv_af2 references left in venv_bc/bin (want 0): $(grep -rls --binary-files=without-match "$ROOT/venv_af2" "$ROOT/venv_bc/bin" 2>/dev/null | wc -l)"

# Second, and the one that actually guarantees it: never source activate. The
# python binary resolves its own prefix from its location, so calling it by
# path cannot land in the wrong env whatever the scripts say.
echo "--- which env are we in ---"
"$PY" -c 'import sys; print("prefix:", sys.prefix); print("exe:", sys.executable)'
PREFIX=$("$PY" -c 'import sys; print(sys.prefix)')
case "$PREFIX" in
  *venv_bc) echo "OK: operating on venv_bc" ;;
  *) echo "*** WRONG ENV ($PREFIX) - refusing to install. Stopping."; exit 1 ;;
esac

# --- install ColabDesign, and nothing else ---------------------------------
# --no-deps so pip cannot resolve away the pins this env was built around
# (jax 0.4.20 against jaxlib 0.4.20+cuda12.cudnn89, dm-haiku 0.0.13 because
# 0.0.14 reads a jax attribute 0.4.20 lacks, numpy 1.26.4 because jax 0.4.20
# dies on numpy 2). --no-build-isolation so the build needs no network.
echo "--- installing colabdesign (local clone, no deps, no network) ---"
"$PY" -m pip install --no-deps --no-build-isolation --no-index "$ROOT/colabdesign" 2>&1 | tail -5

# --- what is present, what is missing --------------------------------------
# Each import is attempted separately and reported. Nothing aborts the sweep,
# because the point is the full list in one job rather than one name per job.
echo "--- import sweep ---"
"$PY" - <<'PYEOF'
import importlib, sys
mods = [
    # the load-bearing four
    "jax", "jaxlib", "haiku", "pyrosetta",
    # ColabDesign and what it leans on
    "colabdesign", "optax", "chex", "ml_collections", "tree", "immutabledict",
    "joblib", "py3Dmol",
    # scientific base BindCraft's own scripts use
    "numpy", "scipy", "pandas", "matplotlib", "seaborn", "Bio", "tqdm",
]
ok, missing, broken = [], [], []
for m in mods:
    try:
        mod = importlib.import_module(m)
        ok.append((m, getattr(mod, "__version__", "?")))
    except ModuleNotFoundError:
        missing.append(m)
    except Exception as e:
        broken.append((m, f"{type(e).__name__}: {e}"))

print("\nPRESENT:")
for m, v in ok:
    print(f"  {m:16s} {v}")
print("\nMISSING (not installed):")
print("  none" if not missing else "\n".join(f"  {m}" for m in missing))
print("\nBROKEN (installed but will not import - version conflict):")
print("  none" if not broken else "\n".join(f"  {m}\n      {e}" for m, e in broken))

# The real entry point. ColabDesign's HEAD carries a commit titled "fix for
# newer JAX", so whether it imports against our deliberately old jax 0.4.20 is
# the open question of this whole job, not a formality.
print("\n--- the entry point BindCraft calls ---")
try:
    from colabdesign import mk_afdesign_model
    print("  mk_afdesign_model OK")
except Exception as e:
    print(f"  FAILED  {type(e).__name__}: {e}")

print("\n--- does jax see the GPU ---")
try:
    import jax
    print("  jax", jax.__version__, "devices:", jax.devices())
except Exception as e:
    print(f"  FAILED  {type(e).__name__}: {e}")
PYEOF

# --- are the weights where BindCraft looks for them ------------------------
echo "--- parameters ---"
ls -1 "$BC/params" 2>/dev/null | grep -c multimer | sed 's/^/  multimer files: /'

echo "=== import check done ==="
exit 0
EOS

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL \
     bash /workspace/rfd/_bc_imports.sh
