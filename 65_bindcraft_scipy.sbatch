#!/bin/bash
# BindCraft step 6: downgrade SciPy inside venv_bc. Short job, one GPU.
#
#   cd ~/rfd && sbatch scripts/65_bindcraft_scipy.sbatch
#
# What job 176xxx (64) reported looked like two independent conflicts:
#
#   colabdesign  AttributeError: module 'scipy.linalg' has no attribute 'tril'
#   optax        ImportError: cannot import name 'stats' from 'jax.scipy'
#
# They are ONE cause. SciPy deprecated scipy.linalg.tril / triu in 1.11 and
# REMOVED them in 1.13; the container ships 1.15.2. jax 0.4.20 predates that
# removal and calls scipy.linalg.tril while importing jax.scipy.linalg. That
# import dies, leaving jax.scipy half-initialised, so the next thing to ask for
# jax.scipy.stats - optax - gets "cannot import name 'stats'". The second error
# is a symptom of the first, and fixing SciPy should clear both.
#
# Note what is NOT the problem: ColabDesign. We suspected its "fix for newer
# JAX" commit would not work against our old jax, and it is innocent here - the
# traceback goes through colabdesign into jax into scipy.
#
# The fix direction matters. Moving jax forward is the obvious instinct and it
# is blocked: the only GPU wheel we have is jaxlib 0.4.20+cuda12.cudnn89 and the
# cluster has cuDNN 8.9, so a newer jax means no GPU. So move SciPy back
# instead. 1.12.0 is the last line that still carries tril/triu (deprecated, not
# removed) and still ships cp312 wheels. If 1.12.0 somehow misses something,
# 1.11.4 is the contemporary of jax 0.4.20 and the fallback.
#
# This only touches venv_bc. venv_af2 keeps the container's SciPy, so stage-3
# scoring is unaffected either way.

#SBATCH --job-name=bc-scipy
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
SCIPY_PIN="${SCIPY_PIN:-1.12.0}"
mkdir -p "$HOME/rfd/logs"
INNER="$HOME/rfd/_bc_scipy.sh"

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

echo "--- before ---"
"$PY" -c "import scipy, scipy.linalg; print('scipy', scipy.__version__, 'at', scipy.__file__); print('has linalg.tril:', hasattr(scipy.linalg,'tril'))"

"$PY" -c "import torch;print('torch=='+torch.__version__)" > /workspace/rfd/pipc_bc.txt
export PIP_CONSTRAINT=/workspace/rfd/pipc_bc.txt

# --no-deps so this cannot drag numpy with it: numpy must stay at 1.26.4,
# because jax 0.4.20 dies at import on numpy 2 (np.trapz, removed there).
echo "--- installing scipy==$SCIPY_PIN into venv_bc only ---"
"$PY" -m pip install --no-deps "scipy==$SCIPY_PIN" 2>&1 | tail -5

echo "--- after ---"
# Check the version AND that the attribute is actually back. A version number
# is what we asked for; the attribute is what jax needs.
"$PY" -c "
import scipy, scipy.linalg, numpy
print('scipy', scipy.__version__, 'at', scipy.__file__)
print('has linalg.tril:', hasattr(scipy.linalg,'tril'))
print('numpy', numpy.__version__, '(must still be 1.26.4)')
"

echo "--- the three imports that failed before ---"
"$PY" - <<'PYEOF'
import warnings
warnings.filterwarnings("ignore")   # scipy.linalg.tril is deprecated here, and
                                    # the DeprecationWarning is not news
for label, fn in [
    ("jax.scipy.linalg",            lambda: __import__("jax.scipy.linalg", fromlist=["x"])),
    ("jax.scipy.stats",             lambda: __import__("jax.scipy.stats",  fromlist=["x"])),
    ("optax",                       lambda: __import__("optax")),
    ("colabdesign",                 lambda: __import__("colabdesign")),
]:
    try:
        fn(); print(f"  {label:22s} OK")
    except Exception as e:
        print(f"  {label:22s} FAILED  {type(e).__name__}: {e}")

print("\n--- the entry points ---")
try:
    from colabdesign import clear_mem, mk_afdesign_model
    from colabdesign.mpnn import mk_mpnn_model
    print("  colabdesign AF + MPNN OK")
except Exception as e:
    print(f"  FAILED  {type(e).__name__}: {e}")

print("\n--- jax device (scipy must not have cost us the GPU) ---")
import jax
print("  jax", jax.__version__, jax.devices())
PYEOF

echo "--- BindCraft's own imports ---"
cd "$ROOT/BindCraft"
"$PY" -c "
import sys, traceback, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '.')
try:
    from functions import *
    print('  functions package OK - the environment is done')
except Exception as e:
    print(f'  FAILED  {type(e).__name__}: {e}')
    traceback.print_exc(limit=4)
" 2>&1 | tail -20

echo "=== scipy job done ==="
exit 0
EOS

export SCIPY_PIN
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     --export=ALL \
     bash /workspace/rfd/_bc_scipy.sh
