#!/bin/bash
# STEP 3 -- build the Python environment. NO GPU: there is deliberately no
# --gres line, so the 10-20 minutes of pip downloading never occupies a GPU.
#
#   sbatch 20_buildenv.sbatch
#   tail -f logs/buildenv_*.out
#
# The venv lands at ~/rfd/venv, which is on /home, so it survives the container
# and every later job reuses it. Self-contained: writes its own inner script.

#SBATCH --job-name=rfd-env
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=logs/buildenv_%j.out
#SBATCH --error=logs/buildenv_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
echo "node: $(hostname)  job: $SLURM_JOB_ID  (no GPU requested)"
test -d "$HOME/rfd/RFdiffusion" || { echo "run 01_fetch.sh first"; exit 1; }

# The payload is written to a file rather than passed as a quoted string, so
# nothing here depends on shell quoting surviving two layers.
cat > "$HOME/rfd/_buildenv_inner.sh" <<'INNER'
set -euo pipefail
ROOT=/workspace/rfd
VENV=$ROOT/venv

echo "python: $(python --version)"
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.version.cuda)"

python -m venv --system-site-packages "$VENV"
source "$VENV/bin/activate"
python -m pip install --quiet --upgrade pip

echo "=== 1/3 straightforward dependencies ==="
pip install hydra-core omegaconf icecream e3nn opt_einsum pyrsistent scipy

echo "=== 2/3 dgl (the risky one) ==="
# RFdiffusion's SE3Transformer is built on DGL. DGL publishes wheels per
# torch/CUDA combination and lags new torch releases, and this container is
# torch 2.7 / cu129 / py3.12 -- newer than anything DGL has published. We try
# the plausible wheel indexes, then verify by IMPORTING, because an ABI
# mismatch installs fine and only fails at import.
DGL_OK=0
for IDX in \
  "https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html" \
  "https://data.dgl.ai/wheels/torch-2.3/cu121/repo.html" \
  "https://data.dgl.ai/wheels/repo.html" \
  "PYPI"; do
  echo "--- trying: $IDX"
  if [ "$IDX" = "PYPI" ]; then
    pip install dgl || continue
  else
    pip install dgl -f "$IDX" || continue
  fi
  if python -c "import dgl; print('dgl', dgl.__version__, 'imports cleanly')" 2>&1; then
    DGL_OK=1; break
  else
    echo "--- installed but failed to import; uninstalling and trying next"
    pip uninstall -y dgl >/dev/null 2>&1 || true
  fi
done

if [ "$DGL_OK" -eq 0 ]; then
  echo
  echo "############################################################"
  echo "# dgl will not work in this container (torch 2.7 / py3.12). #"
  echo "# This is the known hard edge. No GPU time was spent.       #"
  echo "# Next step is an older container -- tell Claude.           #"
  echo "############################################################"
  exit 3
fi

echo "=== 3/3 SE3Transformer + RFdiffusion ==="
pip install --no-deps --no-build-isolation "$ROOT/RFdiffusion/env/SE3Transformer"
pip install --no-deps --no-build-isolation -e "$ROOT/RFdiffusion"

python -c "import dgl, e3nn, hydra, torch, rfdiffusion; print('ALL IMPORTS OK')"
touch "$VENV/.complete"
echo "=== environment built: $VENV ==="
INNER

srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash /workspace/rfd/_buildenv_inner.sh
