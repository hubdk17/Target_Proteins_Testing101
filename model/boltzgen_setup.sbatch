#!/bin/bash
# Install BoltzGen (MIT, Hannes Stark / MIT J-Clinic) into its own venv on the Thapar
# cluster. Mirrors the WORKING AF2 install recipe: own venv with
# --system-site-packages (reuse the container torch, tiny download) + a pip constraint
# pinning the container's torch AND numpy<2.
# History of this file:
#  (1) torch-pin only -> boltzgen upgraded numpy to 2.0.2, which broke the numpy-1-built
#      container torch ("import torch" crash).
#  (2) fresh self-contained venv -> "Network is unreachable" pulling the ~GB torch/CUDA
#      stack (compute network flaky over IPv6).
#  (3, this) reuse container torch + PIN numpy to the container's 1.26.x so boltzgen's
#      deps resolve against numpy 1.x and numpy vanishes from the install list (proven
#      with AF2 / jax). Small download, no ABI break.
# No GPU needed for install -> no --gres (queues immediately). pip hits PyPI (reachable
# from compute). Model WEIGHTS (HuggingFace ~40 GB) are staged separately from login.
#
#   cd ~/rfd/scripts && git checkout claude/project-thread-ho43or && git pull \
#     && sbatch model/boltzgen_setup.sbatch
# Read the .out: the constraint, the install tail, then a `numpy 1.26 torch 2.x` line
# and `boltzgen --help` with NO traceback.

#SBATCH --job-name=bgsetup
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --exclude=dgxanode03
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

IMAGE=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh

echo "node: $(hostname)  job: $SLURM_JOB_ID"
srun --container-image="$IMAGE" \
     --container-mounts="/home/$USER:/workspace" \
     --container-workdir=/workspace \
     bash -lc '
set -uo pipefail
ROOT=/workspace/rfd
VENV=$ROOT/venv_boltzgen
python --version
echo "=== recreate venv reusing container torch ==="
rm -rf "$VENV"
python -m venv --system-site-packages "$VENV"
source "$VENV/bin/activate"
python -c "import torch;print(\"torch==\"+torch.__version__)"  > "$ROOT/pipc_bg.txt"
python -c "import numpy;print(\"numpy==\"+numpy.__version__)" >> "$ROOT/pipc_bg.txt"
echo "--- pip constraint (keep container torch + numpy<2):"; cat "$ROOT/pipc_bg.txt"
echo "=== install boltzgen (numpy kept at container 1.x) ==="
PIP_CONSTRAINT="$ROOT/pipc_bg.txt" pip install --no-cache-dir --retries 6 --timeout 90 boltzgen 2>&1 | tail -30
echo "=== verify (import must succeed; no GPU on this node so cuda-avail is irrelevant) ==="
which boltzgen
python -c "import numpy,torch;print(\"numpy\",numpy.__version__,\"torch\",torch.__version__,\"cuda\",torch.version.cuda)"
boltzgen --help 2>&1 | head -15
echo "=== done ==="
'
