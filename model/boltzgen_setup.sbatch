#!/bin/bash
# Install BoltzGen (MIT, Hannes Stark / MIT J-Clinic) into its OWN venv on the
# Thapar cluster, mirroring the Boltz install that worked (job 171581): own venv
# with --system-site-packages so it reuses the container's torch (no 2.5 GB
# re-download), and DISCARD the container's pip constraint file keeping only a torch
# pin so BoltzGen's deps resolve instead of hitting ResolutionImpossible.
# No GPU needed for install -> no --gres (queues immediately). Network: pip hits
# PyPI (reachable from compute, per the Boltz precedent). Model WEIGHTS are a
# separate problem (HuggingFace, ~40 GB) handled after we confirm reachability.
#
#   cd ~/rfd/scripts && git checkout claude/project-thread-ho43or && git pull \
#     && sbatch model/boltzgen_setup.sbatch
# Read the .out: resolved pip plan, then `boltzgen --help` + torch version if OK.

#SBATCH --job-name=bgsetup
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --exclude=dgxanode03
#SBATCH --time=01:00:00
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
echo "=== create venv (reuses container torch) ==="
python -m venv --system-site-packages "$VENV"
source "$VENV/bin/activate"
python -c "import torch;print(\"container torch==\"+torch.__version__)"
python -c "import torch;print(\"torch==\"+torch.__version__)" > "$ROOT/pipc_bg.txt"
echo "=== DRY-RUN install plan (catches torch/constraint conflicts) ==="
PIP_CONSTRAINT="$ROOT/pipc_bg.txt" pip install --dry-run boltzgen 2>&1 | tail -40
echo "=== REAL install ==="
PIP_CONSTRAINT="$ROOT/pipc_bg.txt" pip install --no-cache-dir boltzgen 2>&1 | tail -25
echo "=== verify ==="
which boltzgen && boltzgen --help 2>&1 | head -15
python -c "import torch;print(\"final torch==\"+torch.__version__,\"cuda\",torch.version.cuda)"
echo "=== done ==="
'
