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
# First attempt reused the container torch (--system-site-packages + torch pin), but
# BoltzGen pulls numpy 2.0.2 and the container torch (NGC 2.7.0) was built against
# numpy 1.x -> "import torch" crashes. So build a SELF-CONTAINED venv: no
# system-site-packages, discard the container pip constraint entirely (PIP_CONSTRAINT
# off), and let BoltzGen install its own consistent torch+numpy2+numba stack from PyPI
# (reachable from compute). Costs a ~GB torch download; thats fine on /home.
echo "=== recreate venv WITHOUT system-site-packages ==="
rm -rf "$VENV"
python -m venv "$VENV"
source "$VENV/bin/activate"
pip install --upgrade pip 2>&1 | tail -1
echo "=== install boltzgen with its own stack (container pip constraint disabled) ==="
export PIP_CONSTRAINT=/dev/null
pip install --no-cache-dir boltzgen 2>&1 | tail -25
echo "=== verify (no GPU on this env-build node, so cuda avail=False is expected) ==="
which boltzgen
python -c "import numpy,torch;print(\"numpy\",numpy.__version__,\"torch\",torch.__version__,\"cuda\",torch.version.cuda,\"avail\",torch.cuda.is_available())"
boltzgen --help 2>&1 | head -15
echo "=== done ==="
'
