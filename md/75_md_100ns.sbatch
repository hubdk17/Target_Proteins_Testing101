#!/bin/bash
#SBATCH --job-name=md100
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=32G
#
# 100 ns MD of one complex, build included, in 10 ns resumable chunks.
# Submit FROM your MD clone, e.g. ~/rfd/scripts_md (sbatch opens md/...
# relative to cwd, and every script resolves its own location):
#
#   sbatch --gres=gpu:1 --time=24:00:00 --exclude=dgxanode03,dgxhnode6 \
#          --job-name=md100_p53 --output=logs/%x_%j.out --error=logs/%x_%j.err \
#          md/75_md_100ns.sbatch mdm2_blind/ref_p53.pdb runs100/p53_mdm2 10
#
# args: <complex.pdb> <outdir> [nchunks, 10 ns each; default 10 = 100 ns]
# Re-submitting the identical line after any failure resumes from the last
# finished chunk.  dgxanode03 is excluded because it reports an NVML driver/
# library version mismatch (no CUDA job can run there); dgxhnode6 because the
# prebuilt pmemd.cuda may not carry H100 kernels.
PDB="$1"; OUT="$2"; NCHUNK="${3:-10}"
IMG=/home/hpc_apps/nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh
# The inner script is read on the COMPUTE NODE at runtime, out of the working
# tree you submitted from -- so hardcoding a clone path lets another
# workstream's `git checkout` in a shared clone delete it between submit and
# execution (job 193870, exit 127). Resolve it from the submit directory
# instead, mapped into the container ($HOME is mounted at /workspace).
SUB="${SLURM_SUBMIT_DIR:-$PWD}"
CSUB="/workspace${SUB#$HOME}"
echo "[submit] $SUB -> container $CSUB"
srun --container-image=$IMG --container-mounts=$HOME:/workspace \
     bash "$CSUB/md/75_md_100ns_inner.sh" "$PDB" "$OUT" "$NCHUNK"
