#!/bin/bash
#SBATCH --job-name=md100
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=32G
#
# 100 ns MD of one complex, build included, in 10 ns resumable chunks.
# Submit FROM ~/rfd/scripts (sbatch opens md/... relative to cwd):
#
#   sbatch --gres=gpu:1 --time=24:00:00 --exclude=dgxhnode6 \
#          --job-name=md100_p53 --output=logs/%x_%j.out --error=logs/%x_%j.err \
#          md/75_md_100ns.sbatch mdm2_blind/ref_p53.pdb runs100/p53_mdm2 10
#
# args: <complex.pdb> <outdir> [nchunks, 10 ns each; default 10 = 100 ns]
# Re-submitting the identical line after any failure resumes from the last
# finished chunk.  dgxhnode6 is excluded because the prebuilt pmemd.cuda may
# not carry H100 kernels; drop the exclude once it is seen to work there.
PDB="$1"; OUT="$2"; NCHUNK="${3:-10}"
IMG=/home/hpc_apps/nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh
srun --container-image=$IMG --container-mounts=/home/dakshcsed:/workspace \
     bash /workspace/rfd/scripts/md/75_md_100ns_inner.sh "$PDB" "$OUT" "$NCHUNK"
