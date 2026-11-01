# RFdiffusion on Thapar AI Foundry — trypsin backbone run

Rewritten 2026-10-16 to follow the platform's documented job pattern (Enroot
containers from `/home/hpc_apps`, home mounted at `/workspace`) rather than a
generic conda install. The earlier conda-based scripts are gone.

## Run it — four small steps, cheapest first

Each step is checkable on its own, and only two of them touch a GPU at all.

```bash
# on the login node (bcmnode1)
module load slurm
mkdir -p ~/rfd/logs
cp 01_fetch.sh rfd_job.sh 10_hello_gpu.sbatch 20_buildenv.sbatch 30_rfd.sbatch ~/rfd/
cd ~/rfd

sbatch 10_hello_gpu.sbatch      # 1. ~30 s of GPU: container loads, torch sees GPU
cat logs/hello_*.out            #    read this before going on

bash 01_fetch.sh                # 2. login node, no GPU: source + weights + trypsin

sbatch 20_buildenv.sbatch       # 3. NO GPU at all: pip installs, 10-20 min
tail -f logs/buildenv_*.out

sbatch 30_rfd.sbatch            # 4. a few minutes of GPU: 1 backbone, 50 aa
ls outputs/
```

Then the trypsin run, same small size:

```bash
MODE=binder sbatch 30_rfd.sbatch
```

And only once that works, scale up:

```bash
NUM_DESIGNS=10 LENGTH=100 DIFFUSER_T=50 sbatch 30_rfd.sbatch   # raise --time too
```

Logs are in `~/rfd/logs/`, backbones in `~/rfd/outputs/`.

## The five files

| file | GPU? | what it does |
|---|---|---|
| `10_hello_gpu.sbatch` | ~30 s | proves the container loads, home is mounted at `/workspace`, torch sees the GPU, and whether the compute node can reach PyPI |
| `01_fetch.sh` | no (login node) | clones RFdiffusion, pulls the two needed checkpoints (~10 GB), cleans PDB 2PTN, prints its residue numbering |
| `20_buildenv.sbatch` | **none** | builds the venv inside the container — no `--gres` line, so pip never holds a GPU |
| `30_rfd.sbatch` | few min | the real run: 1 design, 50 aa, 20 timesteps, 20-min cap |
| `rfd_job.sh` | — | the body that runs inside the container for steps 3 and 4 |

Deliberately small defaults in `30_rfd.sbatch`: `NUM_DESIGNS=1`, `LENGTH=50`,
`DIFFUSER_T=20` (RFdiffusion's own default is 50 timesteps). Cost scales as
designs × timesteps × residues², so these three knobs are where GPU time goes.

## What the platform documentation actually specifies

<https://aifoundry.thapar.edu/documentation/> has five sub-pages — SLURM,
PyTorch/TensorFlow, Gromacs, Quantum Espresso, Amber. The PyTorch page is the
one that matters; it defines the house style for a GPU job, and these scripts
follow it:

- **Containers, not conda.** Jobs run in Enroot images with
  `--container-image=/home/hpc_apps/<image>.sqsh`. Available images:
  `nvidia+pytorch+25.04-py3.sqsh`, `nvidia+tensorflow+25.02-tf2-py3.sqsh`, and
  (from the Gromacs page) `nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh`.
  The docs say users don't upload or build their own.
- **Always mount home at `/workspace`** (`--container-mounts=/home/$USER:/workspace`)
  — the docs warn about permission errors otherwise.
- **Resource shape** from the documented example: `--nodes=1 --ntasks-per-node=1
  --cpus-per-task=20 --mem=64G --gres=gpu:1`, optional `--nodelist=dgxhnode6`.
- **Interactive debugging** is `srun --pty ... --container-image=... /bin/bash`,
  which is how to poke at a failure instead of resubmitting blind.
- `module load slurm` first, then `sinfo` / `squeue` / `scancel <jobid>`.

Confirmed from this account on 2026-10-16: `defq` is the only partition and the
default; GRES is `gpu:8` so `--gres=gpu:1` is correct syntax; nodes are
`dgxanode03`, `dgxanode04` (A100) and `dgxhnode6` (H100) — no B200, so Blackwell
kernel support is not a concern. `/home` is a 97 TB share with 12 TB free and no
per-user quota visible, so the ~10 GB of weights is not a problem. `enroot` is at
`/usr/bin/enroot`; no `apptainer`, `singularity` or `docker`. There is no conda
module — hence the venv-inside-container approach.

## Where this can still fail

1. **`dgl`** — the real risk. RFdiffusion's SE3Transformer is built on DGL, and
   DGL ships wheels per torch/CUDA combination, lagging new torch releases. The
   PyTorch 25.04 image is torch 2.7-era, newer than DGL's published wheels.
   `rfd_job.sh` tries three wheel indexes and stops cleanly with a pointer here
   if all fail, rather than dying later in a confusing place.

   **If dgl fails**, don't patch the script — change the container. Two options,
   in order of preference:

   ```bash
   # a) import an older NGC PyTorch image yourself (enroot is available)
   cd /home/$USER/rfd
   enroot import docker://nvcr.io/nvidia/pytorch:23.10-py3   # torch 2.1, dgl wheels exist
   # produces nvidia+pytorch+23.10-py3.sqsh -- point IMAGE= at it in the sbatch
   ```

   b) Ask the AI Foundry admins (same ticket system) to add an RFdiffusion image
   to `/home/hpc_apps`. The upstream repo ships `docker/Dockerfile`, so this is a
   routine request, and it's the right long-term answer if the group will use
   RFdiffusion repeatedly.

2. **Compute-node internet.** The venv build runs `pip install` from inside the
   job. The login node clearly has outbound access; if the compute nodes don't,
   the env build fails on the first `pip install`. Check it with the documented
   interactive pattern before burning a queue slot:

   ```bash
   srun --pty --gres=gpu:1 --cpus-per-task=20 --mem=64G \
     --container-image=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh \
     --container-mounts=/home/$USER:/workspace \
     --container-workdir=/workspace /bin/bash
   # then: pip download hydra-core -d /tmp/x    (does it reach PyPI?)
   ```

3. **`files.ipd.uw.edu`** (model weights) must be reachable from the login node.
   If `01_fetch.sh` stalls there, download the two checkpoints elsewhere and scp
   them into `/home/$USER/rfd/models`.

4. **No `--account` flag** is set because none is documented. If submission is
   rejected for a missing account, add it.

## On "trypsin without specificity"

I read this as no hotspot residues specified, and the scripts give both readings
because they are very different runs:

- **`MODE=unconditional`** (default) is a true smoke test. A 100-residue backbone
  is generated from noise with no target at all — **nothing about the output
  relates to trypsin**. Its only job is proving container, GPU, weights and
  install work end to end.
- **`MODE=binder`** conditions on the trypsin structure but passes no
  `ppi.hotspot_res`, so the model picks its own surface site. This is the real
  "trypsin, unspecified site" run. Expect designs scattered over the surface and
  many of them poor; RFdiffusion's own docs treat 3–6 hotspots as the normal way
  to steer PPI design.

**Numbering caveat.** 2PTN (bovine β-trypsin) uses chymotrypsinogen numbering:
chain A runs about 16–245, not 1–223, with gaps and insertion codes. The contig
`A16-245` in the sbatch is a default, not a verified fact — `01_fetch.sh` prints
the actual first and last residue and any gaps it finds, so set `TARGET_CONTIG`
from that output. If there's a gap mid-chain, split the contig across it
(`TARGET_CONTIG="A16-120/A125-245"`) rather than spanning it.

Whole-target binder design is also slower than it needs to be — RFdiffusion
scales as O(N²) in residue count. Once the smoke test passes, trimming trypsin to
~10 Å around the S1 pocket is the standard next move, and it's also where
hotspots would go.

## After a backbone exists

Outputs are `*.pdb` plus a `*.trb` record per design. Stage two of the pipeline is
ProteinMPNN sequence design on those backbones; stage three is the AF2 or
Chai/Boltz re-prediction filter scored on interface ipTM/PAE. Both can reuse the
same container pattern.
