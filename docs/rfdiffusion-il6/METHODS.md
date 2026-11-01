# Methods — how to run each stage, and how to measure it

Everything here has been run successfully on the Thapar AI Foundry cluster.
Index: `README.md`. Results: `STAGE3-boltz-scoring.md` and `RESULTS.md`.

Live working directory on the cluster is `~/rfd/` (not this shared folder):
`venv/`, `venv_boltz/`, `RFdiffusion/`, `ProteinMPNN/`, `models/`, `inputs/`,
`logs/`, `outputs/`, `mpnn_in/`, `mpnn_out/`, `boltz_in/`, `boltz_out/`,
`boltzhome/`.

## Cluster basics

```bash
module load slurm          # required on every fresh login, or sbatch is not found
```

One partition, `defq`, no time limit. Nodes `dgxanode03`, `dgxanode04` (A100),
`dgxhnode6` (H100), GRES `gpu:8`, so `--gres=gpu:1`.

There is **no conda and no CUDA module**. The pattern is Enroot containers from
`/home/hpc_apps/`, with home mounted at `/workspace`:

```
--container-image=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh
--container-mounts=/home/$USER:/workspace
--container-workdir=/workspace
```

Python packages live in venvs created inside the container with
`--system-site-packages`, so they reuse the container's torch rather than
downloading 2.5 GB of their own.

**Use `sbatch`, never interactive `srun`, for pipeline steps.** An SSH drop kills
an interactive job and leaves a truncated output file.

Three traps that each cost a job:

- `dgxanode03` can fail with `pyxis: failed to create container filesystem` /
  `No space left on device`. That is the **node's local disk**, where Enroot
  unpacks the 20 GB image. Not user-cleanable. Resubmit with
  `--exclude=dgxanode03`.
- The container sets `PIP_CONSTRAINT`, pinning every package it ships, which
  makes most `pip install`s fail with `ResolutionImpossible`. Do not patch the
  constraint file line by line. Discard it and keep only a torch pin:
  ```bash
  python -c "import torch;print('torch=='+torch.__version__)" > pipc.txt
  PIP_CONSTRAINT=pipc.txt pip install --no-cache-dir <package>
  ```
- **Compute nodes reach PyPI but cannot resolve `huggingface.co`.** Download
  model weights on the login node, into the cache directory the job will read,
  then run the job offline. Assume any new tool cannot fetch its own weights.

## Stage 1 — RFdiffusion backbones

Script `30_rfd_il6.sbatch`, parameterised by environment variables so the same
file serves a smoke test and a real batch:

```bash
cd ~/rfd && NUM_DESIGNS=10 DIFFUSER_T=50 BINDER_LEN=40-55 sbatch 30_rfd_il6.sbatch
```

The inference call it builds:

```bash
python scripts/run_inference.py \
  inference.output_prefix=$ROOT/outputs/il6_binder_site3 \
  inference.model_directory_path=$ROOT/models \
  inference.input_pdb=$ROOT/inputs/il6_clean.pdb \
  inference.ckpt_override_path=$ROOT/models/Complex_base_ckpt.pt \
  inference.num_designs=10 \
  diffuser.T=50 \
  'contigmap.contigs=[A19-51/A61-184/0 40-55]' \
  'ppi.hotspot_res=[A47,A155,A157,A158]'
```

- `A19-51/A61-184` is two segments because **52–60 are genuinely absent from the
  1ALU coordinates**. One range would ask the model to bridge a gap that is not
  there.
- `/0 40-55` is the chain break followed by the binder length range.
- `ppi.hotspot_res` is what steers the binder to site III. Omitting it is legal
  and the "no hotspots given" warning is harmless.
- DGL is the one real environment trap: DGL 2.4 dies on import with
  `No module named 'torchdata.datapipes'`. **Use `dgl==2.5.0+cu121`** and set
  `DGLBACKEND=pytorch`.

Output is **two chains and no sequence**: A is the target and does not move
(motif RMSD 0.10), B is the binder numbered from 1, every residue a glycine
placeholder. About one minute per design.

**Give every run a distinct `output_prefix`.** Two jobs submitted with the same
prefix silently overwrote each other.

**Better than filtering afterwards:** enforce compactness during diffusion.

```
potentials.guiding_potentials=["type:binder_ROG,weight:3"]
potentials.guide_scale=2
potentials.guide_decay="quadratic"
```

## Stage 2 — ProteinMPNN sequences

Plain clone, no environment build needed: it ships its weights in
`vanilla_model_weights/` and needs only PyTorch.

Runner `run_mpnn.sh`, three steps: parse the structures, mark chain A as fixed
and chain B as designable, then design. **The header comment in the saved copy
is out of date on one point: request a GPU.** ProteinMPNN takes about 11 minutes
per structure on CPU in this container and about 5 seconds on a GPU, so submit
with `--gres=gpu:1`.

The chain assignment step is not optional. Without it ProteinMPNN rewrites the
target's sequence too:

```bash
python helper_scripts/assign_fixed_chains.py \
  --input_path=$ROOT/mpnn_out/parsed.jsonl \
  --output_path=$ROOT/mpnn_out/assigned.jsonl \
  --chain_list "B"
```

Design settings used: `--num_seq_per_target 8 --sampling_temp 0.1 --seed 37
--batch_size 1`. Low temperature is the usual choice for binder design.

**Reading the output:**

- Score is lower-is-better. `short_9`'s eight samples came in at 1.30–1.37.
- **`seq_recovery` is meaningless here** and will be ~0. The "native" sequence of
  an RFdiffusion backbone is poly-glycine, so there is nothing to recover.
- **Poly-alanine output is a diagnosis, not a bad roll.** It is the signature of
  a backbone with no hydrophobic core, and it means the backbone should be
  rejected rather than redesigned.

## Stage 3 — Boltz co-folding

Boltz lives in its **own venv** `~/rfd/venv_boltz`, not the RFdiffusion one,
because about 29 packages overlap and a downgrade there could silently break
stage 1.

Input YAML, deliberately flow-style so it survives being pasted into a terminal:

```yaml
version: 1
sequences:
  - protein: {id: A, sequence: <target>, msa: empty}
  - protein: {id: B, sequence: <binder>, msa: empty}
```

Check the lengths before running, rather than after:

```bash
grep -o 'sequence: [A-Z]*' ~/rfd/boltz_in/short9_wt.yaml | awk '{print length($2)}'
```

The runner `run_boltz.sh`:

```bash
source /workspace/rfd/venv_boltz/bin/activate
export WANDB_MODE=offline HOME=/workspace/rfd/boltzhome
mkdir -p $HOME
boltz predict /workspace/rfd/boltz_in/short9_wt.yaml \
  --out_dir /workspace/rfd/boltz_out --output_format pdb \
  --diffusion_samples 50 --no_kernels --write_full_pae --override
```

Four things in there are load-bearing:

- **`HOME=/workspace/rfd/boltzhome`** so the `~/.boltz` weight and CCD cache
  lands on `/home` and persists, instead of on node-local disk where it is
  re-downloaded every job and can fill the filesystem.
- **`--no_kernels` is required in this container.** Boltz defaults to a CUDA
  kernel for the triangular multiplicative update that lives in the optional
  `cuequivariance_torch` package, which it does not install. Without the flag
  the run dies mid-prediction with `ModuleNotFoundError`.
- **`--diffusion_samples N` is how you get replicates**, because Boltz has no
  seed flag. It writes `model_0 … model_{N-1}`, each with its own confidence
  JSON. This is cheaper than N separate runs and is the right way to get a
  spread.
- `--override` is needed to rerun into an output directory that already has
  results, otherwise Boltz reuses the cache and reports "All inputs are already
  processed".

Submit line:

```bash
sbatch --job-name=boltzrun --partition=defq --cpus-per-task=8 --mem=48G \
  --gres=gpu:1 --time=01:30:00 --exclude=dgxanode03 \
  --output=logs/boltzrun_%j.out --error=logs/boltzrun_%j.err \
  --wrap="srun --container-image=/home/hpc_apps/nvidia+pytorch+25.04-py3.sqsh --container-mounts=/home/$USER:/workspace --container-workdir=/workspace bash /workspace/rfd/run_boltz.sh"
```

**Weights must be pre-staged on the login node** (compute nodes cannot resolve
the host). Into `~/rfd/boltzhome/.boltz/`, about 5.9 GB total:

```
https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar
https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt
https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_aff.ckpt
```

Use `wget -c` so an SSH drop resumes rather than restarts. Boltz unpacks
`mols.tar` itself, with no network.

**Output path.** Results land in
`out_dir/boltz_results_<input stem>/predictions/<name>/`, **not** directly under
`out_dir/predictions/`:

```bash
for f in boltz_out/*/predictions/*/confidence_*.json; do echo "--- $f"; cat $f; echo; done
```

Metrics in each confidence JSON: `confidence_score`, `ptm`, `iptm`,
`protein_iptm`, `complex_plddt`, `complex_iplddt`, `complex_pde`,
`complex_ipde`, `chains_ptm`, `pair_chains_iptm`. The full PAE matrix is **not**
in the JSON; it goes to `pae_<name>_model_N.npz` and needs `--write_full_pae`.
Per-token pLDDT goes to `plddt_*.npz`.

## Measurement rules

These are the rules that survived being wrong at least once.

### Compactness is a gate

Expected radius of gyration for a compact chain of N residues:

```
Rg_compact = 2.2 · N^0.38
```

Compute the ratio of the measured Rg to that, and **reject above about 1.3**
before looking at any interface metric. Contacts and compactness are
anti-correlated here, because a ribbon lying flat on the target touches more of
it, so ranking on contacts rewards the failure mode.

### Contact cutoffs on backbone-only output

- **CA–CA 10 Å** to survey the footprint.
- **CA–CA 8 Å** to separate grip from grazing.
- All-atom 5 Å is wrong on RFdiffusion output: there are no side chains.
- **Caveat that matters:** 10 Å CA–CA is wider than a helix is thick, so for a
  helix lying on a surface, residues on *both* faces fall inside the cutoff. It
  cannot tell you which face points at the target.

### The two design-triage measurements, both verified against known answers

Run from `~/rfd`. Each is a single line, deliberately, because a long pasted
heredoc mangled a job file once in this project.

**Compactness gate.** Radius of gyration of chain B, against the compact
expectation, per design:

```bash
for f in outputs/il6_n30_*.pdb; do awk 'BEGIN{n=0} substr($0,1,4)=="ATOM" && substr($0,13,4)==" CA " && substr($0,22,1)=="B" {n++; x[n]=substr($0,31,8)+0; y[n]=substr($0,39,8)+0; z[n]=substr($0,47,8)+0} END{if(n<2){printf "%-34s NO chain B CA atoms\n", FILENAME}else{for(i=1;i<=n;i++){sx+=x[i];sy+=y[i];sz+=z[i]} mx=sx/n;my=sy/n;mz=sz/n; for(i=1;i<=n;i++){s+=(x[i]-mx)^2+(y[i]-my)^2+(z[i]-mz)^2} rg=sqrt(s/n); e=2.2*exp(0.38*log(n)); printf "%-34s n=%-3d rg=%6.2f exp=%6.2f ratio=%.2f\n", FILENAME, n, rg, e, rg/e}}' $f; done
```

*Verification.* On four CA atoms at (0,0,0), (10,0,0), (0,10,0), (0,0,10) it
returns `n=4 rg=7.50 exp=3.73 ratio=2.01`, and 7.50 is the hand-computed answer
(mean at 2.5,2.5,2.5; summed square deviations 225; sqrt(225/4) = 7.5). Two
coincident atoms give rg 0.00. A chain A CA at (500,500,500) and a chain B CB at
(-500,-500,-500) are both correctly ignored. Cross-check: on `short_9` it should
reproduce the ratio 1.21 already recorded in `RESULTS.md`.

**Site III contact counter.** Per design: how many of the four conditioned
hotspots the binder reaches at 10 Å CA–CA, the closest approach to any of them,
and how many target residues are contacted in total. The last number separates
"grips site III" from "drapes over the whole protein".

```bash
for f in outputs/il6_n30_*.pdb; do awk 'BEGIN{na=0;nb=0;split("47 155 157 158",h," ");for(i in h)is3[h[i]]=1} substr($0,1,4)=="ATOM" && substr($0,13,4)==" CA " {ch=substr($0,22,1); r=substr($0,23,4)+0; px=substr($0,31,8)+0; py=substr($0,39,8)+0; pz=substr($0,47,8)+0; if(ch=="A"){na++; ax[na]=px; ay[na]=py; az[na]=pz; ar[na]=r} else if(ch=="B"){nb++; bx[nb]=px; by[nb]=py; bz[nb]=pz}} END{if(na<1||nb<1){printf "%-34s MISSING CHAIN na=%d nb=%d\n",FILENAME,na,nb; exit} tot=0; hot=0; minhot=9999; for(i=1;i<=na;i++){best=9999; for(j=1;j<=nb;j++){d=sqrt((ax[i]-bx[j])^2+(ay[i]-by[j])^2+(az[i]-bz[j])^2); if(d<best)best=d} if(best<=10)tot++; if(is3[ar[i]]){if(best<minhot)minhot=best; if(best<=10)hot++}} printf "%-34s nb=%-3d hotspots_hit=%d/4 min_hot=%5.2f target_contacts=%d\n",FILENAME,nb,hot,minhot,tot}' $f; done
```

*Verification.* On a synthetic file with one binder CA at the origin and target
CAs placed at 5.0, 9.9, 10.1 and 50.0 Å from it (A47, A155, A157, A158) plus
two non-hotspot residues at 3.0 and 200.0 Å, it returns
`nb=1 hotspots_hit=2/4 min_hot=5.00 target_contacts=3`, which is the
hand-counted answer. The 9.9/10.1 pair confirms the cutoff is inclusive at
exactly 10 Å, and a chain B CB decoy is ignored (`nb=1`, not 2).

Both carry an explicit `BEGIN{na=0;nb=0}`. That is not decoration: see
`CAUTION-origin-atom-bug.md`.

### Noise floor before any difference

Mean and SD straight off the confidence files. **This helper was checked against
five hand-computed values before use** and reproduces mean 0.581, SD 0.116,
SEM 0.052 exactly:

```bash
grep -h '"iptm"' boltz_out/*/predictions/*/confidence_*.json | tr -d ' ,"' | cut -d: -f2 > /tmp/iptm.txt
awk '{s+=$1;q+=$1*$1;n++} END{m=s/n;printf "n=%d mean=%.3f sd=%.3f sem=%.3f\n",n,m,sqrt((q-n*m*m)/(n-1)),sqrt((q-n*m*m)/(n-1))/sqrt(n)}' /tmp/iptm.txt
```

`grep '"iptm"'` does not also match `"ligand_iptm"` or `"protein_iptm"`, because
both have an underscore rather than a quote before `iptm`. Verified on a mock
file.

### Practices adopted after being burned

- **Write the predicted number down before running the measurement.** It has
  caught two wrong mental models so far, including both predictions for the
  first Boltz run.
- **Treat any measurement that is identical across independent samples as a bug
  until proven otherwise.** Ten independent designs reporting the same distance
  to a decimal place was the symptom of a phantom atom at the origin, caused by
  an uninitialised awk array index. In awk, an unset subscript is the empty
  string, not 0. Always `BEGIN { na=0; nb=0 }`. Full write-up in
  `CAUTION-origin-atom-bug.md`.
- **Sanity-check a measurement script against a known answer before trusting
  it**, which is why the awk helper above carries its verification.
