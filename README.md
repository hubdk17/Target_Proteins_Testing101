# Multi-target binder design scripts

Scripts for running the binder design pipeline across many protein targets on the
Thapar AI Foundry cluster. They expect to live in `~/rfd/` next to the existing
`RFdiffusion`, `ProteinMPNN` and `dl_binder_design` clones and the `venv*`
environments.

## Setup on the cluster

```
cd ~/rfd
git clone https://github.com/hubdk17/Target_Proteins_Testing101.git scripts
cp scripts/* .
```

Afterwards, `cd ~/rfd/scripts && git pull && cp * ..` picks up any changes.

## The order, and where the gates sit

```
Stage 0  prepare targets        02_prep_targets.py          login node, CPU
Stage 1  backbones              31_rfd_multi.sbatch         GPU, array over targets
GATE 1   compactness            32_gate1_compactness.sh     login node, CPU
Stage 2  sequences              33_mpnn_multi.sbatch        GPU, array over targets
GATE 2   self-consistency       thread_multi.py             login node, CPU
                                40_gate2_selfconsistency.sbatch   GPU
                                41_gate2_report.sh          login node, CPU
         rank by grip           42_grip_rank.py             login node, CPU
Stage 3  complex inputs         43_build_complex.py         login node, CPU
```

**Gates throw designs out; rankings only sort what survived.** A design must pass
every gate before any ranking is applied. Ranking before gating once carried a
broken design through this whole project.

| Gate | Measure | Keep |
|---|---|---|
| 1 | Rg / (2.2 · N^0.38) | below ~1.30 |
| 2 | Cα RMSD of the sequence folded alone vs its designed backbone | below ~2 Å |
| 3 | AF2 `pae_interaction` | below ~10 |

**Never rank sequences by pLDDT.** In the IL-6 run the two sequences that folded
into straight rods had the *highest* confidence (93–95) while the six correct folds
sat at 75–79: a simple wrong fold is easy to predict. Only RMSD against the
designed backbone can pick a sequence.

## Running order, one target at a time

```
cd ~/rfd
python3 02_prep_targets.py                              # read the output before submitting
sbatch --array=1 31_rfd_multi.sbatch                    # or --array=1-6 for the batch
bash 32_gate1_compactness.sh
sbatch --array=1 33_mpnn_multi.sbatch
python3 -m py_compile thread_multi.py && python3 thread_multi.py
sbatch --array=1 40_gate2_selfconsistency.sbatch
bash 41_gate2_report.sh
python3 42_grip_rank.py
python3 43_build_complex.py
```

`module load slurm` once per fresh login, or `sbatch` reports "command not found".

## What each script guards against

Every one of these cost a real job or a round of wrong reasoning before it was
written down.

- **`02_prep_targets.py`** — tells a *real* chain break (CA–CA above ~4.5 Å) from a
  numbering gap, which look identical until measured. Derives hotspots from the
  partner in the crystal instead of from recited residue numbers, and warns when
  they do not form a single patch.
- **`31_rfd_multi.sbatch`** — carries the radius-of-gyration potential, which is
  load-bearing: without it 176 of 300 designs came out as single straight helices.
- **`32_gate1_compactness.sh`** — flags a rod-like *median*, which points at the run
  rather than at luck. Computes its median without gawk's `asort`, which is not
  portable.
- **`33_mpnn_multi.sbatch`** — feeds only gate-1 survivors, writes 16 sequences per
  backbone because the gate-2 failure rate is per sequence, and asks for a GPU
  (CPU is ~11 min per structure against seconds).
- **`thread_multi.py`** — asserts every residue lookup and every line length.
  Fixed-column PDB edits have silently corrupted data twice here; it refuses a
  sequence whose length does not match its backbone and checks for blank residue
  names.
- **`40_gate2_selfconsistency.sbatch`** — fresh timestamped score file every run
  (the score file is append-only), no `-debug` so one bad input cannot kill the
  batch, and judges success by output count rather than exit code.
- **`43_build_complex.py`** — enforces the two `predict.py` requirements that fail
  *silently*: the binder must be the first chain, and residue numbers must be
  unique across the whole file regardless of chain letter.

## Operational rules

- A job's success is its **output files and elapsed time**, never its Slurm state.
  One job here reported FAILED on success; another COMPLETED having produced
  nothing.
- **Distinct output names per run.** A reused name has overwritten a batch, and it
  poisons both the log glob and the `--name` lookup.
- **Never interactive `srun`** for a pipeline step; an SSH drop kills the job.
- Write the predicted number down **before** measuring it.
- Treat any invariant across independent samples as a bug until proven otherwise.
