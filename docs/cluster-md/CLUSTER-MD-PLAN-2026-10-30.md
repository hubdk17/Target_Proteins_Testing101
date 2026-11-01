# Running 100 ns of MD on the Thapar AI Foundry cluster

Written 2026-10-30, from the cluster knowledge earned in the other threads plus the institute's
own Amber documentation page. **Decisions are flagged; nothing here has been submitted.**

---

## 1. How the cluster works (the operating sheet)

These are the facts that have each cost a real job at some point.

**Every session**
- `module load slurm` **first**, on every fresh login to `bcmnode1` (user `dakshcsed`). Without it
  `sbatch`, `squeue` and `sacct` are "command not found" and bash suggests installing a package.

**Resources**
- One partition, `defq`, and it has **no time limit**.
- Three nodes × 8 GPUs = **24 GPUs**: `dgxanode03` and `dgxanode04` are **A100**, `dgxhnode6` is
  **H100**. `sinfo` truncates its node column, so it hides nodes — check nodes directly.
- **Do not `--exclude=dgxanode03`.** Its old enroot disk failure was confirmed cleared 2026-10-28;
  carrying the exclusion silently cut our pool by a third for five days.
- `--gres=gpu:1` is the correct form (GRES is `gpu:8`), **not** `--gpus=1`.
- Pending reason `(Resources)` = no GPU free; `(Priority)` = outranked. A fresh job often shows
  `(Priority)` either way, so check the nodes rather than the reason code.

**Software model — this is the part that differs from every generic HPC guide**
- **No conda. No CUDA modules.** Everything runs inside **enroot** containers from
  `/home/hpc_apps/*.sqsh`. Home is mounted into the container with
  `--container-mounts=/home/$USER:/workspace`.
- `$USER` does **not** expand in `#SBATCH --container-*` lines (they are comments) — container
  flags belong on the `srun` line inside the script.
- The NGC containers ship a `$PIP_CONSTRAINT` file pinning everything with `==`, which is
  unsatisfiable. If you must pip-install: write your own constraint file holding only the
  container's torch/numpy versions, and `export PIP_CONSTRAINT=` to that. Always
  `pip install --dry-run` first as a cheap no-GPU job.
- **Compute nodes are host-filtered** (no `huggingface.co`); the login node has wider access, so
  downloads happen there.

**Habits**
- **`sbatch`, never interactive `srun`**, for anything that matters — an SSH drop has killed a run
  mid-way and left a 0-byte output.
- **Distinct output name per run.** A reused job name or output file poisons log globs, and score
  files in this project are append-only.
- **A job's success is its output files and elapsed time, never its Slurm state.** We have seen
  `FAILED` on a successful run and `COMPLETED` on a total failure.
- **`AllocTRES` is the only honest GPU count:** `sacct -j <id> --format=JobID,AllocTRES%60`. A job
  submitted without `--gres` runs anyway, on 2 CPU cores, with no warning — that is how an AF2 job
  went from 43 s to 11.5 min per complex. Check it in the first minute of any long run.
- Env builds need **no** `--gres`, so they queue immediately. Do CPU work while GPUs are busy.
- Cap environment or queue debugging at about **two rounds**, then step back and say so.

---

## 2. The engine decision: use Amber's `pmemd.cuda`, not OpenMM

**Recommendation: Amber 24.** The institute's docs page documents it and it is already on the
cluster.

| | OpenMM (what the laptop platform uses) | **Amber 24 `pmemd.cuda`** |
|---|---|---|
| On the cluster? | no — would need pip install | **yes**, `/home/hpc_apps/amber/amber24.tar.gz` |
| CUDA? | **unverified here.** Our own README says "pip wheels give CPU and OpenCL only" — but that was measured on the Windows laptop and may not hold on Linux | **native**, that is what `pmemd.cuda` is |
| conda needed? | that is the documented CUDA route, and there is no conda here | no |
| Force field | ff14SB / TIP3P-FB | **the same force field, natively** — this is Amber's own ff |
| Institute support | none | documented sbatch example on the docs site |

The cost of choosing Amber is that the system has to be built as Amber `prmtop` + `inpcrd` with
`tleap` instead of an OpenMM `System`. **The analysis layer survives either way** — MDTraj reads
Amber topologies and trajectories natively, so the existing stage 04–21 scripts do not need
porting. That is what makes this a cheap switch rather than a rewrite.

**What the docs page gives us**
- Install is a copy-and-extract, not a module: `cp /home/hpc_apps/amber/amber24.tar.gz` and
  `cp /home/hpc_apps/amber/amber_full_libs.tar.gz` into your own `~/amber/`, then `tar -xvf` both.
  The support libraries are expected to land under `~/usr`.
- Container: `/home/hpc_apps/nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh`
  (the nvhpc image, **not** the pytorch one RFdiffusion uses).
- Inside the container: `export AMBERHOME=/workspace/amber24`, `source "$AMBERHOME/amber.sh"`, and
  prepend `/workspace/usr/lib/...` (OpenBLAS and x86_64 dirs) to `LD_LIBRARY_PATH`.
- Their example pins `--nodelist=dgxanode03`, requests `--gres=gpu:1`, 4 CPUs, 32 GB.
- The page warns that **`pmemd.cuda` must be compatible with the GPU and CUDA environment**.

**Two things to verify before committing a long job** — both cheap, and both are checks that can fail:
1. **Does the archive include `tleap`?** The page names `cpptraj`, `sander`, `pmemd` and
   `pmemd.cuda` but never mentions `tleap`. Without it we cannot build the topology on the cluster
   and would have to build it elsewhere. One `ls amber24/bin` settles it.
2. **Does `pmemd.cuda` actually run on the GPU here?** A prebuilt binary is compiled for specific
   architectures. **Start on an A100 (`dgxanode03`/`04`), not the H100** — the A100 is the older,
   safer target and is what the institute's own example pins. A one-minute job running
   `pmemd.cuda --version` plus `nvidia-smi` settles it.

---

## 3. What 100 ns will and will not tell you

Worth being straight about, because our own platform is the evidence.

**It will** tell you whether the complex **stays assembled** under a real force field with explicit
water and thermal motion. Nothing currently in the binder pipeline does this: RFdiffusion, MPNN and
AF2 all see a single static structure with no solvent and no temperature. A designed binder that
drifts off its site in 100 ns is dead no matter what its pae_interaction said. **That is a genuine
new filter, and the first time anything physical touches the complex.**

**It will not** give binding affinity, ΔG or k_off, and it is important not to let it look as if it
does:
- Barnase–barstar's bound lifetime is **~35 hours** against 60 ns simulated — about 10⁹×. MD at
  this length does not sample dissociation for a real binder, so "it stayed bound" is weak
  evidence and "it fell off" is strong evidence. **The test is one-sided.**
- Contact occupancy is a geometric quantity, not a free energy.
- MM-GBSA cannot rescue it into an affinity: it returns a **precisely wrong, positive** answer for
  trypsin–BPTI (Kd ~10⁻¹⁴ M) while being just as reproducible as when it is right, and our designed
  binders skew net-negative, which is the regime that breaks it.

Also note this **skips the step the project had queued ahead of MD** — PyRosetta
`InterfaceAnalyzer` plus a mutate/repack ΔΔG, which would be the first side-chain-aware energy in
the pipeline and costs hours rather than GPU-days. Doing MD now is a reasonable call, but it is a
reordering, not a continuation.

---

## 4. The run I would actually propose

**Not 1 × 100 ns. Two arms, same protocol, same batch.**

A single run on the best design is n = 1 and, worse, it is a check that cannot fail: if the complex
stays put we learn nothing, because we have no idea what a *bad* design does under the same
protocol. The house rule here is to run the negative control in the same batch.

| arm | structure | role |
|---|---|---|
| **A** | `design_103_s5` — pae_interaction **6.353**, target_aligned_rmsd 1.213 | does a gate-passing design hold up physically? |
| **B** | a **failing sequence off the same backbone 103** (that backbone's sequences spread 6.35 → 25.90) | does the test discriminate? |

**Arm B should come off backbone 103, not a different backbone.** Same geometry, different
sequence, so the comparison isolates the interface chemistry instead of confounding it with a
different shape — a much tighter-matched negative. And **do not pick arm B by any gate-2 number**:
the backbones that fail gate 3 hardest (0, 106, 111 at 25–27) actually fold *better* at gate 2
than 103 does, so gate 2 would select the wrong control.

(Arm A's numbers are confirmed from the replicate run — score file `~/rfd/cplxRep3.sc`, complex
pushed to the private repo as `design_103_s5.pdb`, interface report 34 binder atoms within 5 Å of
Trp157, nearest 2.54 Å, 83 contacts. Relevant caveat from that thread: pae_interaction is **not**
reproducible across nodes — the same inputs on two different A100s moved docked designs by ~7% —
so a number within ~1 of the gate is not a verdict. 6.353 is far enough under that this does not
threaten it.)

If both hold for 100 ns, MD at this length has no discriminating power for us and we have learned
that cheaply. If A holds and B comes apart, we have a second independent filter — and, unlike
pae_interaction, one that is physical.

**Cost estimate**, to be replaced by a measured rate from the first chunk: ~200 residues solvated
is roughly 35–45k atoms; `pmemd.cuda` on an A100 at 2 fs should give a few hundred ns/day, so
100 ns is plausibly **well under a day per arm**, and the two arms run on separate GPUs in
parallel. `defq` has no time limit, so an overnight submit is fine.

**Write the prediction down before running** (house style). Mine: arm A stays bound for the full
100 ns with the binder shifting but not leaving; arm B separates or slides off the site. If B also
holds, that is the informative failure and it kills the filter.

**Run it in chunks with a checkpoint per chunk.** That pattern is what has made long runs
survivable on this project — restart overhead measured at ~1.7% of wall clock over 30 chunks, and
nothing has ever been lost.

---

## 5. Open questions for you

1. ~~Which structure is arm A?~~ **Settled** — `design_103_s5`, pae_interaction 6.353, the
   project's first design comfortably under the calibrated gate. Reported from the IL-6 thread
   rather than checked by me against the score file.
2. **One arm or two?** I recommend two, for the reason in §4. **This is the open decision.**
3. **Getting the structures to the cluster.** `design_103_s5.pdb` is already in the private repo,
   so a `git pull` reaches it. The arm B sequence needs its complex built and pushed the same way —
   transfers here go through the repo, not scp.

---

## 6. First step, when you want it

Install Amber and answer the two verification questions in §2 — one copy-and-extract plus a
one-minute GPU job. No long run until `tleap` is confirmed present and `pmemd.cuda` has printed a
version on an A100.
