# 100 ns MD of p53–MDM2 (and any other two-chain complex)

`73_md_all.sbatch` runs one 10 ns production chunk. These scripts run 100 ns as
ten resumable 10 ns chunks, and read the result per window instead of only
start-vs-end.

## Submit (from `~/rfd/scripts`, after `module load slurm`)

```
sbatch --gres=gpu:1 --time=24:00:00 --exclude=dgxhnode6 \
  --job-name=md100_p53 --output=logs/%x_%j.out --error=logs/%x_%j.err \
  md/75_md_100ns.sbatch mdm2_blind/ref_p53.pdb runs100/p53_mdm2 10
```

Re-submitting that identical line after a crash, a timeout or a node failure
**resumes** from the last finished chunk — build, equilibration and completed
chunks are skipped. A half-written chunk is redone, because a chunk only counts
as finished once both its `.nc` and `.rst` are renamed off `.tmp`.

`--exclude=dgxhnode6` because the prebuilt `pmemd.cuda` may carry no H100
kernels ("no kernel image"); drop it once an H100 run is seen to work.

## Read it

```
srun --cpus-per-task=4 --mem=16G --time=00:40:00 \
  --container-image=/home/hpc_apps/nvidia+nvhpc+25.3-devel-cuda_multi-ubuntu22.04.sqsh \
  --container-mounts=/home/dakshcsed:/workspace \
  bash /workspace/rfd/scripts/md/76_md100_report.sh p53_mdm2
```

Prints, per 10 ns window: binder heavy-atom RMSD against the start, total
interface contacts under 4.5 Å, and contacts for each anchor residue
(default p53 Phe19 / Trp23 / Leu26). The per-anchor rows are the point: on the
PPI platform a mutation that destroyed an interface moved the **total** contact
count by 0.6 against a noise floor of 2.84, while the partner-by-partner
occupancy changed completely. A total-contact screen can miss the event.

## Gates in the build

`70_md_build.sh` aborts before any GPU time if either gate fails:

1. **Chirality** — `tleap` rebuilds hydrogens and side chains, and Amber
   protein force fields carry no improper term on the Cα centre, so a D-residue
   is invisible to energy, RMSD, Rg, DSSP and every contact metric.
2. **Chain separation** — `ref_p53.pdb` and the design complexes carry no `TER`
   between chain A and chain B. If `pdb4amber` fails to insert one, `tleap`
   bonds the last target residue to the first binder residue, giving one
   covalent chain that can never come apart — and the trajectory still reads
   as a perfectly stable complex.

## What 100 ns can and cannot say

It tests whether the complex stays **assembled** under explicit solvent and
thermal motion. It is one-sided: *fell off* is strong evidence, *stayed on* is
weak — barnase–barstar's bound lifetime is ~35 hours against 100 ns simulated.
It does not give affinity, ΔG or k_off, and MM-GBSA cannot rescue it into one
(it returns a precisely wrong **positive** ΔG for trypsin–BPTI, Kd ~10⁻¹⁴ M).
