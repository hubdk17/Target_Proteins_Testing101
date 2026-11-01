# The per-target procedure — moved

The method now lives in one place, written up from the IL-6 run:

**[../pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md](../pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md)**

It holds the five stages, the three gates, the thresholds and the commands. Read
that, not this file. An earlier version of this page restated the same method and
was deleted rather than kept, because two copies of a procedure drift apart and
this folder is written by several threads.

## What lives here instead

The target-specific machinery for running that procedure across many targets:

| File | Stage | What it does |
|---|---|---|
| `02_prep_targets.py` | 0 | Measures each PDB: chain, real breaks, contig, hotspots from the bound partner |
| `31_rfd_multi.sbatch` | 1 | Backbones, one array task per target, ROG potential on |
| `32_gate1_compactness.sh` | Gate 1 | Rho per design, per target; writes the survivor list |
| `33_mpnn_multi.sbatch` | 2 | MPNN sequences for gate-1 survivors only, 16 per backbone, per target |
| `thread_multi.py` | Gate 2 | Threads MPNN sequences onto their own backbones, per target |
| `40_gate2_selfconsistency.sbatch` | Gate 2 | AF2 `-force_monomer`, fresh score file per target, best RMSD per backbone |
| `41_gate2_report.sh` | Gate 2 | Reads the scores, applies the 2 Å cut, flags undesignable backbones |
| `42_grip_rank.py` | ranking | Grip per residue and hotspot coverage, among gate-2 survivors only |
| `43_build_complex.py` | 3 | Stage-3 inputs: binder chain first, numbering unique across chains, verified |
| `BATCH-PLAN.md` | — | Which targets are in batch 1 and why the other six wait |
| `TARGET-TRIAGE-2026-10-29.md` | — | All 19 targets sorted into tiers with PDB entries |

The split is deliberate: the procedure document says *what* the gates are, these
scripts handle *many targets at once* without name collisions, which is the part
that bites at twelve targets and not at one.
