# `TAG=m1` three-arm comparison - results 2026-11-01

Jobs 175791 and 176028, 150 items per cell. Binder length range 40-55, T=50,
ROG potential on. Inputs verified from `.hydra/overrides.yaml` and the in-script
count guard (`hotspots received: ... (5 of 5 expected)`).

Built by `50_build_mutation_arms.py`, read out by `52_mutation_report.py`.

## The three input variants

| arm | name field | side-chain atoms kept | varies |
|---|---|---|---|
| `wt` | original | all | baseline |
| `mut` | `ALA` | N, CA, C, O, CB | name field and atom count |
| `ctl` | original | N, CA, C, O, CB | atom count only |

`mut` and `ctl` carry identical atoms and differ only in columns 18-20, so that
contrast isolates the name field. The builder refuses unless that holds, atom
counts match, and every file has one unique ATOM line length.

**Statistic:** chain-B CA atoms within 8 A of the selected residue's CA, per
output file. Measured to the CA because the arms differ in side-chain atom count
by construction, so a nearest-any-atom statistic would measure the edit instead.
The name field is checked in every file before it is counted.

## 3BPS, residue A379 (`PHE`) - complete

| arm | n | mean 8A | sd | sem | mean 10A | nearest |
|---|---|---|---|---|---|---|
| `wt` | 150 | 2.49 | 1.24 | 0.101 | 6.67 | 6.56 |
| `mut` | 150 | 2.93 | 1.21 | 0.099 | 7.09 | 6.27 |
| `ctl` | 150 | 2.35 | 1.24 | 0.102 | 6.41 | 6.57 |

Welch t: **`wt` vs `mut` -3.16**, `ctl` vs `wt` -0.93, **`ctl` vs `mut` -4.10**.
Control lands -30% of the span from `wt` to `mut`.

`ctl` and `wt` are indistinguishable, so removing the atoms while keeping the
name field changes nothing measurable; the separation is carried by the
name-field contrast. Same direction and structure as the earlier 1ALU result, on
an independently built input.

## 2J30, residue A206 (`TRP`) - complete

| arm | n | mean 8A | sd | sem | mean 10A | nearest |
|---|---|---|---|---|---|---|
| `wt` | 150 | 2.06 | 1.54 | 0.126 | 4.79 | 6.70 |
| `mut` | 150 | 2.33 | 1.63 | 0.133 | 4.97 | 6.60 |
| `ctl` | 150 | 1.83 | 1.53 | 0.125 | 4.69 | 6.96 |

Welch t: `wt` vs `mut` **-1.49**, `ctl` vs `wt` -1.32, `ctl` vs `mut` **-2.78**.

**The primary contrast does not clear |t| = 2**, so the readout refuses to
attribute the span to the control; dividing by a span that is itself noise
manufactures a number. Recorded rather than discarded: the ordering matches 3BPS
(`mut` > `wt` > `ctl`) and the name-field contrast is the largest of the three at
-2.78, so plausibly the same effect at a smaller size. That needs more items per
cell, not a different reading of these.

## File counts per output directory

| directory | files of 150 |
|---|---|
| `bcl2_wt` / `bcl2_mut` / `bcl2_ctl` | 150 / 150 / 150 |
| `casp3_wt` / `casp3_mut` / `casp3_ctl` | 150 / 150 / 150 |
| `pcsk9_wt` / `pcsk9_mut` / `pcsk9_ctl` | 150 / 150 / 150 |
| `pd1_wt` / `pd1_mut` / `pd1_ctl` | 150 / 150 / 150 |
| `glp1r_wt` / `glp1r_mut` / `glp1r_ctl` | 28 / 150 / 150 |
| `brca1_wt` / `brca1_mut` / `brca1_ctl` | 95 / 86 / 68 |

Four of six complete. `glp1r_wt` (task 13) and the three `brca1` tasks were still
running at the time of the listing. **Unequal n between arms is the one way this
comparison misleads, so read a target only once all three directories hold 150.**

`pd1_*` is 4ZQK chain A, which `COMPND` identifies as the ligand, not the
receptor. The files keep the old name; rename at a clean boundary, never mid-run.

## Scheduling notes

- `--exclude=dgxanode03` throughout: its enroot cache filled and tasks died at
  container start in under two minutes. Across two jobs, 23 failures to 2
  successes on that node against 0 failures elsewhere. Intermittent, so it may
  recover; take the exclusion off when it does.
- Rate at 3BPS (435 residues), measured from files on disk on a *running* task:
  **~165 s per item**. Not the ~379 s obtained by dividing a cap by a killed
  job's output, which charges container setup and model load to the items and
  freezes one hour's contention as if it were constant. **Count files on the
  running job.**
- Judge completeness by file count per directory, never by `sacct` State. This
  project has seen a job report FAILED on success and COMPLETED on total failure.

## Open

- A fourth variant (name field `ALA`, side chain retained) would complete the
  2x2 and test whether the name-field effect survives when the atoms are
  present. Probably not constructible from a file edit: if `parse_pdb_lines`
  keys atoms to the name field, the extra atoms are dropped at load and the
  variant silently becomes `mut`. Entry points are
  `RFdiffusion/rfdiffusion/inference/utils.py:519` and `:526` - read the atom
  handling there before building anything.
- A size ladder at one position, every arm a real residue so all parse natively,
  is the better alternative and needs no override. Its cost is genuine
  side-chain rebuilding, so every input must clear `44_chirality_check.py`
  first. The 1ALU thread owns the builder.
- Stage 2 on the four complete targets, 94 gate-1 survivors from the separate
  `h5` pool.
