# Documentation index

Write-ups produced alongside the scripts in this repository. Collected here
2026-11-01 from the project's shared folder, which several working threads write
into; this directory is the version-controlled copy.

Each entry carries the date the document itself was last written. A document is
**current** when nothing measured since its date contradicts it, and
**superseded** when a later document on the same subject exists. Section 6 lists
results that exist as job output on the cluster but have no write-up in this
repository yet.

## 1. Start here

| File | Date | What it is |
|---|---|---|
| [pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md](pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md) | 2026-10-29 | The method: five stages, three gates, thresholds, commands. The one procedure document; everything else defers to it. **Current.** |
| [report/plain-language-project-report.html](report/plain-language-project-report.html) | 2026-11-01 | Whole-project report in plain language, for a reader who has not seen the code. **Current.** |
| [reports/what-the-pipeline-misses.pdf](reports/what-the-pipeline-misses.pdf) | 2026-11-01 | Same evidence, external-reader version. **Current.** |

The repository's top-level [README.md](../README.md) covers the scripts and the
order they run in.

## 2. IL-6 (1ALU) — the lead case

`rfdiffusion-il6/`. The first target taken through all three stages, and where
most of the pipeline's behaviour was characterised.

| File | Date | What it is |
|---|---|---|
| [rfdiffusion-il6/README.md](rfdiffusion-il6/README.md) | 2026-10-28 | Index and state of play for the IL-6 run. Predates section 6's results. |
| [rfdiffusion-il6/METHODS.md](rfdiffusion-il6/METHODS.md) | 2026-10-23 | How to run and measure each stage on this target. |
| [rfdiffusion-il6/RESULTS.md](rfdiffusion-il6/RESULTS.md) | 2026-10-28 | Chronological run log, deliberately keeping conclusions that were later overturned, with the corrections marked in place. Read as history, not as summary. |
| [rfdiffusion-il6/REPORT-2026-10-28.md](rfdiffusion-il6/REPORT-2026-10-28.md) | 2026-10-28 | The long-form report as of that date. |
| [rfdiffusion-il6/HOTSPOTS-COMPACTNESS.md](rfdiffusion-il6/HOTSPOTS-COMPACTNESS.md) | 2026-10-28 | Gate 1: the Rg-ratio compactness gate and the straight-helix failure mode it catches. |
| [rfdiffusion-il6/STAGE2-mpnn-audit.md](rfdiffusion-il6/STAGE2-mpnn-audit.md) | 2026-10-28 | What ProteinMPNN does and does not read from the target. |
| [rfdiffusion-il6/STAGE3-af2-initial-guess.md](rfdiffusion-il6/STAGE3-af2-initial-guess.md) | 2026-10-28 | The scoring stage in use: `pae_interaction`, the gate threshold, and the noise measurement. **Current for stage 3 mechanics.** |
| [rfdiffusion-il6/STAGE3-boltz-scoring.md](rfdiffusion-il6/STAGE3-boltz-scoring.md) | 2026-10-23 | The Boltz ipTM scoring attempt and why it was dropped. **Superseded as the scoring route** by the file above; kept because the negative result stands. |
| [rfdiffusion-il6/LADDER-RESULT.md](rfdiffusion-il6/LADDER-RESULT.md) | 2026-10-28 | The stage-3 residue ladder with its negative control. |
| [rfdiffusion-il6/CAUTION-origin-atom-bug.md](rfdiffusion-il6/CAUTION-origin-atom-bug.md) | 2026-10-28 | A check that could not fail, and how it was found. Read before writing a new check. |
| [rfdiffusion-il6/BATCH-2026-10-25.md](rfdiffusion-il6/BATCH-2026-10-25.md) | 2026-10-25 | The unattended batch left running over that weekend. |

Run scripts for this target are beside the documents
(`01_fetch_il6.sh`, `30_rfd_il6.sbatch`, `run_mpnn.sh`), with the ChimeraX view
script in `rfdiffusion-il6/chimerax/`.

## 3. Many targets

`multi-target/`, alongside the numbered scripts at the repository root.

| File | Date | What it is |
|---|---|---|
| [RESULTS-m1.md](RESULTS-m1.md) | 2026-11-01 | `TAG=m1` three-arm comparison, jobs 175791 and 176028: per-arm tables, t-values, per-directory file counts, the `dgxanode03` exclusion and the corrected per-item rate. **Current.** Added by the thread that ran it; sits at the top of `docs/` rather than in `multi-target/`. |
| [multi-target/TARGET-TRIAGE-2026-10-29.md](multi-target/TARGET-TRIAGE-2026-10-29.md) | 2026-10-29 | The 19 candidate targets, tiered, with the reason for each. |
| [multi-target/BATCH-PLAN.md](multi-target/BATCH-PLAN.md) | 2026-10-29 | The batch that widened the run past IL-6. |
| [multi-target/PARSER-ATOM-NAME-LIMITATION.md](multi-target/PARSER-ATOM-NAME-LIMITATION.md) | 2026-11-01 | `inference/utils.py: parse_pdb_lines` matches atom records by name against the canonical tuple for the label on the record. No positional fallback, no warning. **Current, and it settles which input variants are constructible at all** — read it before designing any file-edit experiment. |
| [multi-target/PROCEDURE.md](multi-target/PROCEDURE.md) | 2026-10-29 | A pointer only; the method moved to section 1. |

## 4. Other workstreams

| File | Date | What it is |
|---|---|---|
| [bindcraft/BINDCRAFT-COMPARISON-PLAN.md](bindcraft/BINDCRAFT-COMPARISON-PLAN.md) | 2026-11-01 | BindCraft read as a pipeline, its config defaults that delete the outputs a fold-back check needs, and the one comparison worth running. Exploration and benchmark, not adoption. **Current; nothing submitted yet.** |
| [cluster-md/CLUSTER-MD-PLAN-2026-10-30.md](cluster-md/CLUSTER-MD-PLAN-2026-10-30.md) | 2026-10-30 | Plan for running 100 ns of MD here, with the cluster operating sheet. **Current; nothing submitted yet.** |
| [ppi-interaction/ANALYSIS-2026-10-30.md](ppi-interaction/ANALYSIS-2026-10-30.md) | 2026-10-30 | Audit of the separate barnase–barstar MD platform: what is in it, what was verified, what it changes. **Current** as of that folder's own last write, 2026-10-20. |

## 5. Cluster environment

`rfdiffusion-trypsin/` holds the first target run on this cluster, and is kept
because the environment build is documented nowhere else.

| File | Date | What it is |
|---|---|---|
| [rfdiffusion-trypsin/README_cluster.md](rfdiffusion-trypsin/README_cluster.md) | 2026-10-16 | Cluster setup, the container and venv pattern, the DGL details. Still the reference for a fresh environment. |

Its scripts (`01_fetch.sh`, `10_hello_gpu.sbatch`, `20_buildenv.sbatch`,
`30_rfd.sbatch`, `chaincheck.sh`) are beside it.

## 6. Results on the cluster with no write-up here yet

These exist as job output and score files on the cluster. The numbers are not
reachable from this repository, so nothing is stated for them here; each needs
its output read out before a document can be written.

- **Stage-1 label experiment on 1ALU** — the three-arm run at the IL-6 position,
  the counterpart to [RESULTS-m1.md](RESULTS-m1.md). Builder
  `50_build_mutation_arms.py`, readout `52_mutation_report.py`.
- **`glp1r_*` and `brca1_*` arms of `TAG=m1`** — incomplete at the time
  `RESULTS-m1.md` was written (`glp1r_wt` 28 of 150; the three `brca1`
  directories 95 / 86 / 68). Unequal n between arms is the one way that
  comparison misleads, so those two targets are read only once all three
  directories hold 150 files.
- **Confidence against fold-back** — `53_confidence_vs_foldback.py`,
  `plddt_binder` against `binder_aligned_rmsd` over the stage-2 pool. The script
  itself is not in this repository.
- **PyRosetta interface energy** — job 176357, `90_interface_energy.py`,
  `dG_separated` and `dSASA` columns with a replicate-derived noise floor. The
  script is not in this repository either.
- **The stage-3 design pool from 2026-10-30** — the `pae_interaction` sweep whose
  survivors the later arms were built from.

Adding any of these means pushing the script with the write-up, so the numbers
and the code that produced them land together.
