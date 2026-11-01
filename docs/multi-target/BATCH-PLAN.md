# Batch 1 — widening the pipeline past IL-6

Decision by Kaila7, 2026-10-29: put the other targets on the cluster too, since the
jobs will finish in time anyway. Triage and the reasons behind the tiers are in
`TARGET-TRIAGE-2026-10-29.md`.

PD-1 goes to the front of the queue regardless of the batch width, because it is the
only target on the list with published lab-confirmed de novo binders and so is the
only one that can tell us whether our AF2 gate is calibrated. It costs nothing to
order it first.

## The 12 targets being prepared

Everything the triage called workable, minus INSR (large, glycosylated, disulfide-
linked dimer — ruled out 2026-10-21 and nothing has changed).

| # | Target | PDB | Which part we design against | Note |
|---|---|---|---|---|
| 1 | PD-1 | 3RRQ | whole IgV domain | diagnostic control, runs first |
| 2 | PCSK9 | 2P4E | catalytic domain | EGF-A groove is the epitope (see 3BPS) |
| 3 | BCL2 | 2XA0 | BH3 groove | easiest surface on the list; intracellular |
| 4 | STAT3 | 1BG1 | SH2 domain only | file is a dimer on DNA — one chain, protein only |
| 5 | Casp3 | 2J30 | p17 chain | heterotetramer; pick one chain pair |
| 6 | pTEN | 1D5R | phosphatase + C2 | flexible C-tail to trim |
| 7 | ERK2 | 4QTB | kinase domain | **provisional** — "MAPK" not yet pinned down |
| 8 | JAK1 | 3EYG | JH1 kinase domain | |
| 9 | JAK2 | 4IVA | JH1 kinase domain | |
| 10 | Beta-RAF | 4MNE | kinase domain | dimerises; use one chain |
| 11 | BRCA1 | 1T29 | BRCT tandem domain only | rest of BRCA1 is disordered |
| 12 | GLP1R | 3IOL | extracellular domain only | whole receptor is 7-TM, out of scope |

## Batch 1 is 6 of the 12, split on one criterion: is the epitope measurable?

Hotspots decide whether a run is aimed at anything real, and reciting epitope
residue numbers from memory is exactly the kind of unverified claim that has
already cost this project a week. So batch 1 is the targets whose PDB entry
**already contains the binding partner** sitting in the site we want to block.
For those, the epitope is measured off the crystal rather than asserted.

| # | Target | PDB | Partner in the file | Site |
|---|---|---|---|---|
| 1 | PD-1 | **4ZQK** | PD-L1 | PD-L1 interface |
| 2 | PCSK9 | **3BPS** | EGF-A domain of LDLR | EGF-A groove |
| 3 | BCL2 | 2XA0 | inhibitor | BH3 groove |
| 4 | BRCA1 BRCT | 1T29 | phosphopeptide | phosphopeptide site |
| 5 | GLP1R ECD | 3IOL | exendin-4 | peptide-binding face |
| 6 | Casp3 | 2J30 | active-site inhibitor | substrate groove |

PD-1 moves from `3RRQ` to `4ZQK` for this reason: the apo structure is a cleaner
file, but it cannot tell us where PD-L1 binds.

**Batch 2 waits on a decision, not on compute.** pTEN (1D5R), ERK2 (4QTB), JAK1
(3EYG), JAK2 (4IVA) and Beta-RAF (4MNE) have no relevant partner in the file, and
STAT3 (1BG1) is bound to DNA, which is not a protein epitope. For all six, the site
has to be chosen from the literature. They are all intracellular anyway, so they are
the lower-value half of the list.

## Why this cannot go straight to sbatch

A contig string encodes real residue numbers and real chain breaks. On IL-6 we
verified the numbering twice and still found a genuine chain break at 51 -> 61 that
would have silently corrupted the contig. We have no numbering for any of these 12
yet, and guessing it is exactly the class of error that has already cost this
project a week once.

So the order is:

1. **Fetch and inspect** (login node, no GPU, about a minute). Download each file and
   print its TITLE, chain letters and CA count. **The TITLE line is the check** — it
   is how we find out if any ID is the wrong protein, the wrong species, or a
   complex we did not expect. Nothing gets queued until those 12 titles have been
   read.
2. **Per target: strip to protein, pick the chain, find chain breaks, write the
   contig.** This is mechanical but it has to be done per file.
3. **Pick hotspots per target from its mapped epitope.** This is the only step that
   is science rather than plumbing, and it is what determines whether a run is
   aimed at anything real. Hotspot conditioning also roughly triples the
   compactness pass rate, so it is not optional.
4. **One Slurm array job across targets**, carrying the ROG potential, which is
   load-bearing: without it two thirds of designs come out as straight rods.
5. **Self-consistency gate on every design, before any interface scoring.** See
   below — this is new, and for a 12-target batch it is the most important step
   on this list.

## The self-consistency gate (added 2026-10-29, corrected the same day)

**Correction first.** An earlier version of this file said ProteinMPNN had walked the
straight-rod failure back into the pipeline, based on one sequence folding into a
rod. **That was wrong, and generalised from a single data point.** Job 174717 ran 5
backbones x 8 MPNN sequences = 40 structures through AF2 with `-force_monomer`, and
**31 of 40 came in under 2 Å**. MPNN is not broken.

What is true is narrower and still serious: **some sequences fail, and nothing in the
pipeline was checking which.** The spread per backbone:

- `short_2`: all eight sequences at **0.34-0.40 Å**, pLDDT 92-95. Excellent.
- `short_9`: six sequences at 1.07-2.21 Å, but **samples 2 and 4 fold into straight
  rods at 14.1 and 14.4 Å**.

**Every downstream result in this project so far used `short_9` sample 2 — one of
the two bad ones — and nothing caught it.** That is the argument for the gate: not
that MPNN is unreliable, but that it has a failure rate and we were sampling blind.

### pLDDT is inverted here — never rank by confidence

Inside `short_9`, the **two rods score pLDDT 93.4 and 94.8** while the **six correct
folds score 75-79**. A straight helix is an easy thing to predict; a hairpin is hard.
So confidence is anti-correlated with correctness in exactly the case that matters.

**Only RMSD against the designed backbone can be trusted to pick a sequence.** Any
ranking by pLDDT would have actively selected the rods.

### The gate, and where it is written up

The method, thresholds and commands live in one place:
**`docs/pipeline-procedure/BINDER-PIPELINE-PROCEDURE.md`**. In short:
make **several MPNN sequences per backbone** (the failure rate is per sequence, so
one sample is a coin flip), predict each **alone** with `-force_monomer`, superpose
on its own backbone, keep under **~2 Å**, and discard the rest whatever their
confidence or interface score.

It comes **before** interface scoring because it is about a fifth the cost —
marginal at one target, decisive at twelve.

    backbones
      -> gate: compactness (Rg ratio under ~1.3)       32_gate1_compactness.sh
      -> several MPNN sequences per surviving backbone
      -> gate: self-consistency RMSD under ~2 A        thread_multi.py
                                                      40_gate2_selfconsistency.sbatch
                                                      41_gate2_report.sh
      -> rank: grip (interface contacts)
      -> score: interface (AF2 complex, pae_interaction < ~10)

**Fold is a gate; grip is a ranking applied only among designs that already passed
it.** Ranking by grip before gating on fold is how this project selected unfoldable
binders once already.

## Expected cost

At IL-6's rate (~1 min per design on an A100) and 50 designs per target, stage 1 for
all 12 is roughly 10 GPU-hours. That is the cheap stage. Stage 3 scoring is what
actually costs, which is the argument for reading PD-1's control number before
scoring all 12.
