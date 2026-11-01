# IL-6 binder design — project index

De novo protein binder design against **interleukin-6**, PDB **1ALU**, on the
Thapar AI Foundry GPU cluster. The point of the project is not to reproduce
AlphaFold; it is to design binders against a specific target and, along the
way, to audit where the standard three-stage pipeline is and is not sensitive
to the chemistry of the target surface.

Last updated 2026-10-28.

## The files

| file | what it is |
|---|---|
| `README.md` | this index, and the state of play below |
| `METHODS.md` | how to run each stage, and the measurement rules, with exact commands |
| **`LADDER-RESULT.md`** | **the project's main result so far: the co-folding gate fails** |
| `BATCH-2026-10-25.md` | design key for the eight-condition ladder: every input, every pre-registered reading |
| `STAGE3-boltz-scoring.md` | the first scoring run and the noise floor, self-contained |
| `STAGE2-mpnn-audit.md` | design for the ProteinMPNN sensitivity test, with its prediction written down |
| `STAGE3-af2-initial-guess.md` | the AF2 environment, the input decisions, and the ladder transferred to the calibrated filter |
| **`HOTSPOTS-COMPACTNESS.md`** | **hotspots triple the compactness pass rate at n=150 per arm** |
| `RESULTS.md` | the chronological run log, 600+ lines, **contains superseded sections** |
| `CAUTION-origin-atom-bug.md` | a measurement bug that invalidated three rounds of reasoning, kept deliberately |
| `01_fetch_il6.sh` | fetch and clean 1ALU |
| `30_rfd_il6.sbatch` | stage 1 job script, parameterised by environment variables |
| `run_mpnn.sh` | stage 2 runner |
| `chimerax/view_short9.cxc` | ChimeraX scene for the designed complex |

`RESULTS.md` is a log, not a summary. It records what was believed at each
point, including conclusions later shown to be wrong, with the corrections
marked in place. Read it for the reasoning; read this file and
`LADDER-RESULT.md` for what is currently true.

## The pipeline

1. **RFdiffusion** generates binder backbones conditioned on the target's
   backbone and a list of hotspot residues. Output is backbone only, no sequence.
2. **ProteinMPNN** writes sequences onto those backbones.
3. **A co-folding model** re-predicts the complex, and interface confidence
   metrics are used to decide which designs are worth making.

Stage 3 is the weakest link and the highest-value place to contribute. As of
2026-10-28 it is also the stage shown to be broken as configured.

## The target

IL-6, 1ALU chain A, 1.90 Å. 157 CA residues, numbered **19–184**.

- **One real chain break, 51 → 61**, a disordered loop. Residues 1–18 are a
  missing N-terminal tail. Hence the contig is two segments: `A19-51/A61-184`.
- Four helices: **21–47 (A), 80–104 (B), 107–129 (C), 156–182 (D)**.
- Three functional surfaces: site I binds IL-6R alpha, site II the first gp130,
  site III the second gp130's D1 domain. **Trp157 is the dominant hotspot** on
  site III, which is the surface being targeted.
- Numbering verified against the coordinates: 19 LEU, 47 SER, 157 TRP, 184 MET.
  W157 is index 139 of the chain A sequence string used as Boltz input.

## State of play

### Established

- **The co-folding gate fails.** At n = 50 per condition, a random 40-mer scores
  ipTM 0.502 against IL-6 while the designed binder scores 0.466, and a
  deliberately alanine-scanned dead interface scores 0.393, *below* the random
  peptide. Single-sequence Boltz with an unpinned target cannot rank binder
  designs. Full tables in `LADDER-RESULT.md`. **This is the main result.**
- **The controls prevented a false positive.** W157A dropped ipTM by 0.072 at
  t = −2.76 and would have read as a result on its own. But E93A, chosen because
  it should be irrelevant, dropped *more* (0.090). The model responds to any
  perturbation of the target sequence, not to interface chemistry.
- **The mechanism is target foldability.** Interface confidence tracks
  whole-complex pLDDT at r = +0.78. IL-6 alone scores pLDDT 0.667 ± 0.090 in
  single-sequence mode, on a 1.90 Å crystal structure we handed the model.
- **`short_9` folds on its own.** `binderonly` pLDDT **0.905 ± 0.028**, the
  highest and tightest number in the table. The stage 2 worry that its
  conserved hydrophobics faced only the interface is resolved in the design's
  favour.
- **Compactness must be a gate, not a tiebreaker.** Ranking backbones by contact
  count actively selected for unfoldable designs, because an extended ribbon
  draped over the target racks up contacts by being spread out. Gate on radius
  of gyration against `2.2 · N^0.38` first, keep ratios under ~1.3, and only then
  look at contacts. ProteinMPNN returned poly-alanine on the top
  contact-ranked design and real amphipathic sequences on the compact ones.
- **`short_9` is a textbook helical hairpin.** 40 residues, Rg ratio 1.21, best
  MPNN score 1.3046. Hydrophobics recur at i, i+4, i+8 on two separate helices
  across all eight samples.
- **A specific RFdiffusion failure mode: single straight helices.** Without the
  ROG guiding potential, 44 of 60 length-50 designs came back at Rg ratio
  2.222 ± 0.023 — the arithmetic value for 50 collinear points at 1.5 Å rise.
  80% attrition against the compactness gate, versus 45% predicted.
- **Hotspot conditioning works.** Specifying `A47,A155,A157,A158` moved site III
  contacts from 2 of 10 to 7 of 7. **n = 1 per arm**, so this is the most
  quotable number in the project and also the weakest.

### Not established

- **No design is known to bind.** Nothing has been validated.
- **Whether a *pinned* target rescues the gate.** The template ladder is built
  (`~/rfd/boltz_batch_t/`) and not yet run. Prediction written down: if pinning
  worked, `targetonly` pLDDT must rise above 0.667 and its spread must shrink.
- **Whether AF2 initial-guess behaves differently.** Untested; this is the
  calibrated filter and the one the conclusion should rest on.
- **Whether ProteinMPNN is sensitive to a target point mutation.** The one
  unaudited stage, and about a minute of GPU time to test.

## The research question

For five weeks this project rested on the claim that **RFdiffusion's binder mode
cannot see a target surface point mutation**, because it conditions on backbone
coordinates and hotspot residue indices rather than side-chain identity, making
W157 and W157A the same input. It was graded the project's one contribution with
legs.

> **That claim is false, established 2026-10-28 by reading
> `rfdiffusion/inference/model_runners.py`.**
>
> ```
> 454:  seq_t = torch.full((1, L_mapped), 21)                        # all mask token
> 455:  seq_t[contig_map.hal_idx0] = seq_orig[contig_map.ref_idx0]   # REAL sequence copied in
> 461:  seq_t[~self.mask_seq.squeeze()] = 21                         # re-masked only where not visible
> 462:  seq_t = one_hot(seq_t, num_classes=22)
> ```
>
> `seq_orig` is `self.target_feats["seq"]`, read from the target PDB. Line 455
> copies it into every mapped position and line 461 re-masks only the diffused
> ones; line 478 names the same tensor `visible=`, which settles the ambiguous
> name `inpaint_seq`. **RFdiffusion is handed IL-6's actual sequence.**
>
> The other half of the claim stands: lines 683 to 704 build `hotspot_tens` as
> zeros with ones at the hotspot positions, appended to `t1d`, so the hotspot
> signal really is one bit per residue.

### What the question becomes

Better, not worse. It stops being architectural and becomes **empirical**: the
model holds the information, does it act on it? That has a real answer either
way, where "it cannot see it" was a statement about an input format.

**And the experiment is now worth running, where before it was pointless by
definition.** Fifty designs against wild-type IL-6 and fifty against a W157A
target PDB, identical settings, compared as distributions. If they are
indistinguishable, that measures a model ignoring information it was given,
which is a far stronger finding.

### The improvement, restated correctly

Not "add the sequence". Side-chain identity arrives as a **bare 22-way one-hot
in a template feature channel**, carrying no side-chain geometry, no chemistry
and no notion of complementarity, while the signal saying "this residue matters"
is **a single bit**. The model is told that residue 157 is a tryptophan and that
it matters, and nothing about where its rings point, how much volume they
occupy, or what would pack against them.

That is the original wide-representation hypothesis in defensible form: the
useful capacity is not MSA width, it is the target-conditioned interface
channel, which today is one one-hot plus one bit.

### The method lesson

A claim graded as the project's main contribution went five weeks without anyone
reading the twelve lines of code that decide it. **Check a load-bearing premise
in the source before building on it**, especially one that is cheap to check and
comfortable to believe.

So the question is where in the pipeline sensitivity enters at all:

| stage | sensitive to W157A? |
|---|---|
| RFdiffusion backbone generation | **Unknown and now testable.** It receives the target sequence as a one-hot, so it *could* respond; whether it does has never been measured |
| ProteinMPNN | sees the target's sequence, so in principle yes — **untested** |
| co-folding re-prediction | **No, not usefully.** It moves, but it moves just as much for an irrelevant mutation, because it is reading target foldability rather than interface chemistry |
| alchemical free energy | yes, by construction (parked, see below) |

The emerging answer is that **no stage of the pipeline ever asks the binding
question directly.** Stage 1 discards target side-chain identity by
construction. Stage 2 sees it but is already committed to a backbone chosen
without it. Stage 3 responds to how easy the target is to fold. That is an
architecture property, not a bug in any one tool.

Two limits on that claim, to state plainly: AlphaFold is untested, and the step
from "the metric does not discriminate here" to "the pipeline cannot see this"
is interpretation, not measurement.

## What is next, in order

1. **The template ladder.** Identical eight conditions with 1ALU supplied as a
   template for chain A, so the target cannot drift. One variable changed.
   Acceptance check before reading anything else: `targetonly` pLDDT must rise
   and its spread must shrink. Then the real question — does `random40` fall
   below `wt`.
2. **AF2 with the initial-guess protocol.** The calibrated filter, and it solves
   target drift by construction since the designed complex is the starting
   structure. Weights are already on the cluster and the AlphaFold source ships
   with `dl_binder_design`; the remaining cost is a JAX-CUDA and PyRosetta
   install, confirmed necessary because PyRosetta is the script's I/O layer on
   every input path.
3. **The ProteinMPNN sensitivity test.** Same `short_9` backbone, chain A as WT
   versus W157 mutated to ALA in the PDB, 48 sequences each, compare residue
   preferences at the positions contacting 157. The cheapest remaining
   experiment in the project.
4. **Triage the 300 new backbones** (`outputs/*il6_n150_*.pdb`, jobs
   172779/172780) through the Rg gate and the site III contact counter.
5. **Check E93 and K120 are genuinely clear of the binder** with the contact
   counter, and discard either control that is not.
6. **Redo stage 1 with the ROG guiding potential**, which is what the 80%
   attrition was telling us.

**Parked, not dropped:** an alchemical free energy arm, tryptophan to alanine.
It is sensitive by construction and would turn the study from "does the model
notice" into "where does sensitivity enter the pipeline", but it is weeks of
setup. A side-chain-aware empirical ddG is the cheap stand-in to try first.

---

# State of play, end of 2026-10-28

## Stage 3 is still the blocker, and both routes to fixing it moved today

**Route A, pin Boltz with a template (job 174205).** Submitted, reported
`COMPLETED` with exit `0:0` in 2:14, and **7 of its 8 conditions had failed.**
Only `binderonly` ran, which is the only condition carrying no template. Cause:
`template_id: [A1]` in the YAML, where `il6_clean.pdb` has a single polymer
subchain in chain A, so the lookup indexed past a one-element list
(`Error: list index out of range`). Fixed by deleting `chain_id` and
`template_id` and letting Boltz match automatically. Not yet resubmitted.

**Route B, AF2 initial guess.** Environment fully verified: `venv_af2` with
JAX 0.4.20, PyRosetta 2026.39, numpy held at 1.26.4, and **JAX sees the GPU**
(`DEVICES [cuda(id=0)]`, `SUM 262144.0`, job 174258). Runner `run_af2.sh`
written with two arms, the real run and a free `-no_initial_guess` control.
**Not yet submitted** — see the chain-order problem below.

## The catch that would have invalidated the AF2 run

`predict.py` takes `binderlen` from `pose.split_by_chain()` and slices
everything by it: `plddt_binder = plddt[:binderlen]`, and at line 135
`residue_mask = [i > binderlen ...]`. **The first chain in the file must be the
binder.** RFdiffusion writes the target first, so submitting
`af2_in/short9_threaded.pdb` as it stands would have made IL-6 the "binder",
inverted the mask that drives the initial guess, and produced a complete score
file with no error of any kind. Caught by reading the code before submitting.
Fix: reorder the ATOM records, binder chain first, letters unchanged.

## Three process findings worth more than the day's science

1. **An exit code lies in both directions.** Job 172280 said FAILED when the
   work had succeeded; 174205 said COMPLETED when 7 of 8 conditions died. Read
   success off the output file count and the elapsed time against a known rate.
2. **A tool's own error handling can defeat your wrapper's.** Boltz catches a
   bad input, prints `Skipping`, and exits 0, so `grep -c "CONDITION FAILED"`
   returned 0 while seven conditions were failing.
3. **Delete a failing field rather than guess its right value.** Two guesses at
   the Boltz template cause were wrong before reading the schema.

## Unblocked and waiting, all CPU-only

MPNN sensitivity test (wild type against W157A, the one stage never audited),
PyRosetta `InterfaceAnalyzer` and a W157A ddG, and the 2x2 over binder-length
policy and the ROG potential. `dgxanode03` is confirmed healthy, so
`--exclude=dgxanode03` should come out of every script and the pool is 24 GPUs
again, not 16.
