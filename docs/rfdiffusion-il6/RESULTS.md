# IL-6 binder design — run log

Target: human interleukin-6, PDB **1ALU**, 1.90 Å X-ray.
Pipeline: RFdiffusion on the Thapar AI Foundry cluster, `~/rfd/`.
Cluster setup and the DGL/venv details are in `../rfdiffusion-trypsin/README_cluster.md`.

> **This file is a chronological log, not a summary.** It records what was
> believed at each point, including two conclusions later shown to be wrong,
> with the corrections marked in place. For what is currently true, read
> `README.md` (index and state of play), `STAGE3-boltz-scoring.md` (the scoring
> result) and `METHODS.md` (how to run and measure each stage).

---

## The target structure

| | |
|---|---|
| PDB | 1ALU |
| Resolution | 1.90 Å |
| Chains | one, chain A |
| Residues resolved | 157, numbered **19–184** |
| Missing | 1–18 (N-terminal tail), **52–60** (surface loop) |
| Non-protein | 120 HOH, 20 SO4 — crystallisation buffer only |
| Cleaned file | `~/rfd/inputs/il6_clean.pdb` |

Residues 52–60 are a **real chain break**, not a numbering convention. (Trypsin's
six gaps in 2PTN were the opposite case — chymotrypsinogen numbering, physically
continuous — so that target stayed a single contig.) The IL-6 target therefore
has to be given as **two segments**.

Helix boundaries, from the `HELIX` records in the deposited file:

| Helix | Residues |
|---|---|
| A | 21–47 |
| (short) | 69–71 |
| B | 80–104 |
| C | 107–129 |
| CD-loop helix | 141–152 |
| D | 156–182 |

---

## Run 1 — unconstrained binder, no hotspots

Jobs 170908 / 170909, 2026-10-22. (Submitted twice by accident; both used the
same `inference.output_prefix`, so the second overwrote the first. Use distinct
prefixes for parallel jobs.)

```
contigmap.contigs=[A19-51/A61-184/0 70-100]
inference.ckpt_override_path=.../Complex_base_ckpt.pt
diffuser.T=20, num_designs=1
```

| | |
|---|---|
| Runtime | 0.38 min on one A100 |
| Output | `~/rfd/outputs/il6_binder_nohotspot_0.pdb`, 67 KB |
| Chain A | 157 residues — target, numbering preserved (19–184) |
| Chain B | **99 residues — designed**, numbered 1–99 |
| Sampled motif RMSD | 0.10 Å — target held essentially exactly in place |

The two-segment contig worked; RFdiffusion raised no complaint about the gap.
The warning *"model trained on complexes and hotspot residues, without
specifying hotspots"* is expected — none were given on purpose.

### Measuring the interface

The output is **backbone only** (N, CA, C, O — no side chains), so normal
interface definitions mislead:

| Cutoff | Result | Verdict |
|---|---|---|
| CA–CA < 12 Å | 83 residues | too loose — reaches into the core |
| all-atom < 5 Å | 3 residues | too tight — no side chains to bridge the gap |
| **CA–CA < 10 Å** | **58 residues** | right scale for a bare backbone |

### Where it landed

Contacts come in seven consecutive runs, which is what a real interface looks
like rather than scattered noise:

| Contacts | Location |
|---|---|
| 24–36 | helix A, first half |
| 46–48 | end of helix A |
| 64–66 | AB loop |
| 87–98 | helix B, middle |
| 107–128 | **helix C, essentially all of it** |
| 155–158 | start of helix D — **includes Trp157** |
| 166–174 | helix D, middle |

**Reading.** With no guidance at all, the model chose IL-6's concave functional
face over the flat back of the bundle, and its contacts include **Trp157, the
known hot spot of site III**, which contributes roughly 25% of the buried
surface there ([Metcalfe et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7378365/)).
Site I binds IL-6Rα, site II binds the first gp130 at its CHR, site III binds
the second gp130 through its Ig-like D1 domain.

**But the design is too greedy.** 58 contact residues is about double a normal
protein–protein interface. A 99-residue binder with no target site specified has
sprawled across sites II and III at once instead of gripping one properly.

---

## Status

Nothing has been scored. This is a backbone with no sequence, so it is a shape
rather than a protein. Still ahead:

1. **ProteinMPNN** — write a sequence for the designed backbone.
2. **AF2 / Chai re-prediction** — fold that sequence back and score the complex
   on interface ipTM / PAE. This is the real filter and most designs fail it.

## Next run (planned)

Hotspots around Trp157 and a shorter binder (55–75 residues), to watch the
interface tighten onto site III alone. Same job script, two changed values.

---

## Run 2 — site III hotspots

Script: `~/rfd/31_rfd_il6_hotspot.sbatch`

```
contigmap.contigs=[A19-51/A61-184/0 55-75]
ppi.hotspot_res=[A47,A155,A157,A158]
diffuser.T=20, num_designs=1
```

Hotspots were picked from the data rather than from literature memory: every
residue within 12 Å of Trp157's alpha carbon, intersected with the run-1 contact
set so that only demonstrably surface-reachable residues survived. That left
46, 47, 48 and 155, 156, 157, 158.

Output: `~/rfd/outputs/il6_binder_site3_0.pdb`, 57 KB. Chain B = **60 residues**.

### Prediction, written down before measuring

> Contacts drop from 58 to roughly 25–35, collapsing onto helix D and the tail
> of helix A, losing helix B and most of helix C.

### Result: half right

| Region | Run 1 (no hotspots) | Run 2 (hotspots) |
|---|---|---|
| helix A, first half | 10 | 4 |
| helix A, tail (hotspot A47) | 3 | 6 |
| AB loop | 3 | 3 |
| helix B | 11 | 11 |
| helix C | 18 | 4 |
| helix D, start (hotspots A155/157/158) | 4 | 7 |
| helix D, middle | 9 | 9 |
| **Total** | **58** | **44** |

**Right:** the interface moved onto the hotspots. Helix C collapsed 18 → 4,
helix A's first half 10 → 4, and both hotspot patches grew.

**Wrong:** the total was 44, not 25–35, and helix B did not move at all.

**The number that matters:** contacts per binder residue, 0.59 → **0.73**, up a
quarter. The binder shrank 39% while its interface shrank only 24%, so hotspots
made it grip harder per residue spent. Raw contact count hides this; use the
ratio when comparing runs.

**Open question:** helix B may simply be spatially adjacent to site III, in which
case those contacts are a fact about IL-6's shape rather than a focusing failure.
Testable with a CA–CA distance between residues 157 and 92.

### Helix B question — answered, and the guess was wrong

CA–CA distance from residue 157 to residue 92: **26.7 Å**. That is most of the
way across IL-6, so helix B is *not* adjacent to site III.

The binder is therefore genuinely elongated, stretching from site III over to
helix B rather than sitting compactly on the hotspots it was given. Hotspots
focused it; they did not make it compact.

Next diagnostic: tighten the contact cutoff from 10 Å to 8 Å and see which
regions survive, to separate real grip from grazing contact.

### 8 Å cutoff — the binder has two footholds

21 residues survive the tighter cutoff:

| Block | Residues | Where |
|---|---|---|
| 46–49 | 4 | helix A tail — **hotspot A47** |
| 65 | 1 | AB loop |
| 90–98 | 7 | **helix B — no hotspot here** |
| 119 | 1 | helix C |
| 155–157 | 3 | helix D start — **hotspots A155, A157** |
| 169–174 | 5 | helix D middle |

Helix B survives with the largest single block, so it is a real grip rather than
grazing contact. The binder holds IL-6 in **two places 27 Å apart**: the site it
was aimed at, and an unrequested foothold on helix B.

The fix is geometric. A 40–55 residue mini-protein cannot span 27 Å, so it has to
choose one site.

**Method note for future runs:** measure at two cutoffs — 10 Å for the overall
footprint, 8 Å for the real grip — and compare *contacts per binder residue*
between runs rather than raw counts.

---

## ⚠️ CONTACT MEASUREMENTS ABOVE ARE UNDER INVESTIGATION (2026-10-22)

**Do not rely on the contact counts in runs 1 and 2 above, or on the
"two footholds" conclusion, until this is resolved.**

Run 3 (10 designs, 40–55 residues, T=50, same site III hotspots) produced ten
binders whose centroids differ substantially (x from +8.7 to −10.8), yet every
one of them reports contacts with the *same seven* helix B residues at the *same
distances*: 90 (7.0), 91 (5.0), 92 (8.0), 93 (7.8), 94 (4.8), 95 (7.1), 98 (7.5).

Ten independent stochastic designs cannot produce identical contact distances.
Something in the measurement is wrong.

Hypotheses tested and **refuted**:
1. Helix B sits behind hotspot A47 — no, residue 47 is 18–27 Å from 90–98.
2. A phantom chain B residue numbered 0 — no such residue exists; chain B is
   numbered 1..N cleanly, and excluding r=0 changed nothing.

Unexplained: an earlier diagnostic printed "nearest B residue=0" when chain B
has no residue 0, which suggests a variable was unset and the distance may have
been measured from the wrong origin.

Outstanding checks: (a) md5 of each file's chain B CA coordinates, to confirm
the ten binders really differ; (b) whether target residue 94 is being found at
all in the parse.

**No design proceeds to ProteinMPNN until this is understood.**

---

# ✅ CORRECTED RESULTS (2026-10-22)

**The bug is found and fixed.** Every contact number earlier in this file was
inflated by a phantom atom at the coordinate origin, caused by an uninitialised
array counter in the awk parser. Full write-up: `CAUTION-origin-atom-bug.md`.

**Everything above this line is superseded.** The table below is measured with
`BEGIN { na=0; nb=0 }` and is the trustworthy record.

`c10` / `c8` = target residues within 10 Å / 8 Å (CA–CA) of the binder.
`siteIII` = of the `c8` contacts, how many are in 44–49 or 153–182.
`grip` = `c8` divided by binder length.

| Design | len | c10 | c8 | helixB | siteIII | grip |
|---|---|---|---|---|---|---|
| run 1, no hotspots | 99 | 27 | 10 | 0 | 2 | 0.10 |
| run 2, hotspots | 60 | 13 | 7 | 0 | **7** | 0.12 |
| short_0 | 45 | 16 | 7 | 0 | 3 | 0.16 |
| short_1 | 55 | 27 | 12 | 0 | 6 | 0.22 |
| short_2 | 43 | 7 | 5 | 0 | 5 | 0.12 |
| short_3 | 49 | 14 | 5 | 0 | 4 | 0.10 |
| short_4 | 55 | 11 | 3 | 0 | 1 | 0.05 |
| short_5 | 45 | 12 | 6 | 0 | 3 | 0.13 |
| short_6 | 42 | 6 | 2 | 0 | 2 | 0.05 |
| **short_7** | 53 | 28 | **16** | 0 | **9** | **0.30** |
| short_8 | 53 | 9 | 2 | 0 | 2 | 0.04 |
| short_9 | 40 | 21 | 8 | 0 | 6 | 0.20 |

## What the corrected numbers actually say

**1. helixB = 0 in every single structure.** The foothold was entirely the
artifact. No binder ever touched helix B.

**2. Hotspots worked, and worked cleanly.** This is clearer than the corrupted
data ever showed:

| | site III contacts | total contacts | fraction on target site |
|---|---|---|---|
| run 1 (no hotspots) | 2 | 10 | 20% |
| run 2 (hotspots) | 7 | 7 | **100%** |

Without hotspots the binder put a fifth of its grip on site III. With them,
*all* of it. That is the real demonstration of what hotspots buy.

**3. The ten designs vary enormously.** `grip` ranges 0.04 to 0.30, site III
contacts from 1 to 9. Four of the ten (`_4`, `_6`, `_8`, and arguably `_0`)
barely touch the target at all. This is exactly why the method is to generate
many and filter — a single design tells you nothing.

**4. The interfaces are small in absolute terms.** Note that CA–CA on a
side-chain-free backbone understates contact, so these are not directly
comparable to published interface sizes. AF2 scoring is the proper judge.

## Ranked for stage 2

By site III contacts, then grip:

1. **short_7** — 16 contacts, 9 on site III, grip 0.30. Clear standout.
2. **short_1** — 12 contacts, 6 on site III, grip 0.22.
3. **short_9** — 8 contacts, 6 on site III, grip 0.20.
4. **short_2** — 5 contacts, 5 on site III (100% on target), grip 0.12.
5. **run 2 original** — 7 contacts, all 7 on site III, grip 0.12.

---

# ⚠️ The contact ranking was selecting for a failure mode (2026-10-22)

ProteinMPNN on `short_7` returned alanine-saturated sequences, e.g.

```
SEDEKRLEEEEKEEEAAAKAGREKNGAATAAAAAALAAAAAAAALAAAAAALA
```

Poly-alanine is a known signature of a backbone with **no hydrophobic core** —
there is nothing to bury, so the model picks the small residue that cannot
clash. Chain assignment was correct (`fixed_chains=['A'], designed_chains=['B']`),
so this is about the backbone, not the setup.

## Radius of gyration confirms it

A folded mini-protein of N residues has Rg ≈ 2.2·N^0.38. Measured:

| Design | n | Rg (Å) | compact ≈ | ratio | verdict |
|---|---|---|---|---|---|
| nohotspot_0 | 99 | 26.0 | 12.6 | **2.07** | extended |
| site3_0 | 60 | 14.7 | 10.4 | 1.41 | loose |
| short_0 | 45 | 19.5 | 9.3 | **2.10** | extended |
| short_1 | 55 | 23.7 | 10.1 | **2.35** | worst |
| **short_2** | 43 | 10.5 | 9.2 | **1.14** | **compact** |
| short_3 | 49 | 12.1 | 9.7 | 1.25 | acceptable |
| short_4 | 55 | 12.9 | 10.1 | 1.28 | acceptable |
| short_5 | 45 | 19.5 | 9.3 | **2.10** | extended |
| short_7 | 53 | 20.4 | 9.9 | **2.06** | extended |
| short_8 | 53 | 13.2 | 9.9 | 1.33 | acceptable |
| **short_9** | 40 | 10.8 | 8.9 | **1.21** | **compact** |
| short_6 | 42 | 11.1 | 9.1 | 1.22 | compact but only 2 contacts |

## The lesson: contacts and compactness are anti-correlated here

| Design | contacts (c8) | Rg ratio |
|---|---|---|
| short_1 | 12 | 2.35 |
| short_7 | 16 | 2.06 |
| short_0 | 7 | 2.10 |
| short_9 | 8 | 1.21 |
| short_2 | 5 | 1.14 |

**A ribbon draped across the target racks up contacts precisely because it is
spread out.** Ranking on contact count therefore selected for the designs least
likely to fold. `short_7`, ranked first on contacts, is 2.06× too extended to be
a protein.

## Corrected ranking

Compactness is a gate, not a tiebreaker — an unfoldable binder scores zero no
matter how much surface it touches.

1. **short_9** — Rg ratio 1.21, 8 contacts, 6 on site III. Compact *and* gripping.
2. **short_2** — Rg ratio 1.14, 5 contacts, all 5 on site III. Most compact.
3. short_3 — ratio 1.25, 5 contacts, 4 on site III.
4. short_4 — ratio 1.28, 3 contacts.
5. short_8 — ratio 1.33, 2 contacts.

**Dropped:** short_7, short_1, short_0, short_5 and run 1 — all extended.

## Fix for future RFdiffusion runs

Add a radius-of-gyration guiding potential so compactness is enforced during
diffusion rather than filtered afterwards:

```
potentials.guiding_potentials=["type:binder_ROG,weight:3"]
potentials.guide_scale=2
potentials.guide_decay="quadratic"
```

---

# Stage 2 — ProteinMPNN on the compact designs (2026-10-23)

Job **171043**, `40_mpnn.sbatch` with `--gres=gpu:1`. Five compact designs
(`short_9`, `_2`, `_3`, `_4`, `_8`), 8 sequences each, `--sampling_temp 0.1`,
`--seed 37`. About **4.7 seconds per structure** on the GPU, against ~11 minutes
per structure on CPU in the same container.

## The test, written down beforehand

> Sequences with a proper mix — leucine, isoleucine, valine, phenylalanine,
> tryptophan scattered through, not a wall of alanine. Those hydrophobic
> residues are the ones that form a core. If they show up, the compactness
> diagnosis was right. If it is still alanine on a compact backbone, the cause
> is something else and we stop before spending anything on AF2.

## Result: the diagnosis was right

`short_9`, 40 residues, all eight samples:

```
NEYFEYLNELVAKAEEIAKQENAEELEKKLKEELEKKLKE
SSYDEYLAQLVAQAEELAKTADAAAKLAALQAELAERLKK
SAYDEKLAETVAKAEELAKQANAEELLAALKKELEEELKK
SEYEEKMAQTVAEAFKLQQTADAEAKLAALKAAIQKELEA
SSYFQALADRVALAEELAKQANAEALLAALQAQLQAELAA
SSYEEYLAELVAKAEEYSKQEDAEKLLEEMKKKLEELLKK
SAYEEELAETVAKAFELAETEDAEEKLAALKAEIEKKLAE
SSYEESLAELVAKAEEIAKQENAEELLKKLKEELEKKLKE
```

Nothing was changed except the backbone. Same runner script, same
`fixed_chains=['A'], designed_chains=['B']`, same temperature and seed as the
run that produced poly-alanine on `short_7`. So the alanine wall was a property
of that backbone, not of the pipeline. Confirmed.

## The sequences are better than merely "mixed"

Consensus across all eight independent samples:

| position | 3 | 7 | 11 | 17 | 26 | 30 | 34 | 38 |
|---|---|---|---|---|---|---|---|---|
| residue | Y (8/8) | L or M (8/8) | **V (8/8)** | I/L/Y (8/8) | L | L | L | L |

Positions 3, 7, 11 are *i*, *i*+4, *i*+8 — a single face of an alpha helix.
Positions 26, 30, 34, 38 are the same pattern on a second helix. Eight
stochastic samples agreeing at those positions means ProteinMPNN located a
hydrophobic core: **two amphipathic helices packing their buried faces against
each other, with the charged residues on the solvent side.** That is a textbook
helical hairpin and exactly the right architecture for a 40-residue
mini-protein.

The surfaces are very E/K rich (`EELEKKLKEELEKKLKE` in several samples). For a
designed binder that is normal and desirable — it keeps the protein soluble —
but it does mean the binding face could be mostly polar, which the sequence
alone cannot tell us.

## Two numbers not to over-read

**`seq_recovery = 0.0000`.** This is meaningless for a de novo backbone. The
"native" sequence of an RFdiffusion output is poly-glycine, visible in the first
line of every `.fa` file, so recovery against it is zero by construction.
`short_7`'s earlier 0.038 was alanine-against-glycine noise, not a quality
signal.

**ProteinMPNN score 1.30 – 1.37** (negative log-likelihood; lower is better,
and the poly-glycine input scores 2.51). Middling rather than good. These are
plausible sequences, not proven ones.

## What this does and does not establish

Established: these backbones **can** fold. A sequence with a clean amphipathic
core is the prerequisite.

Not established: that they **bind**. The core points inward and IL-6 is outward,
so the folding evidence says nothing about the interface.

## Next

1. Which binder residues lie within 10 Å of IL-6, and what ProteinMPNN put at
   those positions. Cheap, and it checks that the designed interface is not all
   charge. Prediction recorded before measuring: 8–16 of the 40 residues.
2. **Stage 3**: re-predict each sequence in complex with IL-6 and score on
   interface ipTM / PAE. Planned on **Boltz rather than AF2**, because Boltz
   installs with plain `pip` into the existing venv and this cluster has no
   conda, which rules out the localcolabfold installer.

---

## Where the binder actually touches IL-6 (`short_9`, 2026-10-23)

Binder residues whose CA is within 10 Å of any IL-6 CA. Prediction recorded
first: 8 to 16 of 40. **Actual: 18.** Slightly under-predicted.

```
3 4 5 6 7 8 9 10 11 12 13 14 15   18   27 31 34 38
```

The shape matters more than the count:

| Block | Pattern | Reading |
|---|---|---|
| 3–15 | thirteen consecutive | helix 1 lies flat along IL-6 |
| 18 | single | end of helix 1 |
| 19–26 | **absent** | the turn, held clear of the target |
| 27, 31, 34, 38 | *i*, *i*+4, *i*+7, *i*+11 | one face of helix 2 |

A helical hairpin lying on site III with helix 1 flat on the surface, helix 2
angled across touching along a single stripe, and the turn pointing away. As
binding geometry that is sensible.

### The tension this creates, unresolved

The positions ProteinMPNN conserved as hydrophobic (Y3, L7, V11, and L at 34 and
38) are **the same positions that are at the interface**. So those residues may
be buried against IL-6 rather than against each other.

- Good for binding: buried hydrophobic surface is what pays for affinity.
- Open for folding: if the greasy residues all point at the target, the binder
  may have no core of its own when IL-6 is absent.

**This measurement cannot settle it, and a tighter cutoff would not either.**
10 Å CA–CA is wider than a helix is thick, so for a helix lying on a surface
both of its faces fall inside the cutoff. That is precisely why all thirteen
consecutive residues scored positive, and it means the method has no power to
say which face points at the target.

The resolution is stage 3, which answers both halves directly: predict the
binder alone (pLDDT, does it fold) and the complex (ipTM and interface PAE,
does it bind).

---

## Stage 3 — first Boltz scoring of `short_9` (job 172278, 2026-10-23)

**Input.** `~/rfd/boltz_in/short9_wt.yaml`. Chain A = IL-6, PDB 19–184, 166 aa.
Chain B = `short_9` ProteinMPNN sample 2 (best of eight, score 1.3046), 40 aa:
`SSYDEYLAQLVAQAEELAKTADAAAKLAALQAELAERLKK`. Both chains `msa: empty`, i.e.
single-sequence mode.

**Run.** Boltz 2.2.1, `--diffusion_samples 5 --write_full_pae --no_kernels`,
one A100-SXM4-80GB. **17 seconds for all five samples**, 0 failed examples.
Boltz has no seed flag, so the five diffusion samples *are* the replicates: one
job gives the estimate and its noise floor together.

### The numbers

| metric | mean | SD | SEM | min | max | range |
|---|---|---|---|---|---|---|
| `iptm` (= `protein_iptm`) | 0.581 | 0.116 | 0.052 | 0.489 | 0.724 | 0.235 |
| `pair_chains_iptm` [0][1] | 0.496 | 0.152 | 0.068 | 0.336 | 0.670 | 0.333 |
| `chains_ptm` A (IL-6) | 0.684 | 0.144 | 0.064 | 0.489 | 0.795 | 0.306 |
| `chains_ptm` B (binder) | 0.875 | 0.032 | 0.014 | 0.835 | 0.913 | 0.078 |
| `complex_plddt` | 0.685 | 0.043 | 0.019 | 0.632 | 0.727 | 0.095 |
| `complex_iplddt` | 0.641 | 0.037 | 0.017 | 0.595 | 0.672 | 0.077 |
| `complex_ipde` | 5.354 | 1.840 | 0.823 | 3.261 | 7.255 | 3.994 |

Note `pair_chains_iptm[1][0]` equals `iptm` exactly in all five samples, while
`[0][1]` is a different, systematically lower number. The matrix is not
symmetric; do not average the two.

### Predictions written down before the measurement, and how they fared

| prediction | actual | verdict |
|---|---|---|
| interface ipTM 0.2–0.5 | 0.581, up to 0.724 | **wrong**, too pessimistic |
| binder confident on its own, >0.8 | `chains_ptm` B = 0.875 | right, but on pTM, not the pLDDT I named |
| spread across samples ±0.05 | SD 0.116 | **wrong, by more than 2x** |

Both misses point the same way: the model is more generous and far less
self-consistent than assumed.

### KEY FINDING — the noise floor is the result, and it blocks the plan as drafted

`iptm` 0.581 ± 0.116 means the five samples span **0.489 to 0.724**, a range of
0.235. The mutation audit asks whether a target hotspot mutation *moves* an
interface metric. A W→A substitution at a real hotspot would plausibly shift
ipTM by 0.05–0.15. At n = 5 the standard error is 0.052, so **effects of that
size are inside the noise and cannot be resolved.**

The remedy is cheap rather than clever. At 3.4 s per diffusion sample:

- n = 15 → SEM 0.030
- n = 34 → SEM 0.020
- n = 50 → SEM 0.016, about **three minutes** of A100

So raise `--diffusion_samples` to 50 before running any comparison. This is the
single most important consequence of the first run.

### SECOND FINDING — the target itself is not stably folded

`chains_ptm` for chain A, which is **IL-6, a 1.9 Å crystal structure we handed
the model**, ranges 0.489 to 0.795. In single-sequence mode Boltz cannot fold
the target reliably, so in the worst samples the "interface" being scored is
between the binder and a mis-folded IL-6.

Target confidence and interface confidence move together, Pearson r = 0.670
(n = 5): suggestive and nothing more. Sample `model_2` breaks the pattern, with
a well-folded target (0.769) and a poor interface (0.489), so target failure is
not the whole story.

`complex_ipde` is strongly anti-correlated with `iptm`, r = −0.854 (n = 5), as
it should be: the two are not independent readouts and should not be treated as
two pieces of evidence.

**Remedy to try:** pin the target. Either give chain A an MSA (generate it on
the login node, which resolves external hosts, and reference the `.a3m` from the
YAML) or supply 1ALU as a template for chain A. If IL-6 folds the same way every
time, a large part of this variance should disappear — and that is a testable
prediction, not an assumption.

### What this does NOT say

It does not say `short_9` binds IL-6. Boltz's calibration as a binder filter is
unpublished, so 0.581 has no threshold to be compared against, and co-folding
models are known to be optimistic on designed interfaces. The negative control,
a random 40-mer scored the same way, is what will give this number a scale.

### Verified measurement helper

Mean and SD straight off the confidence files, checked against the five values
above before use (it reproduces mean 0.581, SD 0.116, SEM 0.052 exactly):

```bash
grep -h '"iptm"' boltz_out/*/predictions/*/confidence_*.json | tr -d ' ,"' | cut -d: -f2 > /tmp/iptm.txt
awk '{s+=$1;q+=$1*$1;n++} END{m=s/n;printf "n=%d mean=%.3f sd=%.3f sem=%.3f\n",n,m,sqrt((q-n*m*m)/(n-1)),sqrt((q-n*m*m)/(n-1))/sqrt(n)}' /tmp/iptm.txt
```

`grep '"iptm"'` does not match `"ligand_iptm"` or `"protein_iptm"`, because both
have an underscore rather than a quote before `iptm`. Confirmed on a mock file.

---

## The n=30 hotspot comparison (jobs 172281, 172282, 2026-10-25)

Both arms COMPLETED, 20 minutes each, 30 designs each, binder length fixed at
exactly 50 so size cannot confound the comparison. **No radius-of-gyration
potential**, by a deliberate shortcut: the plan called for it in both arms, and
it was dropped to avoid a nested-quoting fight in the job script, with the
intention of applying the same Rg gate afterwards instead.

### FINDING — 44 of 60 designs are single straight helices

The compactness gate returned Rg ratios clustered at **2.222 with a standard
deviation of 0.023** across 44 of the 60 independent designs.

By rule 2 of `CAUTION-origin-atom-bug.md`, a quantity that is near-invariant
across independent stochastic samples is a bug until proven otherwise, so it was
checked rather than reported. For *n* equally spaced collinear points at spacing
*d*, Rg = sqrt(d²(n²−1)/12). For n = 50 and an α-helix rise of d = 1.5 Å that is
**21.65 Å**. The cluster's mean Rg is **21.62 Å**, range 20.82–21.89.

So this is not a measurement artefact. Those designs are **single straight
50-residue α-helices**, and their Rg is tightly reproducible precisely because a
straight rod of fixed length is geometrically determined. A real failure mode,
diagnosed by the arithmetic rather than by eye.

### Yield, and the cost of the shortcut

| arm | designs | pass Rg ≤ 1.30 | pass rate |
|---|---|---|---|
| hotspots `A47,A155,A157,A158` | 30 | 6 | 20% |
| no hotspots | 30 | 5 | 17% |

Passing ratios, hotspot arm: 1.10, 1.20, 1.23, 1.23, 1.26, 1.28.
No-hotspot arm: 1.08, 1.23, 1.23, 1.29, 1.29.

**Prediction written down beforehand: 45% attrition, about 14 usable per arm.
Actual: 80% attrition, 6 and 5.** Wrong by a factor of four, and the cause is
identifiable: without a compactness potential during diffusion, a fixed binder
length of 50 residues makes a single long helix the model's default answer.
Equal pass rates across the two arms is expected and reassuring, since hotspot
conditioning has no reason to affect compactness.

**Consequence: the shortcut turned n = 30 into n = 6**, which is barely better
than the n = 1 the run existed to fix. The radius-of-gyration potential is not a
refinement, it is load-bearing. Two ways out, and the brute-force one was taken
because it carries no new failure modes: 150 designs per arm on the identical
scripts yields roughly 30 usable per arm for about 100 minutes of A100 each.
Doing it properly with `potentials.guiding_potentials=["type:binder_ROG,weight:3"]`
remains the better answer and is the next cleanup.

### Also on 2026-10-25: job 172280, and a false alarm worth recording

Job 172280, `boltzrun`, showed FAILED with exit 1 after 1:58, and `boltz_out`
still reported n = 5. The first reading was that fifty diffusion samples had
broken something five did not, and that no high-sample batch should be
submitted until it was understood.

**That reading was wrong, and the log said so.** The error file contained no
traceback, no memory error and nothing from the model — only a failing `cat`,
and the path it tried was `boltz_out/predictions/short9_wt/...`, the *old*
incorrect path that exists only in the unmodified `run_boltz.sh`. So 172280 was
the original 5-sample script run a second time (1:58 rather than 172278's 4:44
because the CCD cache was already warm), and its exit code came entirely from
that harmless final command. The 50-sample edit was never submitted.

Two lessons. A non-zero exit code from a script whose **last** line is a
diagnostic `cat` says nothing about whether the real work succeeded; put the
summary command behind `|| true`, or do not let it be last. And the path in an
error message identifies *which version of a script ran*, which is often faster
evidence than any job record.

The eight-condition ladder was submitted at **50 samples per condition** rather
than 100: SEM 0.016 is already four times better than the n = 5 that blocked
us, and it halves the exposure if high sample counts do turn out to be a
problem. The loop uses `|| echo "CONDITION FAILED: $y"` and no `set -e`, so one
bad condition cannot kill the remaining seven.

---

## Length audit of run 3 — the confound was already handled (2026-10-28)

Prompted by `Sequence init` in the job logs: `rfd_il6_short_170960.out` reports
`masked=45` while `short_9` is 40 residues. Resolution: the run-3 contig was a
**range, 40–55**, so `Sequence init` reports design 0's sampled length (45,
which is `short_0`) rather than the contig. The log line describes one design,
not the job.

Measured directly from the output PDBs, chain B CA count per design:

| design | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| nB | 45 | 55 | 43 | 49 | 55 | 45 | 42 | 53 | 53 | 40 |

**These match the `len` column of the corrected table above, design for
design.** The varying length was recorded correctly at the time and nothing in
the stage-1 record is wrong.

### The concern that was raised, and why it does not survive

Concern: with lengths spread 40–55, ranking on raw contact counts partly ranks
by chain length. Checked by normalising site III contacts by binder length:

| design | siteIII | len | per residue |
|---|---|---|---|
| short_7 | 9 | 53 | 0.170 |
| short_9 | 6 | 40 | 0.150 |
| run 2 | 7 | 60 | 0.117 |
| short_2 | 5 | 43 | 0.116 |
| short_1 | 6 | 55 | 0.109 |
| short_3 | 4 | 49 | 0.082 |

The top of the order is unchanged; `short_9` gains one place over `short_1`.
And the decision that actually mattered never used a count at all: `short_7`
was rejected on **Rg ratio 2.06**, which is length-normalised by construction
(`2.2·n^0.38`). The `grip` column was already `c8/len`. **No revision required.**

### What the audit did turn up: pass rate against binder-length policy

| run | length policy | pass Rg ≤ 1.30 | rate |
|---|---|---|---|
| run 3 (job 170960) | range 40–55, T=50 | 6/10 | 60% |
| job 172779 | fixed 50, hotspots | 47/150 | 31% |
| job 172780 | fixed 50, no hotspots | 13/150 | 9% |

60% against 31% gives z = 1.9, p ≈ 0.06. **Suggestive and underpowered**, and
confounded: run 3 also used T=50 and differs from the n=150 arms in more than
the length policy. But it is consistent with the straight-helix mechanism — a
model free to choose its own length need not stretch a fixed residue count into
a rod.

**Cheap test, and it belongs in the ROG cleanup already planned:** one job, 2×2
over length policy (fixed 50 against range 40–55) and the ROG potential on or
off, 50 designs per cell. Rg ratio is already length-normalised, so the cells
are directly comparable.
