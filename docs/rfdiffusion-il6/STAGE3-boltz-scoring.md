# Stage 3 — first scoring of the `short_9` IL-6 binder

Boltz 2.2.1, job 172278, 2026-10-23. Self-contained; the chronological version
is in `RESULTS.md`. Index and project framing: `README.md`.

## What was scored

| | |
|---|---|
| chain A | IL-6, PDB 1ALU, residues 19–184, **166 aa** |
| chain B | `short_9`, ProteinMPNN sample 2, best of eight at score 1.3046, **40 aa** |
| chain B sequence | `SSYDEYLAQLVAQAEELAKTADAAAKLAALQAELAERLKK` |
| MSA | `msa: empty` on both chains, i.e. **single-sequence mode** |
| sampling | 5 diffusion samples, `--write_full_pae --no_kernels` |
| hardware | one A100-SXM4-80GB |
| wall time | **17 seconds for all five samples**, 0 failed examples |

Boltz has no seed flag, so the diffusion samples *are* the replicates. One job
gives the estimate and its noise floor together, at about 3.4 s per sample.

## The numbers

| metric | mean | SD | SEM | min | max | range |
|---|---|---|---|---|---|---|
| `iptm` (= `protein_iptm`) | **0.581** | **0.116** | 0.052 | 0.489 | 0.724 | 0.235 |
| `pair_chains_iptm` [0][1] | 0.496 | 0.152 | 0.068 | 0.336 | 0.670 | 0.333 |
| `chains_ptm` A, IL-6 | 0.684 | 0.144 | 0.064 | 0.489 | 0.795 | 0.306 |
| `chains_ptm` B, binder | 0.875 | 0.032 | 0.014 | 0.835 | 0.913 | 0.078 |
| `complex_plddt` | 0.685 | 0.043 | 0.019 | 0.632 | 0.727 | 0.095 |
| `complex_iplddt` | 0.641 | 0.037 | 0.017 | 0.595 | 0.672 | 0.077 |
| `complex_ipde` | 5.354 | 1.840 | 0.823 | 3.261 | 7.255 | 3.994 |

Per-sample `iptm`: 0.724, 0.690, 0.489, 0.497, 0.504.

**`pair_chains_iptm` is not symmetric.** `[1][0]` equals `iptm` exactly in all
five samples; `[0][1]` is a different and systematically lower number. Do not
average the two halves.

## Predictions written down before the measurement

The project rule is to state the expected number first, adopted after a
measurement bug went unnoticed for three rounds.

| prediction | actual | verdict |
|---|---|---|
| interface ipTM 0.2–0.5 | 0.581, up to 0.724 | **wrong**, too pessimistic |
| binder confident on its own, above 0.8 | `chains_ptm` B = 0.875 | right, but on pTM, not the pLDDT named |
| spread across samples ±0.05 | SD 0.116 | **wrong, by more than 2x** |

Both misses point the same way: the model is more generous, and far less
self-consistent, than assumed.

## Finding 1 — the noise floor is the result, and it blocks the plan as drafted

`iptm` 0.581 ± 0.116 means the five samples span 0.489 to 0.724, a range of
0.235 on a metric bounded by 1.

The mutation audit asks whether a target hotspot mutation *moves* an interface
metric. A tryptophan-to-alanine substitution at a real hotspot would plausibly
shift ipTM by 0.05 to 0.15. At n = 5 the standard error is 0.052, so **an effect
of that size sits inside the noise and cannot be resolved.**

The remedy is cheap rather than clever. At 3.4 s per sample:

| samples | SEM on ipTM | GPU time |
|---|---|---|
| 5 | 0.052 | 17 s |
| 15 | 0.030 | ~1 min |
| 34 | 0.020 | ~2 min |
| 50 | 0.016 | ~3 min |

**Raise `--diffusion_samples` to 50 before running any comparison.** This is the
single most important consequence of the first run, and it is the reason to
measure a noise floor before reporting any difference at all.

## Finding 2 — the target itself is not stably folded

`chains_ptm` for chain A ranges **0.489 to 0.795**. Chain A is IL-6: a 1.90 Å
crystal structure that we handed to the model. In single-sequence mode Boltz
cannot fold the target reliably, so in the worst samples the interface being
scored lies between the binder and a mis-folded IL-6.

Target confidence and interface confidence move together, Pearson **r = 0.670**
at n = 5, which is suggestive and nothing more. Sample `model_2` breaks the
pattern outright, with a well-folded target at 0.769 and a poor interface at
0.489, so target failure is not the whole story.

`complex_ipde` is strongly anti-correlated with `iptm`, **r = −0.854**, as it
should be. The two are **not independent readouts** and must not be cited as two
pieces of evidence for the same conclusion.

### Two remedies, and a prediction

1. **Pin the target in Boltz**, either by giving chain A an MSA (generate it on
   the login node, which resolves external hosts, and reference the `.a3m` from
   the YAML) or by supplying 1ALU as a template for chain A.
2. **AF2 with the initial-guess protocol**, which hands the model the designed
   complex as its starting structure. The target is then never predicted from
   sequence and cannot drift. This is the same lever, in the version the field
   has calibrated against experimental hit rates.

Either way the prediction is testable and worth writing down now: **if IL-6
folds the same way every time, most of the 0.116 spread should disappear.** If
it does not, the variance is coming from the binder pose rather than the target,
and that is a different and more interesting problem.

## What this does not say

It does not say `short_9` binds IL-6.

- Boltz has **no published calibration as a binder filter**, so 0.581 has no
  threshold to be compared against. AF2's thresholds are published; they do not
  transfer.
- Co-folding models are known to be optimistic on designed interfaces.
- There is no control yet. A random 40-mer scored the same way is what will give
  this number a scale, and it is the next thing to run after the 50-sample
  rerun.

The honest summary is that the pipeline now produces a number, that the number's
uncertainty is known, and that the uncertainty is the part that needed fixing
first.
