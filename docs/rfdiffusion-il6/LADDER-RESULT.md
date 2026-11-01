# The sensitivity ladder: result

Job 172783, 2026-10-28. Eight conditions, 50 diffusion samples each, 24 minutes
on one A100, zero failed conditions. Design and pre-registered interpretations:
`BATCH-2026-10-25.md`. Index: `README.md`.

**Headline: the go/no-go gate fails. A random 40-mer scores as well as the
designed binder, so interface ipTM in this configuration carries no information
about whether `short_9` binds IL-6, and the mutation audit cannot be run on it
as it stands.**

## The table

| condition | ipTM mean | SD | SEM | diff from `wt` | Welch t | complex pLDDT |
|---|---|---|---|---|---|---|
| `wt` | 0.466 | 0.157 | 0.022 | — | — | 0.665 ± 0.068 |
| `random40` | **0.502** | 0.161 | 0.023 | **+0.036** | +1.13 | 0.638 ± 0.043 |
| `k120a` | 0.414 | 0.195 | 0.028 | −0.052 | −1.46 | 0.655 ± 0.084 |
| `w157a` | 0.394 | 0.097 | 0.014 | −0.072 | −2.76 | 0.651 ± 0.046 |
| `alaif` | 0.393 | 0.146 | 0.021 | −0.073 | −2.40 | 0.616 ± 0.068 |
| `e93a` | 0.376 | 0.143 | 0.020 | −0.090 | −3.03 | 0.600 ± 0.071 |
| `targetonly` | 0.000 | — | — | — | — | 0.667 ± 0.090 |
| `binderonly` | 0.000 | — | — | — | — | **0.905 ± 0.028** |

ipTM is exactly 0.000 for both single-chain conditions, which is the correct
behaviour and a useful check that the pipeline did what was asked.

## 1. The gate fails

`random40` is a sequence drawn at random. It scores **0.502 against the
design's 0.466** — numerically higher, and statistically indistinguishable
(t = +1.13). The interface-dead mutant `alaif` scores 0.393, i.e. **lower than
the random peptide**.

A filter whose negative control outranks its positive cannot be used to choose
designs, and it cannot be used to ask whether a target mutation matters. This
was the pre-registered first branch: *"If `alaif` and `random40` do not drop
clearly below `wt`, the metric is inert on designs of this quality and the audit
cannot run on it."*

`alaif` does drop significantly (−0.073, t = −2.40), so the metric is not
completely insensitive to the binder's sequence. But direction is what matters
for a filter, and a random 40-mer beating the design settles it.

## 2. The target mutations move the metric, but not as a hotspot would

| mutation | chosen to be | diff from `wt` |
|---|---|---|
| W157A | **the site III hotspot** | −0.072 |
| E93A | irrelevant, surface helix B | **−0.090** |
| K120A | irrelevant, surface helix C | −0.052 |

Mean drop across the three: −0.071. **The mutation chosen precisely because it
should not matter produces the largest drop.** This is the second
pre-registered branch: the model responds to *any* perturbation of the target
sequence rather than to the chemistry of the interface.

Taken alone, `w157a` at t = −2.76 looks like a positive result. The two control
mutations are the entire reason we are not reporting it as one. That is what
controls are for, and it is the strongest methodological point in the project so
far.

## 3. The mechanism is visible in pLDDT

Across the five conditions other than `random40`, interface confidence tracks
whole-complex confidence at **r = +0.78**. `e93a` is the lowest on both metrics.

So the chain of causation appears to be: a point mutation makes IL-6 harder for
the model to fold from sequence alone → whole-complex confidence falls →
interface confidence falls with it. The drop is not about the interface at all.

`targetonly` confirms the premise independently: **IL-6 by itself scores
0.667 ± 0.090**, no better than the complex. The target is the weak link
throughout, exactly as the first scoring run suggested.

## 4. Good news: `short_9` folds on its own

`binderonly` gives **pLDDT 0.905 ± 0.028** — the highest mean and the tightest
spread in the table.

This resolves the worry left open at stage 2. The conserved hydrophobics at
Y3, L7, V11, L34 and L38 coincide with the interface positions, which raised the
possibility that the binder was burying its core against IL-6 rather than
against itself and would not fold alone. It folds alone, confidently and
reproducibly. The design is a real helical hairpin.

## 5. The n = 5 estimate was a fluke on the high side

The first run gave ipTM 0.581 ± 0.052 at n = 5. At n = 50 it is 0.466 ± 0.022,
so the original estimate was **high by 0.115**, about two standard errors, and
the true SD is 0.157 rather than 0.116. Insisting on the sampling before
reporting any difference was correct, and would have been correct even if the
means had agreed.

## What this does and does not establish

**Does:** in single-sequence mode with an unpinned target, Boltz interface ipTM
does not separate a designed binder from a random peptide, and its response to
target point mutations is unrelated to whether the mutated residue is at the
interface. Also, `short_9` folds on its own.

**Does not:** that co-folding models cannot filter binders, or that W157 does
not matter. One design, one model, one configuration — and a configuration we
have now identified as broken in a specific way.

## The next step is no longer optional

The diagnosis points at one cause: the model cannot hold the target still. Two
fixes, both already on the list, now become the critical path rather than
refinements.

1. **AF2 with the initial-guess protocol.** It hands the model the designed
   complex as its starting structure, so the target is never predicted from
   sequence and cannot drift. Its thresholds are calibrated against
   experimental hit rates. `params_model_1_ptm.npz` is already downloaded; what
   remains is JAX with CUDA in the container.
2. **Pin the target in Boltz**, by giving chain A an MSA or 1ALU as a template,
   and rerun this same eight-condition ladder unchanged.

Either way the test is the same and the prediction is written down now: **if the
gate passes once the target is pinned, the audit becomes runnable and the
earlier failure was configuration. If the gate still fails with the target
pinned, that is the stronger and more interesting result** — that the filter
the field relies on cannot tell a designed interface from an accidental one.
