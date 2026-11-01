# Hotspots make binders foldable, not just well-placed

Jobs 172779 and 172780, analysed 2026-10-28. 300 designs, the largest sample in
the project. Project index: `README.md`. The underpowered run this corrects, and
the straight-helix diagnosis: `RESULTS.md`.

## The result

Identical scripts, one flag changed, binder length fixed at 50, no radius-of-
gyration potential in either arm. Compactness gate at Rg ratio ≤ 1.30 against
`2.2 · N^0.38`.

| arm | pass | rate |
|---|---|---|
| hotspots `A47,A155,A157,A158` | **47 / 150** | **31%** |
| no hotspots | **13 / 150** | **9%** |

Two-proportion test: pooled p = 0.20, SE = 0.046, **z = 4.9, p ≈ 1 × 10⁻⁶**.

Distribution of Rg ratios across all 300:

| ratio | count |
|---|---|
| 1.0 | 4 |
| 1.2 | 56 |
| 1.3 | 16 |
| 1.4 | 6 |
| 1.5 | 5 |
| 1.6 | 5 |
| 1.7 | 2 |
| 1.8 | 6 |
| 1.9 | 1 |
| 2.0 | 2 |
| 2.1 | 21 |
| **2.2** | **176** |

## What was predicted, and what was wrong

The project rule is to write the expected numbers down first.

| prediction | actual | verdict |
|---|---|---|
| 300 designs measured | 300 | right |
| both arms pass at ~20% | 31% and 9% | **wrong** |
| hotspots do not affect compactness | 3.6-fold difference, p ≈ 10⁻⁶ | **wrong, and wrong about the mechanism** |
| a large spike at ratio 2.2 | 176 of 300 | right |

The n = 30 run gave 6/30 against 5/30 and was read as "indistinguishable, as
expected, since hotspots should not affect compactness". That reading was a
statement about power, not about biology, and it was reported as though it were
the latter.

## The mechanism, stated as a claim that can fail

**Hotspot conditioning confines the contact region to a small patch, so a
50-residue chain has to turn back on itself to reach it.** That produces
hairpins. Without hotspots the model may drape anywhere on the target's surface,
and for a chain of fixed length the lowest-resistance answer is a single long
helix lying across it.

If that is right, hotspots are doing double duty: steering *where* the binder
binds and, as a side effect, making it *foldable*. The flag is not advertised
that way.

### The test, and it held

Cross-tabulating arm against Rg class, with straight defined as ratio ≥ 2.05:

| arm | pass (≤1.30) | middle | straight (≥2.05) |
|---|---|---|---|
| hotspots | 47 | 27 | **76 (51%)** |
| no hotspots | 13 | 15 | **122 (81%)** |

The rod rate falls from 81% to 51%, **z = 5.6**, in exactly the predicted
direction. The confounder was excluded in the same pass: **every design in both
arms has n = 50**, so binder length cannot account for any of it.

### The competing explanation, still open

Site III may simply be a curved patch that happens to suit hairpins, which would
make this a fact about this surface rather than about confinement in general.
**The cheap test is to run the same comparison with hotspots on a flatter site**,
site I or site II. Until that is done the claim stays bounded to one target, one
binder length and one patch.

## Why this matters beyond stage 1

**Two thirds of everything stage 1 produces at this length is a rod**, 197 of
300 at ratio 2.1 or above. That is the strongest available argument for enabling
the ROG guiding potential during diffusion rather than filtering for compactness
afterwards. Brute force recovered the yield; it did not fix the generator.

**Hotspots help but do not fix it.** Half the hotspot arm is still rods, so this
finding changes how the flag should be described, not whether stage 1 needs
repairing.

**And 60 designs now pass the gate**, against the single `short_9` that every
stage 2 and stage 3 result in this project rests on. The weakest thing about the
work so far has been n = 1; that is now fixable.

## The method lesson worth keeping

This rerun was submitted to recover sample size after a shortcut cut n = 30 to
n = 6. It was not designed to test anything about hotspots. It produced the
project's clearest effect anyway.

**When a comparison looks null, check its power before believing it.** The
discipline of writing the predicted number down first is well established here;
this adds the matching discipline for the *reading* of a result, not just the
measurement of one.
