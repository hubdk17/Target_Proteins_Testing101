# BindCraft vs our pipeline — feasibility and the one comparison worth running

Written 2026-11-01. Sources: BindCraft README and wiki
(<https://github.com/martinpacesa/BindCraft>, wiki page "De novo binder design
with BindCraft"), read today. Everything below marked **[doc]** is from those
pages; **[ours]** is our own measured work; **[predicted]** is written down
before measuring, per our own rule.

---

## 1. What BindCraft actually is, in pipeline terms

Ours is three separate programs in a line: shape generator → sequence writer →
scorer. BindCraft is **one optimisation loop**: it backpropagates through the
AlphaFold2-multimer network and minimises a composite loss over the binder
sequence, in single-sequence (no-MSA) mode, annealing from a soft sequence
representation to a one-hot. One run of the loop is a "trajectory". **[doc]**

Then it post-processes: ProteinMPNN redesign of core and surface positions
(soluble weights by default) → final validation, where AF2 re-predicts the
complex, predicts the binder **alone** and compares it back to the designed
backbone (`Binder_RMSD`), and PyRosetta computes interface terms. Outputs:
`trajectory_stats.csv`, `mpnn_design_stats.csv`, `final_design_stats.csv`. **[doc]**

**The structural point: its objective function is the scorer.** It optimises the
binder directly against AF2's own confidence, and ranks finalists by `i_pTM`.
**[doc]** That is the same family of number our gate uses. So BindCraft does not
fix the weakness we have been documenting — it is the weakness, industrialised.

---

## 2. The free result: BindCraft's filter list is a map of which of our
## loopholes someone has already closed

This cost zero GPU hours and it changes our ranked list. BindCraft's documented
filter set includes `Binder_Energy_Score`, `dG`, `dSASA`, `dG/dSASA`,
`PackStat`, `ShapeComplementarity`, `Unrelaxed_Clashes`, `Relaxed_Clashes`,
`MPNN_score`, the AF2 confidence set (`pLDDT`, `pTM`, `i_pTM`, `pAE`, `i_pAE`,
`i_pLDDT`, `ss_pLDDT`), interface composition terms (H-bonds, unsatisfied
H-bonds, hydrophobicity, SASA fraction, secondary-structure fractions),
`HotspotRMSD` / `Target_RMSD` / `Binder_RMSD`, and binder-alone AF2 terms. **[doc]**

Against our four:

| Our loophole | Status in BindCraft |
|---|---|
| **1. No fold-back check** (15% of our sequences never adopt the designed shape, invisible until synthesis) **[ours]** | **CLOSED.** `Binder_RMSD` from a binder-alone prediction is a standard filter. **[doc]** This is our pipeline's gap, not the field's. |
| **2. The confidence gate is non-reproducible exactly where it decides** **[ours]** | **OPEN, and worse.** It mitigates overfitting by randomly switching among the five multimer models per optimisation step **[doc]**, but the final verdict is still one confidence number per design, with no replicates — and the design was *optimised against that number*. |
| **3. No interaction energy inside the pipeline** | **PARTLY closed.** PyRosetta `dG`, `PackStat`, `ShapeComplementarity` exist — but only as a *post-hoc filter* on a design already committed by a confidence-driven optimiser, never as a term the generator optimises. Also: thresholds ship as fixed numbers, and we have measured that the identical input file gives dG **+364.58 / −19.41 / −57.56** across three configurations **[ours]**, so a shipped threshold without its settings is not a measurement. |
| **4. Targeting channel is chemically empty** (one bit per residue) **[ours]** | **Not applicable / closed by construction.** There is no hotspot one-hot channel; hotspots enter as a contact term in the loss, evaluated through AF2's own representation of the target, which is chemistry-aware. Hotspots are also optional, and the docs admit the binder may ignore them if another site scores better. **[doc]** |

**Commercial reading: loophole 2 is the lead, not loophole 1.** A competitor
already ships the fold-back check. Nobody ships replicate-aware gating, and
BindCraft's design procedure makes the circularity sharper, not softer.

Also worth quoting back to anyone who claims otherwise: the BindCraft docs
themselves state **none of their computational metrics predict binding
affinity**, and that `i_pTM` is a binary binder/non-binder call, not an affinity
estimate. **[doc]** Our core claim survives intact and is now
competitor-corroborated.

---

## 3. What a run would cost here

- **Size.** IL-6 chain A is 157 residues; binder 40–55. So ~200–212 residues per
  trajectory. The docs give ~5 min for a 250-residue trajectory on an H100, with
  an A100 ~4× slower, and memory the real limit (32 GB ≈ 550 residues). **[doc]**
  We are comfortably inside memory on either node.
- **Volume.** Docs: expect a few hundred trajectories for accepted designs on an
  easy target, 1,000–10,000 on a hard one. **[doc]**
- **So 300 trajectories at IL-6 size ≈ 25 GPU-h on the H100**, ~100 GPU-h on an
  A100. As an array throttled to `%4` on `dgxhnode6` that is one overnight
  submission — the same order as one of our existing arms. **Affordable.**
- Storage: 5.3 GB of AF2 params **[doc]** plus trajectory output.

**The install is the risk, not the compute.** The shipped installer is
`bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'conda'` **[doc]**, and we
have **no conda on this cluster**. Three things to settle before committing:

1. Can the login node reach the AF2 multimer parameter host at all (compute nodes
   are host-filtered)?
2. Does our existing AF2 container's JAX satisfy ColabDesign, or does it need its
   own venv? Our AF2 env is pinned at JAX 0.4.20 / dm-haiku 0.0.13.
3. PyRosetta is already working in that container, which removes the hardest
   dependency.

Cap this at ~2 rounds of environment debugging, as usual.

---

## 4. The comparison to run — and the one that would be rigged

**Do not run a success-rate bake-off.** We have no binding measurements, so
"which tool makes better binders" is unanswerable and any number we produce
would be decoration.

**Also do not compare our raw pool against BindCraft's accepted designs.** Their
accepted set is filtered *on* fold-back and *on* confidence. Comparing a filtered
output to our unfiltered one would manufacture a win for them out of nothing but
pipeline stage. **Compare at matched stages:** their `mpnn_design_stats.csv`
(post-MPNN, pre-final-filter) against our `scN150` pool.

**The informative experiment.** Both pools, both scorers, one table:

|  | our AF2 initial-guess `pae_interaction` | fold-back `binder_aligned_rmsd` | PyRosetta interface dG, one fixed configuration |
|---|---|---|---|
| our pool | measured | measured | pending |
| BindCraft pool | to run | to run | to run |

Everything goes through **our** scripts, not theirs, so the numbers are
commensurable: `43_build_complex.py` to assemble complexes, our pinned-node AF2
scoring, `90_interface_energy.py` for dG. Target: IL-6 site III, hotspots
`A47,A155,A157,A158`, binder length band 40–55 — identical to our arms.

### Predictions, recorded before measuring **[predicted]**

1. **BindCraft designs beat ours on the confidence gate by a wide margin** —
   median `pae_interaction` ~4–8 against our pool's ~10–18, with a far higher
   pass rate under the ~10 cut. This is *not* evidence they are better binders;
   it is what optimising against the metric looks like.
2. **Fold-back, measured at the matched pre-filter stage, is roughly a tie** —
   their pre-filter pool within about 10 points of our 85% under 2 Å, and
   plausibly *worse*, because an AF2-hallucinated sequence has no independent
   backbone it was written for. If prediction 1 holds and 2 holds, that is the
   result: **optimising against the gate buys gate score and does not buy
   foldability.**
3. **PyRosetta dG will not separate the two pools cleanly** at any single
   configuration, consistent with our own dG settings-sensitivity.
4. **Their advantage on our gate shrinks under replicates.** Both pools scored
   ~3× on a pinned node; given our measured heavy-tailed noise, a meaningful
   fraction of their sub-10 designs sit within the ~3-unit resolvability band.

If prediction 1 fails — if BindCraft does *not* dominate our confidence gate —
that is interesting in the opposite direction and means the two confidence
scorers are less interchangeable than assumed.

### Honest statement of what this settles

It settles **whether a design optimised against a confidence score is thereby
more foldable**. It does not settle which pipeline produces binders. Any claim
beyond the first sentence is not supported by the data we can collect.

---

## 5. Phasing

- **Phase 0 — done, free.** Section 2. Changes the loophole ranking today.
- **Phase 1 — three cluster checks** (host reach for params, JAX/ColabDesign fit,
  disk). One paste-and-go block. Stop here if the env fight exceeds ~2 rounds.
- **Phase 2 — 100 trajectories** at IL-6 site III, throttled `%4`,
  `--exclude=dgxanode03`, distinct output name. Read `mpnn_design_stats.csv`.
- **Phase 3 — the 2×2 table** through our own scripts, with replicates on a
  pinned node.

Phase 2 and 3 only if Phase 1 is clean. Phase 0 already pays for the thread.
