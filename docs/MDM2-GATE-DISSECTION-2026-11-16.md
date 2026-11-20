# MDM2 gate dissection, 2026-11-16

CPU only, no GPU used. Target: MDM2 N-terminal p53-binding domain.

## Headline

**Buried apolar surface area — a term that exists nowhere in our pipeline,
and costs two seconds of CPU per structure — recovers the experimentally
known hotspot order on MDM2. Contact counting, which our gates do use, does
not.**

Experimental hotspot order for the p53 peptide against MDM2 is
Trp23 >> Phe19 > Leu26, with Asp21 and Lys24 irrelevant.

| arm | Δ buried apolar (Å²) | significance | Δ heavy contacts | Δ gap volume (Å³) |
|---|---|---|---|---|
| W23A | **−113.2** | 21 σ | −35 | −66.5 |
| F19A | **−65.1** | 12 σ | −35 | +11.3 |
| L26A | **−33.0** | 6 σ | −9 | +12.0 |
| L22A | −17.9 | 3.4 σ | −13 | +25.7 |
| E28A | −5.7 | 1.1 σ | −4 | −45.6 |
| L25A | −2.8 | 0.5 σ | 0 | −4.5 |
| T18A | −2.8 | 0.5 σ | 0 | −4.8 |
| **K24A** (control) | −2.8 | 0.5 σ | 0 | −46.3 |
| **D21A** (control) | **0.0** | 0 σ | 0 | 0.0 |
| E17A | +5.7 | — | −18 | −82.0 |
| N29A | +8.5 | — | −15 | −70.7 |
| S20A | +15.8 | — | −2 | −25.4 |

Buried apolar area gives **W23 > F19 > L26 > L22**, all above 3 σ, in the
right order, with both negative controls at the noise floor and D21A at
exactly zero.

Contact counting **ties W23A with F19A at −35**, cannot separate the two
most important residues in the complex, and ranks the non-hotspot **N29A
(−15) above the real hotspot L26A (−9)**.

Gap volume as implemented is **worse than useless here**: its three largest
responses are E17A, N29A and K24A — two rim residues and a negative
control — and it moves the *wrong way* for F19A, L26A and L22A. Do not
ship it as a gate without reworking it.

## Significance is measured, not assumed

SASA is rotation invariant, so the spread over six random rigid rotations
of the same complex is pure sphere-point discretisation error. At 128
points, `dsasa_apolar` has **sd = 5.3 Å²**. Everything above is quoted
against that. Anything under ~16 Å² (3 σ) is not resolved — L22A at 3.4 σ
is real but weak, and L25A, T18A, K24A are nulls, not small effects.

## What the structure actually looks like

Chain A is MDM2 **25–109**, 85 residues, **no chain breaks**. Chain B is
p53 **17–29**, `ETFSDLWKLLPEN`. Every side chain is complete. All 94
chiral residues are L: the N–C–CA–CB improper is **−122.3 ± 3.0°**, which
is now the calibrated reference for this test.

Pocket, any MDM2 heavy atom within 4.5 Å of the peptide — 22 residues:

```
A25,A26,A50,A51,A54,A57,A58,A61,A62,A67,A70,A71,A72,A73,A75,A91,A93,A94,A96,A99,A100,A104
```

Leu54, Ile61, Met62, Tyr67, Val75, Val93, Ile99, Tyr100 all confirmed.
**Gln72 is the third-largest contributor on the MDM2 side at 35 contacts**,
which the usual description of this pocket underplays.

Interface size: dSASA total **1464 Å²**, so ~730 Å² buried per side, 67 %
apolar. That agrees with the published interface area for this complex.

The canonical Trp23 NE1 → Leu54 backbone O hydrogen bond is present at
**2.83 Å**.

**Asp21 and Lys24 make exactly zero contacts with MDM2.** They point at
solvent. That hands us free negative-control positions inside every batch.

## Two gates that cannot work on this ladder, by construction

- **Gate 1, compactness** is shape only. Every arm of a fixed-backbone
  alanine scan has the same backbone, so the Rg ratio is identical in all
  13 arms. It cannot register a 1000-fold affinity difference.
- **Gate 2, self-consistency** is target-blind for the same reason.

That is worth stating as a finding rather than an aside: **two of our
three gates are structurally incapable of responding to the thing we are
asking them about.**

## A gate that is actively misleading on an alanine scan

Helix propensity. Alanine is the most helix-favouring residue on the Pace
& Scholtz scale, so **every** alanine mutant scores as *easier* to
pre-organise than wild type. On this ladder the gate is anti-correlated
with affinity by construction.

It is still worth keeping, but only for a **sequence-diverse** comparison —
p53 17–28 versus PMI versus pDI — where the peptides differ in intrinsic
helicity. The underlying point stands: **nothing anywhere in our pipeline
prices the conformational entropy a flexible peptide pays to fold before
it binds**, and for a 12-mer that is a large fraction of the measured Kd.
That is a second missing-capacity claim alongside the target-chemistry one.

## A zero that is a measurement, not a dead check

`unsat_buried_polar` returns 0 on all 13 arms. Thirteen identical values
is exactly the invariant-across-samples signature that cost us three rounds
of conclusions on the phantom-atom bug, so it was tested: deleting the two
MDM2 polar atoms nearest Trp23 NE1 makes the gate return 1 and name
`B/TRP23:NE1`. **The gate fires when it should.** Zero unsatisfied buried
polars in a natural, evolved complex is the right answer.

## Honest limits

1. **Fixed backbone, no relaxation.** The cavity left by a deleted side
   chain never closes, so the loss from a large buried residue is
   systematically over-read. This is standard computational alanine
   scanning, not a shortcut we invented, but it means the absolute numbers
   are not free energies.
2. **`dumb_dG_apolar` is −24.4 kcal/mol for wild type**, against a real
   ΔG of roughly −8.5 kcal/mol. Buried area with no desolvation penalty
   overestimates by about 3×. It is a **ranker, not an energy**, and must
   never be quoted as one.
3. **Hydrogen bonds are an upper bound.** With no hydrogens present there
   is no donor–H–acceptor angle, so a close-but-badly-oriented pair is
   indistinguishable from a real bond.
4. **No measured affinities in this table yet.** Everything above is
   ranking against *published hotspot identity*, not against Kd. The
   quantitative test still needs the SKEMPI 2.0 1YCR rows, which this
   container cannot pull.
5. The 1YCR coordinates came from the RFdiffusion repository's example
   inputs, not from RCSB, because RCSB is unreachable from the build
   container. The file was verified (chains, numbering, sequence,
   completeness, chirality) rather than trusted, but an authoritative copy
   should replace it.

## What this does and does not say about MDM2

It is a **positive control**. 1YCR is from 1996 and is in the training set
of every tool we use; the bound solution is one helix with three
hydrophobic anchors. None of the above is evidence that a model has
learned anything. What it *is* evidence for is that a two-second classical
term outperforms our contact-based gate at recovering known energetics —
and that is a statement about our gates, which is what we were testing.

MDM2 is also intracellular with no delivery route. Method target, never a
program target.

## Next

1. SKEMPI 2.0 1YCR rows → the quantitative Kd correlation.
2. PyRosetta REF15 interface energy on the same 13 arms, in the same
   batch, for a head-to-head against buried apolar area.
3. Rework or drop gap volume.
4. Sequence-diverse ladder (p53 / PMI / pDI) to make the helix-propensity
   gate testable.
5. AF2 `pae_interaction` on the same arms when a GPU is free — the only
   remaining gate we have not pitted against this.

## Reproducing

```
python3 61_make_ala_scan.py 1YCR.pdb --chain B --outdir arms/
python3 60_gate_battery.py arms/*.pdb --target-chain A --binder-chain B \
    --out gates_scan.tsv
```

13 arms, about 25 s total on one core.
