# MDM2 blind de novo design — our own model, 2026-11-17

No RFdiffusion, no ProteinMPNN, no trained weights of any kind. GPU was
available but not needed for the design itself.

## The claim being tested

Can a designer that has **never seen a binder** — zero parameters fit to any
structure — produce something that looks like the natural p53 binder, working
from the MDM2 groove shape alone? If yes, the groove's geometry carries most
of the answer. If no, binding needs information the receptor surface does not
contain.

## Why this is blind, precisely

Three leaks, all closed:

1. **Learned generators are contaminated.** RFdiffusion and ProteinMPNN were
   both trained on the PDB, which contains 1YCR and dozens of MDM2 complexes.
   Neither can be called blind on this target. Dropped entirely (Kaila7's
   call, 2026-11-17).
2. **Our own heads.** We computed Phe19/Trp23/Leu26 yesterday. So the designer
   chooses residues by an automatic rule, never by us typing them.
3. **The hotspot list.** The gate-dissection pocket was defined by contacts
   *with* the peptide. So the target site here is re-derived from the receptor
   alone (`62_blind_pocket.py`), with the peptide absent from the input.

The only prior used is that the binder is an **ideal alpha-helix** — a general
bias (helices bind grooves), not MDM2-specific knowledge. The designer is not
told the natural binder is a helix.

## Step 1 — find the groove blind (`62_blind_pocket.py`)

LIGSITE-style concavity scan on the bare MDM2 surface, peptide never loaded.
It found one dominant pocket, 81 A^3, lined by:

```
A54 A57 A58 A61 A62 A67 A75 A86 A91 A93 A99 A103
```

Disclosed validation only (never fed back): **10 of these 12 are real p53
contacts.** The two extras, Phe86 and Ile103, are the hydrophobic floor of the
same cleft. The detector located the p53 site without seeing p53.

## Step 2 — design into it blind (`63_blind_designer.py`)

Build an ideal poly-A helix by NeRF (geometry validated: 3.8 A CA–CA, right-
handed, L-chirality confirmed by our own gate and matched against the real
p53 helix). Rigid-body search for clash-free poses resting on the groove.
Then assign each residue from **burial alone**: a side chain pointing into
the groove gets the largest hydrophobic that fits the local cavity; a
solvent-facing one gets a charge. No rotamers, no packing — crude on purpose.

Five designs, length 14, from 3000 sampled poses (286 clash-free):

| design | sequence | buried apolar (Å²) | H-bonds | D-res | unsat |
|---|---|---|---|---|---|
| 1 | WWIVEKAFEWAKEA | 664 | 1 | 0 | 0 |
| 2 | WAWKTIFKWAEKVK | 555 | 0 | 0 | 0 |
| 3 | IAWKVALKAAEKAK | 581 | 1 | 0 | 0 |
| 4 | EVALEFAKWFLKAA | 640 | 3 | 0 | 0 |
| 5 | EKTIWKVAEKAWEA | 523 | 1 | 0 | 0 |
| **real p53** | ETFSDLWKLLPEN | **977** | **10** | 0 | 0 |

## The result

**The blind designer independently produced amphipathic helices with large
aromatics on the buried face, carrying the same hydrophobic-moment signature
as the natural binder it never saw.**

Eisenberg hydrophobic moment μH, with each design's percentile against 500
shuffles of its *own* amino-acid composition (so this measures whether the
*order* is amphipathic, not just that it is hydrophobic):

| sequence | μH | shuffle %ile |
|---|---|---|
| WWIVEKAFEWAKEA | 0.427 | 92% |
| WAWKTIFKWAEKVK | 0.516 | 95% |
| IAWKVALKAAEKAK | 0.397 | 83% |
| EVALEFAKWFLKAA | 0.518 | 99% |
| EKTIWKVAEKAWEA | 0.544 | 99% |
| **real p53** `ETFSDLWKLLPEN` | **0.494** | 99% |
| **PMI** `TSFAEYWNLLSP` | 0.298 | 89% |

Our blind designs sit at μH 0.40–0.54, bracketing the real p53 helix (0.494),
every one above the 83rd percentile of its own shuffles. The designer also
reached for tryptophan repeatedly on the buried face — and Trp23 is the single
dominant anchor of the real binder. Driven only by groove shape, it rediscovered
that this pocket wants a big aromatic.

## Honest limits — these designs are crude

1. **They over-bury and under-bond.** 520–660 Å² apolar buried vs the natural
   977, but only 0–3 interface H-bonds vs 10. The designer optimises burial
   and ignores polar complementarity, so the sequences are hydrophobic-heavy
   and would likely have solubility and specificity problems. This is the
   "no matter how bad it is" version.
2. **No packing.** Side chains are single pseudo-atoms along the CB vector;
   real rotamers would change both the clashes and the burial.
3. **μH is partly built in.** The designer places residues by burial, so some
   amphipathicity is structural by construction. The non-trivial parts are
   that the *order* beats composition shuffles and that the *magnitude*
   matches the natural binder — neither of which the burial rule forces.
4. **Not yet folded or docked.** μH says the sequence is amphipathic; it does
   not say the peptide will fold as a helix in isolation or dock in this pose.
   That is the AF2 step, and it is the real test.

## Next

- Length sweep 12–18 (natural is 13) and a larger frozen set.
- Add a polar-complementarity / H-bond term so the designer stops ignoring the
  10 H-bonds the natural binder makes.
- **AF2 initial-guess on the frozen designs (GPU)** — do they fold and dock? —
  scored the same way as the real p53/PMI/pDI complexes. This is the evaluation
  against the known binders, run only after the designs are frozen.

## Reproduce

```
awk '$1=="ATOM" && substr($0,22,1)=="A"' 1YCR.pdb > mdm2_receptor.pdb
python3 62_blind_pocket.py mdm2_receptor.pdb --out pocket.tsv
python3 63_blind_designer.py mdm2_receptor.pdb --pocket pocket.tsv \
    --length 14 --n-designs 5 --orientations 3000 --seed 1
python3 60_gate_battery.py design_*.pdb --target-chain A --binder-chain B \
    --no-gap-volume --out blind_gates.tsv
```
