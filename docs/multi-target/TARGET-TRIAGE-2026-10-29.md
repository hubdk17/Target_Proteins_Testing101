# Triage of Diptiman Sir's 19 targets for the binder pipeline

Date: 2026-10-29. List received 20 Oct 2026. Judged against the same criteria used
to pick IL-6 (1ALU) on 2026-10-21.

## What our pipeline actually needs from a target

RFdiffusion -> ProteinMPNN -> Boltz/AF2 can only be pointed at a target that is:

1. **One soluble chain, or one domain that folds on its own.** A heterodimer or an
   assembly means the binder has to be diffused against the right chain and the
   interface we care about may be buried.
2. **No essential cofactor, lipid, metal or nucleic acid.** RFdiffusion sees only
   protein atoms, so anything held together by a ligand comes in as a different
   shape than it really has.
3. **Not heavily glycosylated.** The sugars are absent from the crystal file, so
   we would design binders onto a surface that is covered in vivo.
4. **Extracellular, if we care about the result being usable.** A binder for an
   intracellular target cannot get into a cell on its own.
5. **A crystal structure with low B-factors and no large missing loops.**
6. **A mapped epitope**, so hotspots can be chosen rather than guessed.
7. **Modest size**, roughly under ~300 residues, so the GPU cost per design stays
   near IL-6's ~1 min.

IL-6 met all seven. That is the bar below.

---

## Tier 1 — can be set up now

| Target | PDB | Size | Why it works |
|---|---|---|---|
| **PD-1** | **3RRQ** (apo), 4ZQK (with PD-L1) | ~115 aa IgV domain | The best fit on the list, and arguably a better learning target than IL-6: single small immunoglobulin domain, secreted/ectodomain, the PD-L1 interface is mapped residue-by-residue, and **de novo binders against it are published**, so it doubles as the positive control our AF2 scoring is missing. |
| **PCSK9** | **2P4E** (full), 3BPS (with EGF-A of LDLR) | ~690 aa, 3 domains | Soluble and secreted, the therapeutic epitope (the EGF-A groove on the catalytic domain) is defined by a co-crystal, and clinically validated. Harder than PD-1 only because it is bigger and multi-domain — we would diffuse against the catalytic domain alone. |

## Tier 1.5 — one unusual case worth calling out

| Target | PDB | Note |
|---|---|---|
| **BCL2** | **2XA0**, 4LVT | Intracellular and membrane-tethered, so a protein binder has no delivery route — it fails criterion 4. But structurally it is the single most favourable surface on the list: the BH3 groove **natively binds a single alpha helix**, which is exactly what RFdiffusion produces most readily. If the goal is to learn what the pipeline can do, this is the easiest win on the list. If the goal is a usable binder, skip it. Needs your call on which. |

## Tier 2 — possible, but only one domain of them, and all intracellular

These are all real folded domains with good structures. Every one fails criterion
4 (inside the cell), and each needs a decision about which domain to use before any
job can be written.

| Target | PDB | What has to be decided first |
|---|---|---|
| STAT3 | 1BG1 | Use the SH2 domain alone; the full-length protein is a multi-domain dimer on DNA. |
| Casp3 | 2J30, 1NME | It is a p17/p12 heterotetramer with a catalytic cysteine; which chain do we target? |
| pTEN | 1D5R | Folded and single-chain, but the C-terminal tail is flexible and would need trimming. |
| MAPK | 1ERK (ERK2), 1P38 (p38a) | **"MAPK" is a family, not a protein.** Needs Diptiman Sir to name ERK, p38 or JNK. |
| JAK1 / JAK2 | 3EYG (JAK1 JH1), 4IVA (JAK2 JH1) | Kinase domain only, and the ATP site is a small-molecule pocket, not a binder epitope. |
| Beta-RAF | 4MNE, 1UWH | Kinase domain only; it dimerises and is disordered without an inhibitor bound. |
| BRCA1 | 1T29 (BRCT), 1JM7 (RING) | **Most of BRCA1 is intrinsically disordered.** Only the BRCT tandem domain is usable, and the RING domain only exists as a complex with BARD1. |
| GLP1R | 3IOL (extracellular domain only) | The whole receptor is a 7-helix membrane protein and is out. The isolated extracellular domain in 3IOL **is** a small soluble target and is genuinely interesting. |
| INSR | 4ZXB (ectodomain), 3W11 (kinase) | Already ruled out on 2026-10-21 as a large glycosylated multi-domain disulfide-linked dimer. An isolated domain could work later. |

## Tier 3 — not suitable, with the reason

| Target | Why not |
|---|---|
| **TLR-1** (2Z7X) | Ruled out 2026-10-21. Large glycosylated leucine-rich-repeat ectodomain, and it only functions as a TLR1/TLR2 heterodimer around a bound lipopeptide. |
| **TLR-4** (3FXI) | Same, plus it needs the MD-2 co-receptor and LPS present to have its real shape. |
| **AKT** (3O96) | Ruled out 2026-10-21. Intracellular, and the PH and kinase domains move relative to each other. |
| **PI3K-alpha** (4OVU) | ~1050-residue p110a/p85a heterodimer, intracellular. Fails size and chain count badly. |
| **Microtubule** (1JFF for the tubulin dimer) | A polymer of ab-tubulin heterodimers requiring GTP. The thing itself has no single structure file; cost per design would be enormous. |
| **Histone** (1KX5, 7XD1) | Comes as an octamer wrapped in DNA. The surface is strongly basic and its real binding partner is a nucleic acid, which our pipeline cannot represent at all. |
| **Fo-F1 ATPase** (6FKF, 5ARA) | A multi-subunit rotary membrane machine. Out on every criterion. |

Score: **2 ready, 1 judgement call, 9 need a domain chosen, 7 out.**

---

## Recommendation: do not queue 19

As of today (job 174676) the pipeline has produced its first calibrated interface
number, and our best IL-6 design `short_9` **failed** it: pae_interaction 26
against a published cut-off near 10. Two independent co-folding models now decline
to endorse it.

So the pipeline is not yet known to produce a binder for the one target we have
spent eight days on. Running it against 19 targets now most likely buys 19
unvalidated failures and a very large GPU bill, and we would not be able to tell a
bad design from a broken pipeline.

**Add exactly one target: PD-1 (3RRQ).** The reason is not throughput, it is
diagnosis. PD-1 has published, experimentally confirmed de novo binders. If we
run one of those through our own AF2 scoring and it scores near 10, our scoring is
sound and `short_9` is simply a bad design. If the published binder also scores 26,
the fault is in our scoring setup, and scaling to 19 targets would have told us
nothing. **That single run is worth more than the other eighteen combined.**

Once PD-1 is calibrated, PCSK9 is the natural second, and it is the one on this
list with a real clinical story to tell Diptiman Sir.

## What to tell Diptiman Sir

The list is strong on biological importance and thin on structural tractability:
most of it is intracellular, and five entries are assemblies rather than proteins.
Worth asking him two concrete questions:

1. Which MAPK does he mean — ERK, p38 or JNK?
2. Is the intent a usable binder, or a computational demonstration? That one answer
   moves nine of these targets between tiers, because intracellular targets are fine
   for the second and not for the first.
