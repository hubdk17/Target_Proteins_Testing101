# Binder design pipeline — the procedure, in the corrected order

Written 2026-10-29 from the IL-6 (1ALU) run, for reuse on any target.
This is the *method*; the multi-target thread already has the target prep
script (`02_prep_targets.py`), the RFdiffusion batch (`31_rfd_multi.sbatch`)
and the batch plan in `docs/multi-target/`. This document adds the
two things that run **after** those, including the new gate found on 2026-10-29.

The one idea behind the whole thing: **gates throw designs out; rankings sort
the survivors.** A design must pass every gate before any ranking is applied.
We previously ranked before gating and carried a broken design through the whole
project. Do not repeat that.

---

## The five stages and where each gate sits

```
Stage 0  Prepare target           ─ not a gate, but get it right or nothing after is valid
Stage 1  RFdiffusion backbones  ──▶ GATE 1: compactness   (Rg ratio < ~1.30)
Stage 2  ProteinMPNN sequences  ──▶ GATE 2: self-consistency (RMSD < 2 Å)   ← NEW, the one we were missing
         ── rank survivors by grip (site coverage), this is a RANKING not a gate ──
Stage 3  AF2 complex, initial guess ──▶ GATE 3: interface  (pae_interaction < ~10)
```

Costs on an A100, per design, so you can budget: RFdiffusion ~1 min, ProteinMPNN
~5 s on GPU, self-consistency AF2 monomer ~70 s, complex AF2 ~2 min. The monomer
gate is a fifth the cost of the complex, which is exactly why it goes first.

---

## Stage 0 — prepare the target (CPU, login node)

Every target needs this done by hand. It cannot be scripted blindly because a
numbering gap and a real chain break look identical until you measure.

1. Pick ONE PDB entry. Prefer a co-crystal with the natural partner so the
   epitope is defined. Note resolution.
2. Strip to protein-only (drop waters, ligands, ions):
   `grep '^ATOM' raw.pdb > clean.pdb` is the crude version; better to keep only
   the one chain you want.
3. **Distinguish a numbering gap from a real chain break.** Walk consecutive CA
   atoms and measure the distance; a true peptide bond is ~3.8 Å. A jump in
   residue *number* with the CAs still ~3.8 Å apart is just numbering convention
   (keep the contig whole). A CA-CA distance above ~5 Å is a real break (split
   the contig there). IL-6 had one real break (51→61) hidden among none; trypsin
   had six numbering gaps and zero breaks. **Measure, do not assume.**
4. Write the contig string from that, e.g. IL-6 was `A19-51/A61-184`.
5. Choose hotspot residues from the known epitope and **verify each against the
   chain sequence** — read the residue name at that number out of the PDB and
   confirm it is the amino acid you expect. We verified W157 was TRP twice.

---

## Stage 1 — RFdiffusion backbones + GATE 1 (compactness)

Run with the radius-of-gyration potential ON. It is load-bearing, not a
refinement: without it a fixed-length binder collapses to a single straight
helix ~75% of the time.

```
potentials.guiding_potentials=["type:binder_ROG,weight:3"]
potentials.guide_scale=2
potentials.guide_decay="quadratic"
```

**GATE 1, compactness.** For each backbone, compute radius of gyration of its CA
atoms and divide by the compact expectation:

  ρ = Rg / (2.2 · N^0.38)        keep ρ < ~1.30

N is the binder residue count. The 2.2·N^0.38 is the Rg a well-packed globular
protein of that length would have, so ρ is length-free: ρ≈1 is compact, ρ≈2.2 is
a straight rod. One-liner (chain B is the binder):

```
awk 'substr($0,1,4)=="ATOM" && substr($0,13,4)==" CA " && substr($0,22,1)=="B"{
  x=substr($0,31,8)+0;y=substr($0,39,8)+0;z=substr($0,47,8)+0;
  sx+=x;sy+=y;sz+=z;q+=x*x+y*y+z*z;n++}
  END{rg=sqrt(q/n-(sx/n)^2-(sy/n)^2-(sz/n)^2);
      printf "N=%d Rg=%.2f rho=%.2f\n", n, rg, rg/(2.2*n^0.38)}' design.pdb
```

Discard anything with ρ above 1.30 before spending another second on it.

---

## Stage 2 — ProteinMPNN sequences

Run MPNN on each surviving backbone with chain A (target) fixed and chain B
(binder) designed. Generate **many** sequences per backbone, not a handful — the
whole point of Gate 2 is to pick among them, so give it a real choice. We used
8 and it was enough to expose the problem but is on the low side; 8–48 is sensible.
MPNN writes binder-only FASTA, one sequence per `>T=` header.

---

## GATE 2 — self-consistency (this is the new one, 2026-10-29)

Thread each MPNN sequence back onto its own RFdiffusion backbone, predict that
chain **alone** with AF2 (`-force_monomer`), superpose the prediction on the
design, and keep only sequences whose Cα RMSD to the design is **under 2 Å**.

RMSD = √( (1/N) Σ dᵢ² ) after best rigid superposition (Kabsch). It asks one
question: does this sequence actually fold into the backbone it was written for?

### 2a. Thread the sequences (PyRosetta or the checked Python below — NOT raw awk)

Threading rewrites the residue name in columns 18–20 of each backbone ATOM line.
Fixed-column edits have silently corrupted data twice on this project, so the
script asserts on every lookup and every line length, and you compile it before
running so a mangled paste fails loudly:

```
cat > thread.py <<'PY'
import glob, os
M = dict(zip("ACDEFGHIKLMNPQRSTVWY", "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()))
os.makedirs("sc_in", exist_ok=True)
for fa in sorted(glob.glob("mpnn_out/seqs/*.fa")):
    tag = os.path.basename(fa)[:-3]
    bb  = [l for l in open("outputs/%s.pdb" % tag) if l[:4]=="ATOM" and l[21]=="B"]
    res = sorted({int(l[22:26]) for l in bb})
    L   = open(fa).read().split("\n")
    seqs= [L[i+1].strip() for i,h in enumerate(L) if h.startswith(">T=")]
    for k, s in enumerate(seqs, 1):
        assert len(s)==len(res), (tag,k,len(s),len(res))
        m = dict(zip(res, s)); out=[]
        for l in bb:
            nl = l[:17] + M[m[int(l[22:26])]] + l[20:]
            assert len(nl)==len(l)
            out.append(nl)
        open("sc_in/%s_s%d.pdb" % (tag,k), "w").write("".join(out))
    print(tag, "nres", len(res), len(seqs), "threaded")
PY
python3 -m py_compile thread.py && echo INTACT || echo "MANGLED - STOP"
python3 thread.py
# integrity checks that can actually fail:
ls sc_in | wc -l                                   # expect #backbones x #seqs
awk '{print length($0)}' sc_in/*.pdb | sort -u     # expect a single value
awk '{print substr($0,18,3)}' sc_in/<one>.pdb | sort -u | tr '\n' ' '  # no blanks; TYR present
```

### 2b. Predict each threaded chain alone and read the RMSD

Use the `dl_binder_design` AF2 initial-guess environment (see the env notes in
project memory). Drop `-debug` for a batch so one bad input cannot kill the rest;
use a **fresh** `-scorefilename` every run (the score file is append-only).

```
python predict.py -pdbdir sc_in -outpdbdir sc_out \
  -scorefilename sc.sc -maintain_res_numbering -force_monomer
```

Then, per backbone, the best (lowest) `binder_aligned_rmsd` across its sequences
is the deciding number:

```
grep "^SCORE:" sc.sc | awk 'NR>1||$2+0==$2{n=$11; sub(/_s[0-9]+_af2pred$/,"",n);
  if(!(n in b)||$2<b[n])b[n]=$2} END{for(k in b) printf "%-32s best=%.2f\n",k,b[k]}' | sort
```

**Keep sequences under 2 Å. Discard the rest.** If *every* sequence for a
backbone fails, the backbone is not designable — drop the backbone, do not try to
rescue it.

### THE TRAP, do not fall in it

**Do not rank sequences by pLDDT (confidence).** On IL-6 the two sequences that
folded into *rods* had the HIGHEST confidence (93–95) while the six correct ones
sat at 75–79. A simple wrong fold is easy to predict and so scores confidently.
Gate on RMSD-to-design, never on confidence.

---

## Rank the survivors by grip (a ranking, not a gate)

Only now, among sequences that passed Gate 2, rank by interface quality: CA-CA
contacts to the target within 8 Å (grip) vs 10 Å (footprint), and coverage of
the intended epitope. This is where the biology comes back in. Carry the top few
forward.

---

## Stage 3 — AF2 complex, initial guess + GATE 3 (interface)

Build the complex input: binder chain FIRST, then target, with residue numbers
**unique across both chains** (offset the binder, e.g. 201+, keep the target's
crystal numbers). Both requirements are enforced by the script and each silently
corrupts or rejects the run if wrong — see the AF2 input notes in memory.

```
python predict.py -pdbdir cplx_in -outpdbdir cplx_out \
  -scorefilename cplx.sc -maintain_res_numbering
```

**GATE 3, interface.** `pae_interaction` is AF2's cross-chain error in Å; the
published bar for a plausible binder is **< ~10**. Treat it as a gate at the
failing end — it saturates near 30, so it sorts good from bad but does not finely
rank among bad ones. Also run `-no_initial_guess` as a free control: if the two
agree, the guess is doing nothing and the number means less.

IL-6 `short_9` sample 2 scored 26 here. But that was the rod — it never should
have reached Stage 3, because it failed Gate 2. That is the whole lesson in one
number.

---

## Checklist to hand each target

- [ ] Stage 0 done: one PDB, cleaned, breaks-vs-gaps measured, contig, hotspots verified against sequence
- [ ] Stage 1 run with ROG potential; Gate 1 applied (ρ < 1.30)
- [ ] Stage 2 many MPNN sequences per surviving backbone
- [ ] **Gate 2 self-consistency; keep RMSD < 2 Å; drop undesignable backbones; never rank by pLDDT**
- [ ] Rank survivors by grip
- [ ] Stage 3 complex, binder-first + unique numbering; Gate 3 pae_interaction < 10; run the no-guess control
- [ ] Fresh job name AND fresh score filename every run; verify success by output-file count and elapsed time, never exit code
