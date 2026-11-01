# Stage 2 audit — is ProteinMPNN sensitive to a target point mutation?

Designed 2026-10-28, not yet run. This is the one unaudited stage of the
pipeline and the cheapest remaining experiment in the project: about a minute of
GPU time. Project framing: `README.md`. The stage 3 result that makes this
necessary: `LADDER-RESULT.md`.

## Why this stage is the interesting one

The three-stage pipeline treats the target's side-chain chemistry differently at
each step:

| stage | sees W157 vs W157A? |
|---|---|
| RFdiffusion | **No.** Conditions on backbone coordinates and hotspot residue *indices*. A surface point mutation that does not move the backbone is literally the same input. |
| ProteinMPNN | **Yes, in principle.** The fixed chain's sequence is supplied as one-hot input and conditioned on, so the network is told that position is a tryptophan. |
| Boltz / co-folding | Measured, and it fails: it moves for any target mutation, interface or not, because it is reading target foldability rather than interface chemistry. |

So ProteinMPNN is the only stage that both *can* see the mutation and is
actually used in the standard pipeline. If it does not respond either, the
project's claim strengthens from "the scoring metric is confounded" to "no stage
of the pipeline asks the binding question".

## Prediction, written down before running

The project rule is to state the expected number first (see
`CAUTION-origin-atom-bug.md` for why).

**Expect a small effect, confined to a few positions.** ProteinMPNN is dominated
by backbone geometry, and the backbone here is byte-identical between the two
runs. Removing a tryptophan's bulk opens a cavity, so the packing-aware response
would be to favour *larger* hydrophobics at the binder positions facing 157.

Quantitatively: a shift in amino-acid preference at **2 to 4** binder positions,
the ones within about 8 Å of residue 157, and no detectable change elsewhere.
Mean MPNN score should move by less than 0.02.

**If nothing changes at all, that is the stronger result**, not a null. It would
mean the one stage with access to the information does not use it.

## The experiment

Two runs on the **identical** `short_9` complex backbone, differing only in the
identity of chain A residue 157. Chain A fixed, chain B designable, 48 sequences
each, same seed.

### 1. Locate the complex PDB

```bash
cd ~/rfd
ls -la outputs/ | grep -i short
grep -c "^ATOM" outputs/*short_9*.pdb
```

The file is `outputs/il6_binder_site3_short_9.pdb`, **confirmed 2026-10-28**:
chain A is IL-6 with 157 CAs spanning 19..184, chain B the binder backbone with
40 CAs spanning 1..40. RFdiffusion preserved the crystal numbering, so residue
157 in this file really is Trp157 and no index arithmetic is needed.

### 2. Build the W157A copy

One line. It keeps only N, CA, C, O and CB of chain A residue 157 and renames
the residue, which is exactly what an alanine is. Everything else passes
through untouched, so the two files differ in nothing else.

```bash
awk '{ if (substr($0,1,4)=="ATOM" && substr($0,22,1)=="A" && substr($0,23,4)+0==157) { an=substr($0,13,4); gsub(/ /,"",an); if (an=="N"||an=="CA"||an=="C"||an=="O"||an=="CB") print substr($0,1,17) "ALA" substr($0,21) } else print }' outputs/il6_binder_site3_short_9.pdb > mpnn_audit_w157a.pdb
```

### 3. Verify before trusting it

**First establish which case we are in.** RFdiffusion writes backbone-only
output, so chain A may carry only N, CA, C, O per residue rather than full side
chains. If residue 157 has 14 atoms it is a full tryptophan and the mutation
removes 9 of them; if it has 4, there is no side chain to remove and the edit is
a rename alone — which is still the right input for ProteinMPNN, since MPNN
reads residue identity rather than side-chain coordinates. Either case is fine;
silently assuming the wrong one is not.

```bash
awk 'substr($0,1,4)=="ATOM" && substr($0,22,1)=="A" && substr($0,23,4)+0==157' outputs/il6_binder_site3_short_9.pdb | wc -l
awk 'substr($0,1,4)=="ATOM" && substr($0,22,1)=="A" && substr($0,23,4)+0==157' mpnn_audit_w157a.pdb
diff <(grep -c "^ATOM" outputs/il6_binder_site3_short_9.pdb) <(grep -c "^ATOM" mpnn_audit_w157a.pdb)
```

The second command must print only `ALA` at A 157 — four lines if the input is
backbone-only, five if it had a CB. If it prints `TRP` anywhere, the column
offsets are wrong for this file and nothing downstream is valid.

### 4. Run both

Use `sbatch` with `--gres=gpu:1`, never interactive `srun` — see `METHODS.md`.
ProteinMPNN is about 11 minutes per structure on CPU and seconds on a GPU.
Copy `40_mpnn.sbatch` and change only the input path and output directory, so
the two runs are identical in every other respect.

### 5. The comparison

Per-position amino-acid frequency across the 48 sequences, for each of the 40
binder positions, in both runs. The positions to look at first are the binder
residues whose CA is within 8 Å of IL-6 residue 157 — get that list with the
site III contact counter from `METHODS.md`, modified to report per-residue
distances to 157 rather than to the whole hotspot set.

Report three things: the set of positions whose modal residue changed, the
magnitude of the largest frequency shift, and the mean MPNN score in each run.

## What this cannot say

It cannot say whether the resulting sequences bind better or worse. It measures
whether the network's output distribution responds to the information at all,
which is a question about the pipeline rather than about IL-6.
