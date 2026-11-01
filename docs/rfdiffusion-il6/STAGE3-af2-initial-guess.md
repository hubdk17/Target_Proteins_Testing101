# Stage 3, second instrument — AF2 initial guess

Environment installed and verified 2026-10-28 (job 174216). This file is the
method and the experiment design. Why AF2 at all, and what Boltz got wrong:
`LADDER-RESULT.md`. Project framing: `README.md`.

## Why this is the instrument the conclusion should rest on

Boltz's failure was diagnosed, not mysterious: in single-sequence mode it cannot
reliably fold IL-6 itself, interface confidence tracks whole-complex confidence
at r = +0.78, and every target point mutation moved the score whether or not it
touched the interface. A random 40-mer outscored the designed binder.

**The initial-guess protocol (Bennett et al. 2023) removes that failure by
construction.** The designed complex is handed to the model as its starting
structure, so the target is never predicted from sequence and cannot drift. It
is the same lever as pinning Boltz with a template, in the version the field has
calibrated against measured hit rates.

## The environment

`~/rfd/venv_af2`, built inside the NGC container. `jax 0.4.20` + `jaxlib
0.4.20+cuda12.cudnn89`, PyRosetta 2026.39, `tensorflow-cpu 2.21.0`, numpy held
at the container's 1.26.4. Repo `~/rfd/dl_binder_design`; it ships the
AlphaFold source, and `params_model_1_ptm.npz` (356 MB) is already in place.

Two pins are load-bearing and were found by reading a `--dry-run` list line by
line rather than by trusting that it resolved:

- **numpy at 1.26.4.** Unconstrained, pip installs numpy 2.4.6, and both jax
  0.4.20 and the 2022-era AlphaFold predate NumPy 2.
- **`nvidia-cudnn-cu12==8.9.7.29`.** jax 0.4.20 declares its CUDA dependencies
  with `>=` and no upper bound, so pip takes cuDNN 9 while the wheel is built
  for 8.9. That mismatch presents as a GPU problem, not a version problem.

**Not yet confirmed:** that jaxlib sees the GPU. The build job requested none.
`jax.devices()` must list a cuda device; `TFRT_CPU` means the driver was not
found.

## PyRosetta is not optional

It is the structure and IO layer of `af2_initial_guess/predict.py` on *both*
input paths, not only the silent-file one: `init(...)` at module level,
`pose_from_pdb` on the `-pdbdir` branch, `pose.split_by_chain()` to decide which
chain is the binder and how long it is, `input_check` over `pose.pdb_info()`,
and `af2_get_atom_positions(pose, ...)` building the AF2 features. Commenting
the imports leaves a script that cannot read its input.

## The input, and the decision behind it

`af2_in/short9_threaded.pdb` — the RFdiffusion complex
`outputs/il6_binder_site3_short_9.pdb` with the `short_9` ProteinMPNN sequence
threaded onto chain B.

Verified in that source file on 2026-10-28:

| | |
|---|---|
| chain A | 157 CAs, numbered 19..184, **all twenty residue types present** |
| chain A tryptophans | **exactly one, at residue 157** — independent confirmation of the numbering anchor |
| chain B | 40 CAs, numbered 1..40, **all glycine** (generated backbone) |

So only chain B needed threading, and renaming residues suffices because
RFdiffusion writes backbone-only coordinates and PyRosetta rebuilds missing
side-chain atoms with ideal geometry on load.

**The first attempt at this corrupted the file**, and the way it was caught is
worth carrying forward — see the second case in `CAUTION-origin-atom-bug.md`.
An awk lookup table built from twenty name pairs was looped to field 38 instead
of 39, so `Y → TYR` was never mapped, every tyrosine lookup returned an empty
string, and eight atom lines were written three characters short with their
coordinate columns shifted. `grep -c "GLY B"` returned 0 as predicted and a
`diff` on the ATOM count was clean, because both tested expected contents rather
than file integrity.

**The verification to use, for any fixed-column edit:**

```bash
awk 'substr($0,13,4)==" CA " && substr($0,22,1)=="B"{printf "%s",substr($0,18,3)}' af2_in/short9_threaded.pdb | wc -c
awk 'substr($0,1,4)=="ATOM"{print length($0)}' af2_in/short9_threaded.pdb | sort -u
```

The first must print **120** (40 residues × 3 characters). The second must print
**exactly one number** — one unique line length means no line was shifted. Plus
a runtime guard in the awk itself, `if(n==""){print "MAP FAIL at residue " r;
exit 1}`, so a missing mapping can never pass silently.

### Why not `dl_interface_design.py`

The repo's own ProteinMPNN wrapper would have been the obvious route, and it is
the wrong one here. It designs its own sequence at `-temperature 0.000001`
(effectively greedy argmax) with `-seqs_per_struct 1`, so it would produce a
*different* 40-mer from the `short_9` sample 2 that Boltz scored. Comparing AF2
on one sequence against Boltz on another would confound the filter comparison
with a design change, and the entire point of adding AF2 is to score the same
object with a better instrument.

### Why not a Boltz-predicted complex

Tempting, because Boltz already produced all-atom structures with the right
sequences. It would import Boltz's bias into the one comparison AF2 exists to
provide independently.

## Flags that matter

| flag | why |
|---|---|
| `-pdbdir` / `-outpdbdir` | directory input, no silent files needed |
| `-maintain_res_numbering` | **essential.** IL-6 has a real chain break at 51 → 61, so `-max_amide_dist 3.0` will judge the input bad and renumber by default. If residue 157 stops being 157, the mutation audit's bookkeeping is silently wrong |
| `-no_initial_guess` | the control Boltz never had: the same input scored with and without the trick, isolating what the initial guess contributes |
| `-force_monomer` | gives the `binderonly` and `targetonly` arms directly |
| `-recycle` | default 3 |
| `-scorefilename` | where ipTM and interface PAE are written |

## The experiment, once a GPU is free

1. **Smoke test.** One structure, default flags plus `-maintain_res_numbering`.
   Confirm it runs and that the output PDB still numbers chain A 19..184.
2. **Noise floor.** AF2 is deterministic given a fixed input and seed, unlike
   Boltz's diffusion sampling, so the replicate structure differs: spread comes
   from recycles and from input perturbation, not from resampling. Measure it
   before reporting any difference — the lesson from the n = 5 Boltz run that
   read 0.115 high.
3. **The same eight-condition ladder**, transferred unchanged, so the two
   instruments are compared on identical arms. The decisive cell remains
   `random40` versus `wt`: a filter that cannot put a designed binder above an
   accidental peptide cannot support the mutation audit either.
4. **Then, and only then, the W157A arm.**

## The pre-registered reading

- **The gate passes and `w157a` drops while `e93a` and `k120a` do not.** The
  project's result: a calibrated filter is sensitive to target interface
  chemistry, Boltz's failure was configuration, and the audit is runnable.
- **The gate passes and all three target mutations drop equally.** The
  sensitivity is to target sequence perturbation in general, not to this
  hotspot. Weaker, and only the controls distinguish it.
- **The gate fails here too.** The stronger and more interesting result: the
  filter the field actually relies on cannot tell a designed interface from an
  accidental one. That claim needs both instruments, which is why AF2 was worth
  the install.
