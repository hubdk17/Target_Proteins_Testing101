# Stage-1 loader: atom records are matched by name, not by position

Verified 2026-11-01 against the upstream RFdiffusion source. The loader is
`rfdiffusion/inference/utils.py`, function `parse_pdb_lines` (line 526); there is
no `rfdiffusion/parsers.py` in the current tree, so any note citing that path is
citing a file that does not exist.

## What the code does

```
526  def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
529      for l in lines:
530          if l[:4] == "ATOM" and l[12:16].strip() == "CA":      # residue list = CA records only
541      xyz = np.full((len(res), 14, 3), np.nan, dtype=np.float32)
542      for l in lines:
543          if l[:4] != "ATOM": continue
554          for i_atm, tgtatm in enumerate(util.aa2long[util.aa2num[aa]][:14]):
558              if tgtatm is not None and tgtatm.strip() == atom.strip():
561                  xyz[idx, i_atm, :] = [...coordinates...]
562                  break
565      mask = np.logical_not(np.isnan(xyz[..., 0]))
566      xyz[np.isnan(xyz[..., 0])] = 0.0
```

`aa` is columns 18-20 of the record: the residue label. `aa2long`
(`rfdiffusion/chemical.py:30`) holds one canonical atom-name tuple per label,
and only the first 14 entries (heavy atoms) are consulted. An atom is written
into the tensor **only if its name appears in the tuple belonging to the label
on its own record**. There is no positional fallback, no warning, and no counter:
a record whose name is absent from that tuple falls out of the loop, its slot
stays NaN, `mask` stays False at line 565, and the coordinate is zeroed at 566.

The ALA tuple's first 14 entries are `N CA C O CB` and nine `None`
(`chemical.py:30`). So under the label ALA, exactly five heavy atoms can ever be
read, whatever else the file contains.

## Consequence 1 — the "label X, atoms of Y" arm is not buildable by file edit

An arm consisting of a record set labelled ALA while physically carrying a full
ring loads as a plain stripped ALA. It is byte-for-byte a different input file
and tensor-for-tensor the same input as the existing `ctl`/`mut` arms. Running
it would produce a fourth set of designs that is a replicate, not a new
condition, and nothing in the logs would say so.

**This is a limitation to document, not to patch.** Treat that cell of the design
as unreachable at stage 1 without a change to the loader.

## Consequence 2 — renaming a label without editing atoms gives a chimera

Because matching is per-name, renaming a label keeps every atom whose name is
shared between the old and new tuple and silently drops the rest. TRP -> TYR
keeps `CG CD1 CD2 CE2`, drops `NE1 CE3 CZ2 CZ3 CH2`, and never supplies `CE1 CZ
OH`: a partial, non-physical side chain, not a substitution. `50_build_mutation_arms.py`
already avoids this by stripping to ALA's five atoms in every non-wt arm; the
rule is that a label change must be paired with an explicit atom set, never
applied on its own.

## Consequence 3 — two smaller edges of the same behaviour

- A residue with no `ATOM ... CA` record never enters `pdb_idx` (line 530) and is
  absent from the loaded structure entirely, side chain and all.
- A non-standard label reaches `util.aa2num[aa]` at line 554 without the
  `if ... in keys() else 20` guard used for `seq` at line 535, so it raises
  `KeyError` at load. The failure is loud; it is the only one here that is.
- `HETATM` records are skipped by the heavy-atom loop (line 543) regardless of
  name.

## Check it on the cluster copy

```
grep -n "aa2long\|aa2num\|== \"CA\"\|tgtatm" ~/rfd/RFdiffusion/rfdiffusion/inference/utils.py | head -20
```
