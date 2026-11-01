# Caution case study: the phantom atom at the origin

**Date:** 2026-10-22
**Cost:** three rounds of wrong conclusions, one unnecessary GPU job, and several
false results written into the project records before it was caught.
**Root cause:** one uninitialised variable in a hand-written awk parser.

---

## What the bug was

The contact-measuring scripts collected coordinates like this:

```awk
if (c=="B") { bx[nb]=x; by[nb]=y; bz[nb]=z; nb++ }
```

`nb` was never initialised. **In awk, an unset variable used as an array
subscript is the empty string, not zero.** So the first atom was stored at
`bx[""]`, and only then did `nb++` make the counter 1. The second atom went to
index 1, the third to 2, and so on.

**Index 0 was never filled.**

The loop that measured distances started at `j=0`:

```awk
for (j=0; j<nb; j++) { d=(ax[i]-bx[j])^2 + ... }
```

`bx[0]`, `by[0]`, `bz[0]` were unset, and awk evaluates unset variables as 0 in
numeric context. So every measurement silently included **an invisible extra atom
at the coordinate origin (0, 0, 0)**, and counted any target residue near the
origin as being in contact with the binder.

The same bug in the target array (`na`) created an unset entry whose residue
number read as `0` — the phantom "residue 0" that appeared in every contact list.

## How it was proved

Residue 94 of IL-6 sits at (1.163, −2.587, 3.933). Its distance to the origin:

```
sqrt(1.163² + 2.587² + 3.933²) = sqrt(23.51) = 4.85 Å
```

The script reported the "nearest binder residue" as being 4.85 Å away. Exact
match. The thing it was measuring to was the origin, not the binder.

## Why it survived three rounds of analysis

This is the part worth internalising.

1. **The numbers were plausible.** 58 contacts, 44 contacts, 21 contacts — all
   in the range a real protein interface occupies. Nothing looked broken.
2. **The artifact was stable, so it read as signal.** RFdiffusion centres its
   output near the origin, which is inside the target protein. The phantom atom
   therefore "contacted" the protein's core residues — including helix B's
   90–98 — in every single structure. A consistent result looks like a finding.
3. **It fed a narrative.** The false contacts appeared on a helix with no
   hotspot, which fit the story "the binder has an unwanted second foothold".
   Two further experiments were designed to chase that non-existent foothold.
4. **A small anomaly was dismissed twice.** A residue numbered `0` kept appearing
   in output lists. It was waved away as "a quirk of the script" and even
   filtered out cosmetically (`if (ar[i] > 0)`) rather than explained. That
   filter hid the symptom while leaving the disease.

## What finally caught it

**A quantity that should have varied across independent samples was identical.**

Ten stochastically generated binders, all different lengths, all in visibly
different positions (centroids spanning 20 Å), reported contacts with *the same
seven residues at the same distances to 0.1 Å*. That is physically impossible,
and it was the first observation that could not be explained away.

The decisive diagnostic was printing the **raw matched line** rather than a
derived number. `nearest line: []` — empty — proved the lookup was returning
nothing, while still producing a confident-looking distance.

## Rules adopted from this

1. **Initialise counters explicitly.** `BEGIN { na=0; nb=0 }`. In awk an unset
   subscript is `""`, not `0`.
2. **Test that a measurement varies when its input varies.** If ten different
   structures give the same number, the number is not measuring the structures.
   This is a cheap and powerful check, and should be run *before* interpreting
   any batch result.
3. **Never filter an anomaly you cannot explain.** A stray value is data about
   your code. Suppressing it cosmetically is how a bug survives.
4. **Print the raw record, not just the derived number.** Deriving 4.85 Å from
   nothing looks identical to deriving it from something.
5. **Validate a hand-written parser against a known answer** before trusting it —
   for example, measure a distance the structure file already implies.
6. **Prefer an established library** (BioPython, MDAnalysis) to a hand-rolled
   column parser when the result is going to carry scientific weight. The awk was
   convenient for one-off looks; it should not have become the measurement
   instrument for three consecutive conclusions.

## What had to be retracted

- "The binder has a second foothold on helix B" — never existed.
- "Two footholds 27 Å apart" — artifact.
- "Helix B survives the tighter cutoff, so it is a real grip" — it was the
  protein's own core, near the centroid.
- Run 3 (ten designs, 40–55 residues) was designed specifically to eliminate the
  phantom foothold. The job itself was cheap and the designs are usable, but the
  reasoning that motivated it was wrong.
- All contact totals in runs 1 and 2 were inflated, run 1's by more than half.

---

# Second case, 2026-10-28 — the unmapped tyrosine

Same class of bug, three weeks of practice later, caught in one round instead of
three. Recording it because the thing that caught it was a deliberate habit, and
the thing that *failed* to catch it was two checks that looked reasonable.

## What happened

Threading the `short_9` sequence onto the RFdiffusion backbone used an awk
lookup table built from pairs:

```awk
split("A ALA C CYS D ASP E GLU F PHE G GLY H HIS I ILE K LYS L LEU M MET N ASN P PRO Q GLN R ARG S SER T THR V VAL W TRP Y TYR", a, " ")
for (i=1; i<=38; i+=2) m[a[i]] = a[i+1]
```

Twenty pairs is forty fields, so the loop bound must be 39. At 38 it stops after
`W → TRP` and **never maps `Y → TYR`**. Every tyrosine lookup returned the empty
string, so those lines were written as columns 1–17 followed by columns 21
onward: three characters short, with every coordinate column shifted left by
three. Eight atom lines silently corrupted, in a file about to be handed to
AlphaFold.

## The two checks that passed for the wrong reason

- **`diff` on the ATOM count was clean.** The line *count* never changed. Only
  the line *contents* did.
- **`grep -c "GLY B"` returned 0, as predicted.** The 40-mer contains no
  glycine, so this was supposed to prove every position had been rewritten. But
  an empty residue name is not GLY either, so the check could not distinguish
  "all forty renamed correctly" from "thirty-eight renamed and two blanked".

Both checks tested an *expectation about the contents*. Neither tested the
file's *integrity*.

## What caught it

The written-down prediction. The expected residue list began `SER SER TYR ASP
GLU TYR LEU …` and the actual output began `SER SER ASP GLU LEU …`, with no TYR
anywhere and 38 entries instead of 40. Writing the first eleven residues down
before running the command is what made a two-residue absence visible
immediately rather than at the point where AlphaFold returned a strange number.

## The check that should have been written first

```bash
awk 'substr($0,1,4)=="ATOM"{print length($0)}' file.pdb | sort -u
```

**One unique line length means no line was shifted.** More than one is the
corruption itself, not a proxy for it. For fixed-column formats like PDB, this
is the cheapest possible integrity test and it is indifferent to what the
contents were supposed to be.

Paired with a guard that makes the failure loud rather than silent:

```awk
n = m[substr(seq,r,1)]
if (n == "") { print "MAP FAIL at residue " r; exit 1 }
```

## The lessons, added to the list

- **A lookup table built from parallel fields has an off-by-one waiting in its
  loop bound.** Count the pairs, then count the fields, then write the bound.
  Better still, let a guard catch the miss at runtime.
- **Prefer a check on structural integrity over a check on expected contents.**
  A content check encodes what you think the answer is; an integrity check
  catches the failure you did not think of.
- **A clean `diff` on a count proves almost nothing** when the transformation
  rewrites lines in place.
- **An empty string is a valid-looking value.** Any lookup that can miss should
  fail loudly, because a blank slots into a fixed-column format without
  complaint and shifts everything after it.
