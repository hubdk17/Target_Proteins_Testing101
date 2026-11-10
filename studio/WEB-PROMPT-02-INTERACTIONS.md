# Add the interactions tab: precise interface analysis, and a reference library of measured dynamics

Follow-up to `WEB-PROMPT.md`. Hand this to the coding agent after week 1 passes its walkthrough.

## Context

A wet-lab user wants to *see* the interaction, not read a table of scores. This task adds a
top-level **Interactions** tab with two halves, and the split between them is the whole point of
the task.

## The one constraint that shapes everything here

**We have never run molecular dynamics on a designed binder.** The pipeline produces a *single
static predicted complex* per design. There is no trajectory, no time axis, no ensemble.

So the tab has two halves and they must never be confused:

| Half | Data | What it is about |
|---|---|---|
| **A. This design's interface** | computed from the predicted complex | the user's own design — static, precise, real |
| **B. Reference dynamics** | real MD from the sibling project, already on disk | *other*, known complexes — a measured contrast, never a result about the design |

**Anything animated, interpolated, morphed or labelled "simulation" on a design's own page would
be invented data.** Do not generate a trajectory client-side. Do not interpolate between
predicted states. Do not animate a static structure and call it dynamics. A plausible-looking
animation of a complex nobody simulated is the single worst thing this app could ship, because a
bench scientist would reasonably believe it.

What makes half A worth building is **precision**, not motion: exact contacts, exact distances,
exact coverage, every number traceable.

---

## Step 1. Routes

```
src/app/interactions/page.tsx                    the library index: reference systems + what they show
src/app/interactions/[systemId]/page.tsx         one reference system, its MD views
src/app/runs/[runId]/designs/[design]/interface/page.tsx   half A, on the design itself
```

Add **Interactions** to the sidebar as a top-level entry, sibling to Runs. On a design detail
page add an **Interface** sub-tab.

## Step 2. What to copy from the sibling app

Everything in half B already exists at `D:\Desktop\PPI_INTERACTION\web\`, built against real
trajectories. Copy components and data; do not re-derive anything from raw trajectories — the
full set is ~367 GB and the app has never needed it, because the importer scripts already reduce
it to JSON.

```
src/components/contacts/*      contact-map, contact-charts, contact-detail,
                               contacts-table, contacts-explorer
src/components/hbond/*         hbond-table, hbond-detail, hbond-explorer,
                               lifetime-scatter, switching-table, classes.ts
src/components/water/*         bridge-network, hydration-viewer, interface-plane,
                               site-detail, water-tables, consensus-panel,
                               crystal-overlap, encoding
src/components/diagnostics/salt-bridge-bars.tsx
src/components/replicas/replica-scatter.tsx
src/lib/contacts.ts  src/lib/series.ts  src/lib/chains.ts
```

And the reduced data, which carries its own provenance:

```
src/data/1brs/hbond_dynamics.json
src/data/4hfz/{hbond_dynamics,contacts,hotspots,effect_sizes,conditions,diagnostics}.json
src/data/4zqk/{hbond_dynamics,structure,system}.json
src/data/{trypsin,sod1,p53tet}/campaign.json   src/data/trypsin/mutation_study.json
src/data/results/*   src/data/runs/*
public/structures/*.pdb
```

These views were written for a project whose chain columns are deliberately neutral
(`chain_a_*` / `chain_b_*`, never a protein name), so they port without renaming. Keep that.

## Step 3. Half A — this design's interface

Everything here is computed from the design's own predicted complex. All of it is cheap, exact
and needs no GPU.

1. **The complex, in Mol\***. Reuse `molstar-viewer.tsx` exactly as copied, including its
   CommonJS import. Colour binder and target from `palette.ts`. **Take the chain split from the
   API's `binder_length`, never from chain labels in the file** — one of our pipeline's outputs
   renumbers both chains from 1 and loses the labels.
2. **Interface contact table.** Every binder residue within the grip cutoff of a target residue,
   with the measured distance and the pair of residue identities. Sortable by residue, by
   distance, by partner. Clicking a row selects that pair in the viewer.
3. **Contact map**, binder residues × target residues, from the copied `contact-map.tsx`. One
   static map, no time axis.
4. **Grip and footprint per residue** — contacts within 8 Å and 10 Å, which
   `42_grip_rank.py` already computes. Render as a per-residue bar strip aligned under the
   sequence.
5. **Hotspot coverage.** For each hotspot the user confirmed on the review screen: does this
   design actually touch it, and at what distance. This is the one panel that answers "did it aim
   where I asked", so put it first.
6. **Buried surface area**, and the chirality result as a per-residue strip so a flagged residue
   is locatable in the structure.
7. **An interface energy, only if the API supplies one, and never as the headline.** If shown it
   must print its full settings next to it. On one unchanged input, three configurations of our
   interface-energy tool gave **+364.58, −19.41 and −57.56** in the same units. An energy without
   its settings is not a measurement. And a negative value is not evidence of binding.

Each panel carries the run chip from week 1 and nothing is sortable by confidence.

## Step 4. Half B — the reference library

Frame it as what it is: **measured interfaces, shown so the user can see what this pipeline does
not measure.** That framing is honest and it is also the most useful thing on the page for
someone from the bench.

Index page: one card per reference system (barnase–barstar, p53–MDM2, PD-1/PD-L1,
trypsin–BPTI, SOD1, p53 tetramer) with its PDB id, what it is, how long was simulated, and which
views exist. Follow the sibling app's convention of listing a view that has not been produced as
**disabled rather than absent**, so the gap is visible in the UI.

Per system, the copied views: contact occupancy over time, H-bond lifetime distributions and
switching, water bridges and hydration sites, salt-bridge occupancy, per-replica agreement.

Every page in half B carries a persistent banner, not a footnote:

> Measured molecular dynamics of a known complex. **Not your design, and not a prediction about
> it.** Shown so you can compare what a simulated interface reveals against what the design
> pipeline measures.

Two things from that project belong on these pages because they are the reason the library is
worth showing at all, and both are measured results:

- **A contact count can miss a change that destroys the interaction.** In one reference system a
  change that moved the key residue 4.3 Å and swapped its entire coordination shell moved total
  interface contacts from **283.2 to 282.6**, against a noise floor of 2.84. Put this on the
  contact pages: it is why our pipeline refuses to use a contact-count tiebreaker.
- **An implicit-solvent binding energy can be precisely wrong rather than noisy.** On a complex
  with a dissociation constant around 10⁻¹⁴ M it returned **+2.98 kcal/mol** — the wrong sign —
  while reproducing across seeds at 5.6%. Put this beside any energy view.

## Step 5. The MD stage as roadmap, disabled

On a design's Interface page, show a disabled **Dynamics** section with one sentence: molecular
dynamics of designed complexes is not part of the pipeline today, and when it is, it will test
whether the complex stays assembled — not how tightly it binds.

This follows the sibling app's own pattern of showing an unproduced analysis as disabled so the
roadmap is visible in the product rather than missing from it. Do not implement it. Do not leave
a button that looks clickable.

## Step 6. Forbidden

- any animation, morph, interpolation or "trajectory" of a design that was not simulated
- the word *simulation* anywhere in half A
- a half-B view rendered on a design page without its banner
- an energy number without its settings
- any claim, in words or by layout, that an interface view predicts binding
- a sort-by-confidence control, here as everywhere

## Step 7. Contracts

Extend `src/lib/types.ts`. The API computes all of this; the UI renders it.

```typescript
export interface InterfaceContact {
  binder_residue: { number: number; label: string };
  target_residue: { number: number; label: string };
  distance: number;                 // Angstrom, CA-CA
  within_grip: boolean;             // <= 8 A
  within_footprint: boolean;        // <= 10 A
}

export interface HotspotCoverage {
  hotspot: string;                  // as confirmed on the review screen
  label: string;
  touched: boolean;
  nearest_distance: number | null;
}

export interface DesignInterface {
  run_id: string; design: string;
  binder_length: number;            // the chain split; never read labels from the file
  contacts: InterfaceContact[];
  grip_per_residue: number[];
  footprint_per_residue: number[];
  hotspot_coverage: HotspotCoverage[];
  buried_surface_area: number | null;
  chirality: Judgement;
  interface_energy: {
    value: number; units: string;
    settings: Record<string, unknown>;   // required; render it, never hide it
  } | null;
  structure_url: string;
  thresholds_version: string; config_hash: string;
}

export interface ReferenceSystem {
  id: string; pdb_id: string; name: string; description: string;
  simulated_ns: number; n_replicas: number;
  views: { id: string; label: string; available: boolean }[];
}
```

Add a new API route `GET /api/runs/:id/designs/:design/interface -> DesignInterface`, mocked like
the rest. Reference systems are static data, so they may be read directly from `src/data/` —
this is the one documented exception to "nothing fetches around the data seam", and the README
must say so and say why.

## Step 8. Verification

```bash
npm run lint && npm run typecheck && npm run build && npm run walk
```

Extend the walkthrough:

- [ ] a design's Interface page shows hotspot coverage first, and a design that missed a hotspot says so
- [ ] clicking a contact row selects that residue pair in the viewer
- [ ] the chain split comes from `binder_length`; relabelling chains in the test PDB changes nothing
- [ ] no half-A page contains the word "simulation", and nothing on it animates
- [ ] every half-B page carries its banner above the fold
- [ ] the Dynamics section is disabled and not clickable
- [ ] an interface energy without `settings` is not rendered at all
- [ ] no sort-by-confidence control exists, in source or in any rendered page
- [ ] a reference view marked `available: false` is shown disabled, not hidden

Note for whoever writes that check: the sibling project's agent found that server-rendered output
is one long line with `<!-- -->` between text nodes, so an unanchored `/a.*b/` matches across the
whole document and a naive "no confidence sort" grep flags the caveat text that explains why
confidence cannot be sorted on. Anchor the pattern to the control's markup.
