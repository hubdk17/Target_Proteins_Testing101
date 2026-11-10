# Build the de novo studio web app, week 1: the whole UI against a mock backend

Hand this file to the coding agent in VS Code. It is written to be executed in order.

## Context

We are building a hosted service that lets a bench scientist run a de novo binder design
pipeline (RFdiffusion backbones → ProteinMPNN sequences → AlphaFold2 scoring) without a
terminal. The pipeline itself already exists as numbered scripts in this repository and runs on
an HPC cluster today. **There is no GPU access for one week**, so this task is the entire front
end plus a mock backend: every screen, every state, every failure mode, demoable end to end with
no GPU and no cloud account.

Two things already exist and must be reused rather than reinvented.

### 1. The verdict layer: `studio/`

> **Correction, 2026-11-09.** This section originally said "in this repository" and told you to
> port `studio/review.py`. On a machine with no clone of the scripts repository, that had no
> referent, and the first agent to run this brief correctly reported the package as missing.
> `studio/` is real but it lives on the branch `claude/project-thread-b81vad` of
> `github.com/hubdk17/Target_Proteins_Testing101` (pull request #3, not yet merged) — it is **not
> on `main`** and **not on the local disk**. To read it:
>
> ```bash
> git clone -b claude/project-thread-b81vad \
>   https://github.com/hubdk17/Target_Proteins_Testing101.git pipeline
> ```
>
> You do not need it to build the front end: write the contracts to the descriptions below and
> keep `src/lib/mock/` swappable, which is what Step 9 asks for anyway. Read it when the real
> backend is wired, and say in the README that the mock stands in for it.

A Python package (see `studio/README.md`) that owns every threshold and every verdict. **The web
app must never compute a verdict, a threshold or a pass rate.** It renders what the API returns.
Read these before writing any UI:

| File | What the UI needs from it |
|---|---|
| `studio/thresholds.toml` | the only place a threshold lives; every verdict carries its hash |
| `studio/gatekeeper.py` | the verdict vocabulary and the refusals |
| `studio/store.py` | the run/design/judgement/event/spend tables the API exposes |
| `studio/provenance.py` | the record every result page must link to |
| `studio/review.py` | the stage-0 review screen, already written as a CLI — port its logic to the UI verbatim |

The verdict vocabulary is fixed and closed: `pass`, `fail`, `borderline`, `refused`,
`saturated`. `pass` and `borderline` carry forward; `refused` and `saturated` do not.

### 2. A sibling Next.js app to seed from: `D:\Desktop\PPI_INTERACTION\web\`

A finished, well-structured Next.js 16 app for a different project (molecular dynamics results).
It is **statically exported with no server**, so it is not the app we are building — but its
design system, shell, chart helpers and Mol\* viewer are exactly what we need, and they encode
details that took real time to get right. Read its `README.md` and `src/lib/data.ts` first.

Its stack, which we match exactly so the copied components compile unchanged:

```
next 16.3.5 (App Router)    react 19.2.8    typescript 5
tailwindcss 4 (@tailwindcss/postcss)        shadcn/ui, style "radix-nova", baseColor neutral
radix-ui 1.6.7   lucide-react   framer-motion 13   molstar 5.11.0
d3-array, d3-color, d3-interpolate, d3-scale, d3-shape
path alias: @/* -> ./src/*
```

**Read `node_modules/next/dist/docs/` before writing routing or data-fetching code.** That app's
`AGENTS.md` warns that Next 16 differs from training data, and it is right.

---

## Step 1. Create the app

Create a new Next.js app at `web/` **in this repository** (the one holding `studio/`). A new app,
not a fork of the PPI one, for three reasons: we need a server (that app sets
`output: "export"`), we need authentication and a database, and the two products will diverge.

```bash
# from the repository root
npx create-next-app@16.3.5 web --typescript --tailwind --eslint --app --src-dir --use-npm
cd web
npx shadcn@4 init      # style: radix-nova, baseColor: neutral, cssVariables: true
```

Then install the same runtime dependencies as the sibling app, pinned to the same versions:
`radix-ui@^1.6.7 lucide-react framer-motion@^13 molstar@^5.11.0 class-variance-authority clsx
tailwind-merge tw-animate-css d3-array d3-color d3-interpolate d3-scale d3-shape`, plus their
`@types/d3-*` as dev dependencies.

`next.config.ts` must **not** set `output: "export"`. Keep `trailingSlash: false`. This app has a
server.

## Step 2. Copy the reusable layer

Copy these from `D:\Desktop\PPI_INTERACTION\web\` and change nothing except import paths:

```
src/components/ui/*              badge, button, card, scroll-area, separator, sheet, switch,
                                 table, tabs, toggle, toggle-group, tooltip
src/components/shell/panel.tsx   the panel primitive
src/components/shell/controls.tsx
src/components/charts/*          export-button, pointer-tooltip, sort-header
src/components/motion/reveal.tsx
src/components/molecule/*        molstar-viewer.tsx, renderer.ts, molecule-canvas.tsx
src/lib/palette.ts               colour tokens, already validated for colour-vision deficiency
src/lib/export.ts                SVG / PNG / CSV export
src/lib/format.ts  src/lib/utils.ts  src/lib/use-element-size.ts  src/lib/geometry.ts
src/lib/residues.ts  src/lib/stats.ts
src/app/globals.css              mirrors palette.ts; keep them in sync
components.json  postcss.config.mjs  eslint.config.mjs
```

**Preserve one detail exactly as it is in `molstar-viewer.tsx`: Mol\* is loaded on demand from
its CommonJS build.** The sibling app's README records why — the ESM build's circular imports
break under production bundling. Do not "modernise" this import.

Adapt `src/lib/export.ts` so every exported file carries the run id, the thresholds version and
the config hash instead of that project's force-field caveat. Keep the principle: **no file
leaves this app without its provenance in it.**

## Step 3. The data seam

The sibling app's single best decision is that `src/lib/data.ts` is **the one place data enters**
and no component imports JSON directly. Keep that, and make it API-backed.

Create `src/lib/api.ts` as the only module that talks to the backend, and
`src/lib/types.ts` holding the contracts below. Every component reads through these. **No
component may construct a URL or fetch directly.**

```typescript
export type Verdict = "pass" | "fail" | "borderline" | "refused" | "saturated";

export const CARRIES_FORWARD: Record<Verdict, boolean> = {
  pass: true, borderline: true, fail: false, refused: false, saturated: false,
};

export interface Judgement {
  gate: "gate1" | "gate2" | "gate3" | "chirality";
  design: string;
  verdict: Verdict;
  value: number | null;
  threshold: number | null;
  reason: string;                      // written by the backend; render it verbatim
  detail: Record<string, unknown>;     // n_replicates, sd, sem, lo, hi, tight_floor, ...
  thresholds_version: string;
  config_hash: string;
}

export interface Rate {                // a pass rate is never a bare percentage
  passed: number; total: number;
  low: number; high: number;           // Wilson interval
  underpowered: boolean;
}

export type RunState =
  | "created" | "awaiting_review" | "queued" | "running"
  | "paused_cost" | "finished" | "failed" | "stopped" | "cancelled";

export type EventKind =
  | "queued" | "started" | "progress" | "finished"
  | "hard_failure" | "silent_failure" | "stall"
  | "refusal" | "ceiling_reached" | "cancelled";

export interface RunEvent {
  at: number; stage: string | null; kind: EventKind;
  message: string; detail: Record<string, unknown>;
}

export interface Run {
  run_id: string; target: string; state: RunState;
  created_at: number;
  thresholds_version: string; config_hash: string;
  cost_ceiling_usd: number; spent_usd: number;
  spec: RunSpec;
}

export interface RunSpec {             // every field here was confirmed by a human
  pdb_source: { kind: "rcsb"; id: string } | { kind: "upload"; filename: string };
  target_chain: string;
  partner_chain: string | null;
  contig: string;
  hotspots: string[];
  binder_length: [number, number];
  n_designs: number;
  n_sequences_per_backbone: number;    // >= 8, enforced by the backend
  n_replicates: number;                // >= 3, enforced by the backend
}

export interface StageProgress {
  stage: "stage1" | "gate1" | "stage2" | "gate2" | "ranking" | "stage3";
  state: "waiting" | "running" | "done" | "blocked";
  done: number; expected: number;
  gpu_hours: number; usd: number;
  started_at: number | null; finished_at: number | null;
}

export interface Funnel {              // per gate, a count of each verdict
  [gate: string]: Partial<Record<Verdict, number>>;
}

export interface TargetReview {        // what studio/review.py produces
  chains: { id: string; n_residues: number }[];
  selected_chain: string | null;
  sequence: string;
  first: number; last: number; n_residues: number;
  real_breaks: { a: number; b: number; distance: number }[];
  numbering_gaps: { a: number; b: number; distance: number }[];
  contig: string;
  atom_coverage: {
    n_residues: number;
    dropped: string[]; missing: string[]; no_ca: string[]; nonstandard: string[];
  };
  chirality: Judgement;
  candidate_partners: { id: string; n_residues: number }[];
  hotspot_readback: { residue: string; label: string | null }[];
  questions: string[];                 // render every one; block submission while non-empty
}
```

Add `src/lib/mock/` holding a canned run that `api.ts` serves when
`process.env.NEXT_PUBLIC_MOCK === "1"`, which is the default for now. Follow the sibling app's
`data:mock` convention: the mock is shaped by real numbers from our own runs, so the UI is
exercised against realistic values.

**The mock must play out over time**, not return a finished run. A run advances stage by stage on
a timer so progress, stalls and failures are all visible in the browser. Include these canned
runs, selectable from a dev-only menu:

| Mock run | What it exercises |
|---|---|
| `happy` | all stages complete; a funnel with survivors at each gate |
| `borderline` | gate-3 designs inside the noise band; must not look like passes |
| `refused` | a backbone sampled 1×, and a gate-3 design with 1 replicate |
| `chirality` | an upload that fails the chirality check before any GPU time |
| `silent` | stage 1 reports success with 12 of 100 structures — run must stop |
| `stall` | no new artifact for 20 minutes, spend still climbing |
| `ceiling` | cost ceiling reached mid-stage-3 |
| `underpowered` | a pass-rate panel whose cells must refuse to be ranked |

## Step 4. Routes

```
src/app/page.tsx                     landing: what it does, what it cannot promise
src/app/runs/page.tsx                the user's runs, newest first, with state and spend
src/app/runs/new/page.tsx            step 1: choose a target (RCSB id or upload)
src/app/runs/new/review/page.tsx     step 2: THE REVIEW SCREEN (see Step 5)
src/app/runs/[runId]/page.tsx        live run: stage progress, funnel, spend, events
src/app/runs/[runId]/designs/page.tsx    per-design table and detail
src/app/runs/[runId]/designs/[design]/page.tsx   one design: structure + every judgement
src/app/runs/[runId]/provenance/page.tsx  the one-page record, with a download
```

Use the sibling app's `AppShell` pattern: a fixed 264px sidebar, a top bar, and a
`max-w-[1320px]` main column. Copy its layout and restyle; do not redesign it.

Carry over its `RunProvider` / `useRun()` / `RunChip` convention and its rule, which matters even
more here: **every page showing a quantitative result carries a chip naming the run it belongs
to, and no component ever hardcodes a data path.** In the pipeline, two runs on one target
silently mixing their pools was a real bug; the UI must make the run visible at all times.

## Step 5. The review screen is the most important page in the app

`studio/review.py` already implements this. Port it faithfully. Four things on it cannot be
defaulted, because a wrong default lets the run finish normally and produce numbers that look
fine:

1. **Which chain is the target.** Show every chain with its residue count, and the selected
   chain's **full sequence** in 60-residue rows. Require an explicit confirmation. Never
   preselect the longest chain: in one of our own targets the longest chain is a different
   protein from the one the file is named after, and every downstream file agreed with the
   mistake.
2. **Which chain defines the site.** List every candidate partner with its contact count.
   Require a choice. Do not rank them by contact count as if that were the answer — on one
   target the highest-contact chain is the wrong domain.
3. **The hotspot set.** For each hotspot the user enters, read the residue label back out of the
   file and show it beside the number. Show the count, with a note that the count changes the
   gate-1 pass rate and so changes the statistics they will compare runs on.

   **Do not put a percentage on this screen.** An earlier version of this brief quoted
   "18.3% at one hotspot, 31.3% at five". Those figures are real and well-powered — equal n=50 in
   all twelve cells, z=+3.69, p≈2e-4 — but they are **pooled over six different targets**, and
   within that batch one target gained 28 points while another gained nothing at all. So as
   guidance about the single target in front of the user they are misleading even though they are
   correct. State it qualitatively: more hotspots aim a run that was aimed badly, and do not lift
   a site that already works.
4. **Anything rebuilt.** Run the chirality check and **refuse**. Never offer a "fix it" button:
   rebuilding is what causes the defect, and no later gate in the pipeline can see it.

Also show, from `atom_coverage`: atoms the pipeline's loader would silently drop, residues with
no CA that would vanish entirely, and non-standard labels that would raise. And show the contig
with each real chain break marked together with the CA–CA distance that decided it, next to the
numbering gaps that are *not* breaks.

**The submit button is disabled while `questions` is non-empty**, and each question is rendered
as a sentence the user can act on. Mirror the CLI's closing line: none of these can be defaulted
safely.

## Step 6. Rendering verdicts honestly

This is where the product earns trust, and where a generic dashboard would quietly lie.

- **Never render a verdict as a boolean or a tick.** Five states, five distinct chips. `refused`
  and `saturated` must not read as either a pass or a fail.
- **Always show the spread.** A gate-3 value is a mean over replicates. Render
  `9.4 ± 2.8 (n=3)` and `4.1 ± 0.6 (n=3)` differently, and draw the threshold with its noise
  band, not as a line. The underlying score is not deterministic: about 1 in 10 repeats flips the
  verdict, so a single number presented as a verdict is a number the measurement cannot support.
- **`detail.tight_floor === false` must be visible** — it means fewer replicates than the number
  that gives a tight noise floor. Say so in words.
- **A pass rate renders as its Wilson interval**, and `underpowered: true` prints the word
  *underpowered* — never as a null or a zero.
- **Two cells whose intervals overlap are not ordered.** Render them as one tier. The backend
  returns tiers; show tiers.
- **There is no sort-by-confidence control anywhere in this app.** Not hidden, not behind a
  toggle, not in a column menu. Pool-wide, the measured correlation between a design's confidence
  and whether it folds the way it was designed to is **r = −0.162** — about 2.6% of the variance,
  which is r² — and the relationship **inverts on 17 of 47 backbones**. Sorting by it is near a
  coin toss. Sortable columns are fold-back RMSD, grip and epitope coverage. If the API ever
  returns a confidence field, display it as context with that caveat beside it, and do not make
  the column sortable.

  (The 2.6% is derived from r, so it appears as a literal in no file; the script that measured r
  is `53_confidence_vs_foldback.py` in the private repository, not in the pipeline folder.)
- Copy the sibling app's honesty mechanism: it shows an amber badge driven by a
  `data_status.illustrative` flag, and a disabled entry for analyses that have not been run so
  the gap is visible in the UI rather than absent from it. Use the same approach for stages not
  reached and for mock data.

## Step 7. Progress, failure and cost

The run page is a monitoring tool, not a spinner. Three failure classes, rendered differently,
because they need different responses from the user:

| Class | What the user sees | Why it differs |
|---|---|---|
| `hard_failure` | "Stage 2 failed to start, retrying (2 of 3)" | transient; the system is handling it |
| `silent_failure` | "Stage 1 reported success but produced 12 of 100 structures — stopped" | **the run stops.** This is never shown as a result |
| `stall` | "No progress for 20 minutes — paused. $18 of your $50 ceiling spent" | the meter is running |

A stage is judged complete by artifact count and elapsed time, never by an exit code — the
backend decides this; the UI must show `done / expected` so the user can see it too.

Show **spend so far against the ceiling, continuously**, on both the run page and the run list.
Rented GPU has no walltime: a stalled job keeps billing. Treat the cost readout as a primary
element, not a footnote.

Render the event log as a plain reverse-chronological list of sentences. Every refusal appears in
it. A quiet run and a run with nine suppressed warnings must not look the same.

## Step 8. Structure viewing

Reuse `molstar-viewer.tsx` as copied. On a design detail page, load the complex and colour the
binder chain and the target chain distinctly, using the tokens in `palette.ts`. Offer a download
of the complex PDB and the binder FASTA.

Our pipeline is known to lose chain labels in one of its outputs, renumbering both chains from 1.
The viewer must therefore take the chain split as data from the API (`binder_length`), not infer
it from chain labels in the file.

## Step 9. Auth, persistence and the API — stub deliberately

Week 1 has no real backend. Define the API surface, implement it against the mock, and leave one
seam:

```
GET    /api/runs
POST   /api/runs                      body: RunSpec        -> Run
GET    /api/runs/:id                  -> Run
GET    /api/runs/:id/progress         -> StageProgress[]
GET    /api/runs/:id/funnel           -> Funnel
GET    /api/runs/:id/judgements       -> Judgement[]
GET    /api/runs/:id/events?since=    -> RunEvent[]
GET    /api/runs/:id/provenance       -> provenance record
POST   /api/targets/review            body: { pdb | rcsb_id, chain?, partner?, hotspots? }
                                      -> TargetReview
```

Implement these as Next route handlers that call into `src/lib/mock/`. Poll `events` and
`progress` on an interval; do not reach for websockets yet. **Do not add authentication,
payments, or a database this week** — put a single `getCurrentUser()` stub in
`src/lib/auth.ts` returning a fixed user, and make every API handler call it, so adding real auth
later is one file.

## Step 10. Verification

The point of this week is a UI that can be judged without a GPU. Finish with:

```bash
npm run lint
npx tsc --noEmit
npm run build
npm run dev
```

Then walk all eight mock runs from Step 3 and confirm each one:

- [ ] `happy` — the funnel adds up; survivors at each gate match the judgement counts
- [ ] `borderline` — a borderline gate-3 design is visually distinct from a pass, and its spread is on screen
- [ ] `refused` — refused designs do not carry forward and are not shown as failures
- [ ] `chirality` — the review screen blocks submission; no run is created
- [ ] `silent` — the run stops, and nothing on the page reads as a result
- [ ] `stall` — the pause is explained, and spend is shown against the ceiling
- [ ] `ceiling` — the run pauses and says what was spent
- [ ] `underpowered` — overlapping cells appear as one tier, and the word *underpowered* is shown
- [ ] no sort-by-confidence control exists anywhere, including column menus
- [ ] every quantitative page carries its run chip
- [ ] every exported file contains the run id, thresholds version and config hash
- [ ] the review screen's submit button is disabled while any question is unanswered

Write `web/README.md` covering: the stack and versions, that the data seam is `src/lib/api.ts`
and nothing fetches around it, how to run in mock mode, the eight mock runs and what each one is
for, and the verdict vocabulary with the rule that the UI never computes one.

---

## What this app must not promise

Put this on the landing page in plain language, and keep it true in every screen. It is not
modesty; it is the product's main credibility asset, and overstating it is the fastest way to
lose a scientific user.

The app can honestly say it runs the stages in the right order with the gates before the
rankings, checks every structure for defects that no other metric in the pipeline can see, shows
the uncertainty on every number, never ranks by confidence, and records exactly what produced
each result.

It must not say, or imply through its layout, that a design will bind. No stage in this pipeline
computes an interaction energy between the binder and the target: confidence that a complex is
*predictable* is standing in for evidence that it *binds*, and those are different claims. So
there is no "best candidates" panel, no score presented as a probability of success, no hit rate,
and no leaderboard. "These passed every gate, here is each gate's reading and its uncertainty" is
the strongest claim the app is entitled to make.
