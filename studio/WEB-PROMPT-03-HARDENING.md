# Harden the studio app, and add an optional ChimeraX bridge

Third brief, after `WEB-PROMPT.md` and `WEB-PROMPT-02-INTERACTIONS.md`.

Two independent pieces of work. Part A is hardening — what stands between a working demo and
software a stranger can use. Part B is ChimeraX, which has a licence constraint that decides its
architecture, so read B2 before writing any of it.

---

# Part A — hardening

The app passes lint, typecheck, build and a 113-check walkthrough. That makes it *correct*,
against inputs you chose. Hardening is about inputs you did not choose. Work in this order; the
first three are the ones that can hurt a user or the company.

## A1. The upload path is the attack surface

A bench user uploads a structure file. Right now that path presumably parses whatever arrives.
Every item below is a real failure mode for PDB handling, not a hypothetical.

- **Size cap, enforced before parsing**, in the route handler, from a constant. A 300 MB file
  must be rejected on the first chunk, not after it is read into memory.
- **Line-length and column guards.** Fixed-column PDB parsing is the single most repeated source
  of silent corruption in this project. Never slice a line without first asserting it is long
  enough. A short `ATOM` record must raise a named error, never yield a coordinate of `0.0`.
- **Reject, do not coerce.** A non-numeric coordinate, a missing CA, a non-standard residue
  label, an insertion code: each is a refusal with a sentence explaining it. A parser that
  repairs its input is a parser that invents data.
- **Atom-count and residue-count sanity bounds.** A file claiming 2 million atoms, or zero, is
  rejected before anything allocates.
- **Never execute or shell out with any part of an uploaded filename.** Store under the run UUID
  with a generated name; keep the original only as a display string, escaped.
- **Path traversal.** `../` in any user-supplied name, and any absolute path, is rejected. Every
  file write resolves inside the run directory and the code asserts that after resolution, not
  before.
- **Decompression.** If you accept `.gz` or `.cif.gz`, cap the *decompressed* size too, or a
  small upload becomes a disk-filling one.

## A2. Every API boundary validates its own input

The handlers currently serve a mock, so they trust their callers. They will not be able to.

- **Parse and validate every request body against a schema at the boundary** — Zod or equivalent
  — and return a 400 with the field name. Do not hand an unvalidated object to anything.
- **`RunSpec` has hard floors in the backend**: `n_sequences_per_backbone >= 8`,
  `n_replicates >= 3`. These are not UI hints; a crafted POST must be refused with the reason.
  The same for a cost ceiling of zero or a negative design count.
- **Every `[runId]` and `[design]` route checks that the run belongs to the current user**
  before returning anything, including in the mock. Add it now, while `getCurrentUser()` is a
  stub, because a missing ownership check is invisible until the day it matters.
- **Treat a run id as opaque.** Validate it is a UUID and look it up; never interpolate it into a
  path or a query.
- **Rate-limit the upload and review endpoints**, even crudely, in memory. Structure parsing is
  the most expensive thing the server does without a GPU.

## A3. Money and state

- **The cost ceiling must be enforced server-side, in the orchestrator**, and the UI must only
  display it. A ceiling that lives in the client is not a ceiling.
- **Make the ceiling check happen before each unit of spend, not after**, and store the decision.
  "Paused at $49.93 of $50" was the right catch; the fix is that the check runs before dispatch
  and the banner states whichever actually happened — *paused below your ceiling* or *ceiling
  reached* — never the wrong one of the two.
- **Every state transition is a write, and illegal transitions are refused.** A run cannot go
  from `finished` to `running`. Encode the allowed transitions as a table and reject the rest;
  this is where a retry loop will otherwise corrupt a run.
- **Idempotent run creation.** A double-clicked submit must not create two runs and two bills.
  Accept a client-generated idempotency key.

## A4. Failure, without losing the user

- **An error boundary per route**, not one for the whole app, each saying what failed and what is
  unaffected. A blank page is the one outcome a monitoring tool may never produce.
- **Polling must survive the server going away.** On a failed poll, back off, show *last updated
  N minutes ago* rather than stale data presented as live, and recover without a reload.
- **Never render stale data as current.** If `events` has not advanced, say so. This is the UI
  version of the pipeline bug where a run that passed nothing printed the previous run's
  survivors.
- **A timeout on every outbound call**, including the eventual provider API. An un-timed-out
  fetch to a rented-GPU provider is how a page hangs forever and a meter keeps running.

## A5. The mock-to-real swap is the highest-risk moment in the project

Plan for it now, because it is the step where a demo quietly becomes wrong.

- **One flag, read in one place.** `NEXT_PUBLIC_MOCK` must be read only inside `src/lib/api.ts`.
  Grep for it as a test: more than one hit is a bug.
- **Mock data must be unmistakable on screen.** Follow the sibling app's amber badge convention.
  A screenshot of mock output must never be presentable as a result.
- **Make the real adapter fail loudly when the flag is off and no backend is configured** —
  a clear error, never a silent fall back to the mock. A hosted app serving mock numbers while
  looking live is the worst outcome available here.
- **Keep the mock in the test suite after the swap.** It is the only way to exercise the stall,
  silent-failure and ceiling paths on demand.

## A6. Tests that can fail

The project's own rule, and the walkthrough already proved its worth twice.

- **Test mechanisms, not words.** Already learned the hard way four times: grepping rendered
  output for a forbidden term matches the sentence forbidding it, and server-rendered HTML is one
  long line with `<!-- -->` between text nodes, so an unanchored pattern matches across the whole
  document. Assert on a control's markup, or on the absence of a handler, never on prose.
- **A fixture that cannot fail is worse than no fixture.** The first chirality fixture written for
  `studio/` had its sign convention backwards, so the mirror-image case tested nothing. For every
  guard, write the input that *should* trip it and assert that it does.
- **Add a property test over the PDB parser**: random well-formed records round-trip; truncated,
  reordered and garbage records raise a named error and never return a coordinate.
- **Pin the chain-convention trap with an assertion.** `42_grip_rank.py` takes chain **A** as the
  target (`GRIP_CUT = 8.0`, `FOOT_CUT = 10.0`, CA atoms only, lines 27–28), while an
  `*_af2pred.pdb` has the binder as chain A because PyRosetta renumbers the pose from 1. Grip and
  footprint counts survive the swap, being symmetric pair counts; **hotspot matching does not**,
  and would silently match almost nothing. Assert the convention at the parse boundary and test
  that relabelling the chains in a fixture changes nothing downstream.
- **Report nearest heavy-atom distance beside any CA–CA figure.** CA–CA understates a contact at
  the edge of a shell: one design reads 7.98 Å CA–CA to a hotspot while carrying 13 binder heavy
  atoms within 5 Å of it. `interface_report.py` in the private repo sets the convention — heavy
  atoms within 5 Å, closest approach, contacts within 4.5 Å.

## A7. Honesty regressions are bugs

Add a test for each, because these are the claims the product is sold on:

- no sort-by-confidence control, asserted on markup
- no gate verdict rendered as a bare boolean anywhere
- no gate-3 value displayed without its replicate count and spread
- no pass rate displayed without its interval
- no interface energy displayed without its settings object
- `refused` and `saturated` never styled as either a pass or a fail
- a single-design delta inside the noise floor is never captioned as *no dependence*; it reads as
  *no measurable dependence*, which is a different claim

## A8. Operational basics

Security headers (CSP, `X-Content-Type-Options`, `Referrer-Policy`), no secrets in
`NEXT_PUBLIC_*`, structured server logs with the run id on every line and never a token, a
`/api/health` that checks its dependencies rather than returning 200 unconditionally, and
dependency versions pinned with a lockfile committed.

---

# Part B — the ChimeraX bridge

## B1. It is technically straightforward, and the catch is CORS

ChimeraX ships a REST server, so a web page really can drive a copy running on the user's own
machine:

```
remotecontrol rest start port 49200 json true cors true
```

Then `GET http://127.0.0.1:49200/run?command=open+design_36_s3.pdb`, and with `json true` the
response carries `python values`, `json values`, `log messages` and an `error` object with `type`
and `message`. Large payloads go by POST as `multipart/form-data`. The port is random unless you
give one; 49152–65535 is preferred.

**The constraint that decides the architecture: CORS is off by default, and when enabled it is
restricted to localhost origins** — UCSF describes that explicitly as a security measure.

So:

| Where the page is served from | Works? |
|---|---|
| `http://localhost:3000` (dev, or a locally-run studio) | **yes**, with `cors true` |
| `https://studio.example.com` (hosted) | **no** — not a localhost origin |

That is the whole story, and it means **the live bridge cannot be a feature of the hosted
product as-is.** Do not build it that way and discover this later.

Three options, in the order I would do them:

1. **Ship a `.cxc` download button. Do this first.** A generated ChimeraX command script: open
   the complex, colour binder and target, select the hotspots, set the view. Zero install, zero
   bridge, zero licence exposure, works from a hosted page today, and the user double-clicks it.
   There is precedent on disk already — `chimerax_session.cxc` in the PPI folder and
   `rfdiffusion-il6/chimerax/view_short9.cxc` in the project folder. Generate it server-side from
   the same data the Interface page renders.
   **Two things the generator must get right:** ChimeraX on Windows accepts forward slashes, and
   a script that opens a structure by relative name needs ChimeraX's working directory set first
   — so emit a `cd` line and keep script and structure together. And respect the chain
   convention: in an `*_af2pred.pdb` the binder is chain A and the target chain B, the opposite
   of the target-only file.
2. **Live control for a locally-run studio only**, behind a capability check: probe
   `remotecontrol rest port`, and if it does not answer, show the download button instead of an
   error. Never block a page on it.
3. **A small local helper**, only if users actually ask for live control in the hosted product.
   Our own process on 127.0.0.1 that serves our origin with its own CORS headers and forwards to
   ChimeraX. It works because the headers are ours to set, but it is an installer, an updater and
   a support burden. Do not start here.

Whatever is built: **never send a ChimeraX command string assembled from unvalidated input.**
`remotecontrol` executes ChimeraX commands, which include file and shell-adjacent operations.
Build commands from a fixed template with validated parameters only.

## B2. The licence position — and it is the same shape as PyRosetta

Read from the UCSF licensing page on 2026-11-10:

- ChimeraX is **"licensed for non-commercial use at no cost by The Regents"**, and you accept the
  UCSF ChimeraX Non-Commercial License Agreement before downloading.
- **"COMMERCIAL USE OF UCSF CHIMERAX IS COVERED BY A SEPARATE WRITTEN LICENSE AGREEMENT."**
  A commercial licence is available; pricing depends on the number of users.
- The page **does not define commercial use**, and **says nothing about redistribution, bundling,
  or offering ChimeraX as part of a service.**

What follows, and what does not:

- **Right now, in academic work, you are fine.** Non-commercial use at no cost covers everything
  this project does today.
- **Three things are clearly safe even commercially**, because they do not use ChimeraX at all:
  generating a `.cxc` file (a text file of commands — we author it, UCSF does not own command
  names), keeping **Mol\*** as the in-browser viewer (MIT, already in the app, no restriction),
  and writing documentation that says "open this in ChimeraX".
- **Three things are not safe without a written answer**: bundling or shipping ChimeraX with the
  product, running it server-side to render images for customers, and making a paying customer's
  use of ChimeraX a required part of a commercial workflow — because their free licence is
  non-commercial, and a commercial customer using it for commercial work needs their own
  commercial licence whether or not we are involved.

**The recommendation.** Keep ChimeraX strictly as something the *user* has installed and licensed
themselves, make the integration optional and never required, and never ship, proxy or
server-side-run it. Mol\* stays the primary viewer so no feature depends on ChimeraX. On that
footing the `.cxc` download is the right first step and carries no exposure worth worrying about.

Then ask UCSF two specific questions in writing, the same way we need to ask UCSF about PyRosetta
and RunPod about hosting — this is now the third item on that list, so send them together:

1. Does a commercial ChimeraX licence, or any licence, permit a hosted commercial service to send
   commands to a copy of ChimeraX that the end user has installed and licensed themselves?
2. Does generating and distributing `.cxc` command scripts, without distributing ChimeraX,
   require any licence?

Ask before it becomes a paid feature, not after. And note the asymmetry that makes this easy:
**Mol\* is MIT and already works, so ChimeraX is a convenience for users who prefer it, never a
dependency.** Keeping it that way is what keeps the licence question cheap.

## B3. Verification

- [ ] a `.cxc` download works from a hosted build, with no ChimeraX present
- [ ] the generated script opens, colours and selects correctly on Windows, forward slashes and
      `cd` line included
- [ ] hotspot selections use the file's own chain convention, verified by residue *name* and not
      by number
- [ ] with ChimeraX absent, the live-control probe degrades to the download button and logs
      nothing alarming
- [ ] no ChimeraX command is ever built by string concatenation from user input
- [ ] Mol\* still renders every structure with ChimeraX uninstalled
