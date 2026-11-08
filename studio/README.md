# `studio/` — the hosted layer over the numbered pipeline

This is the software architecture for the de novo studio: the layer that lets
someone run the pipeline without a terminal. The product plan it implements is
`/mnt/project-files/studio/DENOVO-STUDIO-PLAN-2026-11-08.md`.

**The numbered scripts in the repository root stay the source of truth.** This
package does not reimplement a measurement that already exists — `checks.py`
shells out to `44_chirality_check.py` rather than writing a second signed-volume
test. A reimplementation that drifts from the scripts would give us two
pipelines and no way to tell which is right.

What this package adds is only what a hosted service needs and a hand-run
pipeline does without.

## Modules

| Module | What it is for |
|---|---|
| `thresholds.toml` | The **only** place a threshold lives. Hashed; every verdict carries the hash. |
| `config.py` | Loads it, validates the sections, computes the hash. |
| `gatekeeper.py` | The **only** thing allowed to write a verdict. |
| `scores.py` | Reads `.sc` tables by column **name**, handles append-only files. |
| `stats.py` | Wilson intervals, replicate spreads, and the unequal-arms refusal. |
| `store.py` | Run state in SQLite, keyed by a run UUID. Survivors, funnel, events, spend. |
| `checks.py` | Chirality, atom-name coverage, break-vs-gap, and stage-output verification. |
| `residues.py` | Canonical heavy-atom names, and what a bare relabel would drop. |
| `provenance.py` | The record a result is auditable from. Refuses to be written incomplete. |
| `review.py` | The stage-0 review screen as a CLI — run this on the login node. |

## Run it

```
python3 tests/test_studio.py                 # 31 tests, stdlib only, no GPU
python3 -m studio.review inputs/il6_clean.pdb --chain A --hotspots 47,155,157,158
```

`review.py` exits 1 when something needs a human answer, so it can gate a
submission:

```
python3 -m studio.review inputs/il6_clean.pdb --chain A --partner B \
  --hotspots 47,155,157,158 && sbatch scripts/31_rfd_multi.sbatch
```

**Not yet validated on a real structure.** The container this was written in
cannot reach RCSB, so every check above is tested against synthetic fixtures
only. The first thing to do with it is run `review.py` over
`inputs/il6_clean.pdb` on the cluster and confirm it reports the known truth:
157 CAs numbered 19–184, one real break 51→61, contig `A19-51/A61-184`, and
`A157` reading back as TRP. If it does not, the tool is wrong, not the file.

## The four rules it is built on

1. **The run UUID is the only directory key.** This structurally removes three
   findings from the 2026-11-07 code review instead of patching six scripts:
   `TAG` never reaching the derived directories, `sc_in/` never being cleared,
   and `gate2/<name>.pass` never being truncated — which meant a run where
   nothing passed printed the *previous* run's survivors as its own.
2. **One gatekeeper writes every verdict**, from one versioned config, and
   stamps its hash on every row. Thresholds stop being awk literals.
3. **Success is output count plus elapsed time**, never an exit code or a
   provider's status.
4. **A verdict is never a bare boolean.** `BORDERLINE` exists because the gate
   is not deterministic where it decides; `REFUSED` exists because "not enough
   evidence" is not a result.

## What it refuses to do

Each of these raises or returns `REFUSED` rather than returning a number, and
none of them has an override flag:

- rank designs by any confidence metric (`plddt`, `iptm`, `mpnn_score`, …)
- write a gate-3 verdict from fewer than 3 replicates
- judge a backbone sampled with fewer than 8 sequences
- compare arms with unequal n
- order two pass-rate cells whose Wilson intervals overlap
- report a score above the saturation point as a number
- fall back to filename order when a ranking key is missing
- write a provenance record with no image digest, weights or seeds
- report a chirality result it could not parse

## Not built yet

The provider driver (RunPod API), the stage executors, the web app, and the
cost-ceiling enforcement loop. Those come after the two questions in §8 of the
plan are answered — above all the PyRosetta one, since `predict.py` imports it
at line 31 and that script runs both gate 2 and stage 3.
