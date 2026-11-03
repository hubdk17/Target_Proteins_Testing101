# Runbook: score the finished m1 arms

Takes the arm groups whose stage-1 pools are complete through gate 1, stage 2
(ProteinMPNN), threading, complex assembly and stage 3 (AF2 initial guess), then
reads the result out.

Groups covered: `casp3`, `pcsk9`. **The 1T29 group is left out** — its pools are
short.

Every step is one paste-and-go block. Run them in order and check the printed
counts before moving on; a step's success is its output files, never its exit
code or its Slurm state.

## 0. Before submitting anything

```bash
module load slurm
squeue -u dakshcsed -o "%.10i %.12j %.4t %.10M %R"
```

Two other threads submit to the same account. If a GPU job of theirs is
running, note which node, and pin step 4 somewhere else.

## 1. Gate 1 on each arm (login node, seconds)

```bash
cd ~/rfd && git -C scripts pull
cd ~/rfd && TAG=m1 bash scripts/32_gate1_compactness.sh \
  casp3_wt casp3_mut casp3_ctl pcsk9_wt pcsk9_mut pcsk9_ctl
```

The printed `rate` is the gate-1 pass rate over the **whole** pool. It is read
here and nowhere later, because the next step cuts the pools down.

## 2. Equal n across the arms of each group (login node, instant)

```bash
cd ~/rfd && python3 -m py_compile scripts/60_equal_n_arms.py && echo INTACT
cd ~/rfd && TAG=m1 python3 scripts/60_equal_n_arms.py casp3 pcsk9
```

Writes `gate1/<group>_<arm>.eqn` and `inputs/arms_m1.tsv`. It prints pool, pass
and kept per arm; the record is in `gate1/equal_n_m1.tsv`. If it says the arms
differ by more than half, top the short arm up with a stage-1 resubmit under the
same `TAG` (existing outputs are skipped, so that tops up rather than redoes)
before going on.

## 3. Stage 2, sequences (array, one task per arm)

```bash
cd ~/rfd && TABLE=$HOME/rfd/inputs/arms_m1.tsv PASSEXT=eqn TAG=m1 \
  sbatch --array=1-6%2 --exclude=dgxanode03 scripts/33_mpnn_multi.sbatch
```

`PASSEXT=eqn` is what makes it read the equal-n lists instead of the full gate-1
lists. `%2` throttles to two at a time; an array left to pile onto one node
costs throughput.

Then thread and assemble, both login-node:

```bash
cd ~/rfd && PASSEXT=eqn python3 scripts/thread_multi.py \
  casp3_wt casp3_mut casp3_ctl pcsk9_wt pcsk9_mut pcsk9_ctl
cd ~/rfd && python3 scripts/43_build_complex.py \
  casp3_wt casp3_mut casp3_ctl pcsk9_wt pcsk9_mut pcsk9_ctl
```

## 4. Stage 3, AF2 on the complex (array, one task per arm)

Pick one node for a whole group and pin it:

```bash
cd ~/rfd && TABLE=$HOME/rfd/inputs/arms_m1.tsv TAG=m1 \
  sbatch --array=1-6%2 --nodelist=dgxanode04 --exclude=dgxanode03 \
  scripts/62_af2_arms.sbatch
```

- `--nodelist` is not optional. AF2 here is not deterministic even on one node,
  so letting the arms land on different hardware adds a between-node term to the
  difference being read.
- `--exclude=dgxanode03` always: its container disk fails intermittently
  unpacking the squashfs, FAILED inside two minutes with a zero-byte log.
- Each task writes its own timestamped `.sc`. The score file is append-only, so
  a reused name silently interleaves two runs.

## 5. Read it out (login node, instant)

```bash
cd ~/rfd && TAG=m1 python3 scripts/63_pae_report.py casp3 pcsk9
```

Prints, per arm: n, the count below `pae_interaction` 10, the mean with its
standard error, and whether an arm-vs-arm difference is resolvable.

What the numbers can and cannot carry:

- **The gate is a floor, not a ranking.** Below 10 passes. Among designs that
  fail it the value does not order them.
- **A single design is not a measurement.** About one in ten flips either side
  of the gate on a re-run, and roughly 3 units separate two single designs. Only
  a mean over about 9 has a floor near 0.5, so the script refuses to call a
  difference when either arm is under 9.
- **Columns are read by name.** `grep -c '^SCORE:'` counts the header too, so it
  overcounts designs by one.
- A score file holding a duplicate design name was appended to twice; the script
  skips that arm rather than printing a mean under a warning.
