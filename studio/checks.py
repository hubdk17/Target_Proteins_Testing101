#!/usr/bin/env python3
"""Structure validation on the way in, and stage-output verification on the way out.

These are the checks a bench user should never have to ask for. All of them are
cheap, and each one catches something no other gate in the pipeline can see.

- `chirality()` shells out to the repo's own `44_chirality_check.py` rather than
  reimplementing the signed-volume test. One implementation of a measurement,
  per the 2026-11-07 review's de-duplication finding.
- `atom_coverage()` reports atoms the stage-1 loader would silently drop, and
  residues with no CA record, which the loader omits entirely.
- `chain_walk()` measures consecutive CA-CA distances so a numbering gap and a
  real chain break are distinguished by measurement, never by assumption.
- `verify_output()` decides whether a stage really succeeded. A stage's success
  is its output count and elapsed time, never its exit code or a provider's
  status: an empty arm that was a container which never started looked like a
  success once already.
"""
from __future__ import annotations

import math
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from .residues import CANONICAL

REPO_ROOT = Path(__file__).resolve().parent.parent
CHIRALITY_SCRIPT = REPO_ROOT / "44_chirality_check.py"


class CheckError(Exception):
    pass


# --------------------------------------------------------------------- chirality

@dataclass(frozen=True)
class Chirality:
    n_checked: int
    n_d: int
    d_residues: list[str]
    mirrored_warning: bool
    raw: str


_TESTED = re.compile(r"tested\s+(\d+)\s+residues")
_D_COUNT = re.compile(r"\*\*\*\s+(\d+)\s+D-AMINO ACID\(S\)")
_D_LINE = re.compile(r"chain\s+(\S+)\s+(\S+)\s+V\s*=\s*([-+0-9.]+)")
_NEGATIVE_MAJORITY = re.compile(r"\*\*\*\s+(\d+)\s+of\s+(\d+)\s+residues read NEGATIVE")


def chirality(path: Path | str, script: Path = CHIRALITY_SCRIPT) -> Chirality:
    """Run the repo's gate-0 script and parse its report.

    Raises rather than returning zero if the output cannot be parsed: a check
    that silently reports "no problems" because its parse failed is exactly the
    failure mode this project catalogues as "a check that cannot fail".
    """
    if not script.exists():
        raise CheckError(f"{script} not found; cannot run the chirality check")
    proc = subprocess.run(
        [sys.executable, str(script), str(path)],
        capture_output=True, text=True,
    )
    out = proc.stdout

    tested = _TESTED.search(out)
    if not tested:
        raise CheckError(
            f"could not parse a residue count out of {script.name}'s output; "
            f"refusing to report a result.\n--- stdout ---\n{out}\n"
            f"--- stderr ---\n{proc.stderr}"
        )
    n_checked = int(tested.group(1))

    majority = _NEGATIVE_MAJORITY.search(out)
    if majority:
        return Chirality(n_checked, int(majority.group(1)), [], True, out)

    d_match = _D_COUNT.search(out)
    n_d = int(d_match.group(1)) if d_match else 0
    residues = [f"chain {m.group(1)} {m.group(2)} V={m.group(3)}"
                for m in _D_LINE.finditer(out)] if n_d else []

    # The script exits 1 on any finding. If it did and we parsed none, the parse
    # is wrong, not the structure.
    if proc.returncode != 0 and n_d == 0:
        raise CheckError(
            f"{script.name} exited {proc.returncode} but no finding was parsed; "
            f"refusing to report a pass.\n--- stdout ---\n{out}"
        )
    return Chirality(n_checked, n_d, residues, False, out)


# ------------------------------------------------------------------ atom names

@dataclass
class AtomCoverage:
    n_residues: int = 0
    dropped: list[str] = field(default_factory=list)       # name not in the tuple
    missing: list[str] = field(default_factory=list)       # tuple name absent
    no_ca: list[str] = field(default_factory=list)         # vanishes at load
    nonstandard: list[str] = field(default_factory=list)   # KeyError at load

    @property
    def clean(self) -> bool:
        return not (self.dropped or self.no_ca or self.nonstandard)

    def summary(self) -> str:
        if self.clean and not self.missing:
            return f"{self.n_residues} residues, all atom names canonical"
        bits = [f"{self.n_residues} residues"]
        if self.nonstandard:
            bits.append(f"{len(self.nonstandard)} non-standard label(s) - these raise at load")
        if self.no_ca:
            bits.append(f"{len(self.no_ca)} residue(s) with no CA - these vanish at load")
        if self.dropped:
            bits.append(f"{len(self.dropped)} atom(s) the loader would drop silently")
        if self.missing:
            bits.append(f"{len(self.missing)} canonical atom(s) absent (incomplete side chains)")
        return "; ".join(bits)


def atom_coverage(path: Path | str, chain: str | None = None) -> AtomCoverage:
    """What would actually arrive in the stage-1 tensor, versus what is in the file."""
    seen: dict[tuple, dict] = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if chain is not None and line[21] != chain:
            continue
        if len(line) < 54:
            raise CheckError(f"{path}: ATOM record shorter than 54 columns:\n{line}")
        label = line[17:20].strip()
        key = (line[21], line[22:26].strip(), line[26])
        name = line[12:16].strip()
        rec = seen.setdefault(key, {"label": label, "atoms": set()})
        rec["atoms"].add(name)

    cov = AtomCoverage(n_residues=len(seen))
    for (ch, num, icode), rec in seen.items():
        where = f"{rec['label']}{num}{icode.strip()} chain {ch}"
        tuple_ = CANONICAL.get(rec["label"])
        if tuple_ is None:
            cov.nonstandard.append(where)
            continue
        if "CA" not in rec["atoms"]:
            cov.no_ca.append(where)
        for name in sorted(rec["atoms"]):
            if name not in tuple_:
                cov.dropped.append(f"{where}: {name}")
        for name in tuple_:
            if name not in rec["atoms"]:
                cov.missing.append(f"{where}: {name}")
    return cov


# ------------------------------------------------------------------ chain walk

@dataclass
class ChainWalk:
    chain: str
    first: int
    last: int
    n_residues: int
    real_breaks: list[tuple[int, int, float]]      # (a, b, CA-CA distance)
    numbering_gaps: list[tuple[int, int, float]]   # number jumps that are NOT breaks
    contig: str


def chain_walk(path: Path | str, chain: str, break_cutoff: float = 4.5) -> ChainWalk:
    """Separate a real chain break from a numbering convention, by measurement.

    A true peptide bond is ~3.8 A between consecutive CA atoms. A jump in
    residue number with the CAs still ~3.8 A apart is only numbering and the
    contig stays whole; a CA-CA distance above the cutoff is a real break and
    the contig splits there. IL-6 has one real break (51->61) and trypsin had
    six numbering gaps and zero breaks.
    """
    cas: list[tuple[int, tuple[float, float, float]]] = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        if line[12:16].strip() != "CA":
            continue
        cas.append((int(line[22:26]),
                    (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    if not cas:
        raise CheckError(f"{path}: no CA atoms in chain {chain}")

    real, gaps, segments = [], [], [[cas[0][0]]]
    for (na, xa), (nb, xb) in zip(cas, cas[1:]):
        d = math.dist(xa, xb)
        if d > break_cutoff:
            real.append((na, nb, d))
            segments[-1].append(na)
            segments.append([nb])
        elif nb - na > 1:
            gaps.append((na, nb, d))
    segments[-1].append(cas[-1][0])

    contig = "/".join(f"{chain}{a}-{b}" for a, b in segments)
    return ChainWalk(chain, cas[0][0], cas[-1][0], len(cas), real, gaps, contig)


# ------------------------------------------------------- stage output verification

@dataclass(frozen=True)
class OutputCheck:
    ok: bool
    produced: int
    expected: int
    elapsed_s: float
    kind: str      # 'finished' | 'silent_failure' | 'hard_failure'
    message: str


def verify_output(expected: int, produced: int, elapsed_s: float,
                  exit_code: int | None = None,
                  min_seconds_per_item: float = 0.0) -> OutputCheck:
    """Did this stage really succeed?

    Judged on artifact count and elapsed time, never on the exit code. Two real
    cases this covers: a container that never started produced an empty arm that
    read as a success, and a per-design rate once came from a walltime cap
    divided by a killed job's output, which is a ceiling and not a rate.
    """
    if produced >= expected and elapsed_s >= expected * min_seconds_per_item:
        return OutputCheck(True, produced, expected, elapsed_s, "finished",
                           f"{produced} of {expected} artifacts in {elapsed_s:.0f}s")
    if produced == 0:
        return OutputCheck(
            False, produced, expected, elapsed_s, "hard_failure",
            f"no artifacts at all after {elapsed_s:.0f}s"
            + (f" (exit {exit_code})" if exit_code is not None else "")
            + " - treat as a stage that never started",
        )
    if produced < expected:
        return OutputCheck(
            False, produced, expected, elapsed_s, "silent_failure",
            f"only {produced} of {expected} artifacts"
            + (" despite exit 0" if exit_code == 0 else "")
            + f" after {elapsed_s:.0f}s - stopping, this is not a result",
        )
    return OutputCheck(
        False, produced, expected, elapsed_s, "silent_failure",
        f"{produced} artifacts but only {elapsed_s:.0f}s elapsed, under the "
        f"{expected * min_seconds_per_item:.0f}s this stage needs - the outputs "
        "are probably stale or copied",
    )
