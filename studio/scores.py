#!/usr/bin/env python3
"""Read an AF2 `.sc` score table by column NAME, never by position.

`41_gate2_report.sh` string-matches the header and then uses `$2` as the RMSD
for every row. A reordered `predict.py` score file would hand it a real number
from a different column with nothing on screen to show it. This module is the
replacement: the header decides the index, and asking for a column that is not
there is an error rather than a silently wrong number.

The files are also append-only, so a re-run leaves several headers and several
runs' rows in one file. `read_sc` handles that and `latest_per` picks the last
row per key, which is what "fresh score filename every run" was working around.

    from studio import scores
    rows = scores.read_sc("sc_pd1_b1_1730000000.sc")
    rmsd = scores.column(rows, "binder_aligned_rmsd", float)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterable

# Columns the pipeline actually depends on. Named here so a missing one is
# reported as a missing column and not as a zero.
GATE2_RMSD = "binder_aligned_rmsd"
GATE3_PAE = "pae_interaction"
DESCRIPTION = "description"


class ScoreFileError(Exception):
    pass


def read_sc(path: Path | str) -> list[dict[str, str]]:
    """Every data row of a (possibly append-only) Rosetta-style score file.

    Each row is a dict keyed by the header that was in force above it. Rows are
    returned in file order; later headers apply to later rows.
    """
    path = Path(path)
    header: list[str] | None = None
    rows: list[dict[str, str]] = []

    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line.startswith("SCORE:"):
            continue
        fields = line.split()[1:]
        if not fields:
            continue
        # A header row is the one whose first field is not a number.
        if _is_header(fields):
            header = fields
            continue
        if header is None:
            raise ScoreFileError(f"{path}:{lineno}: data row before any header")
        if len(fields) != len(header):
            raise ScoreFileError(
                f"{path}:{lineno}: {len(fields)} fields against a "
                f"{len(header)}-column header"
            )
        rows.append(dict(zip(header, fields)))

    if header is None:
        raise ScoreFileError(f"{path}: no SCORE: header found")
    return rows


def _is_header(fields: list[str]) -> bool:
    try:
        float(fields[0])
    except ValueError:
        return True
    return False


def column(
    rows: Iterable[dict[str, str]],
    name: str,
    cast: Callable[[str], Any] = str,
) -> list[Any]:
    """One column, by name. Raises if any row lacks it."""
    out = []
    for i, row in enumerate(rows):
        if name not in row:
            raise ScoreFileError(
                f"row {i}: no column {name!r}; have {', '.join(sorted(row))}"
            )
        out.append(cast(row[name]))
    return out


def latest_per(
    rows: Iterable[dict[str, str]], key: str = DESCRIPTION
) -> dict[str, dict[str, str]]:
    """Last row per key, which in an append-only file is the newest run's."""
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        if key not in row:
            raise ScoreFileError(f"no column {key!r} to key rows on")
        out[row[key]] = row
    return out


def backbone_of(description: str) -> str:
    """`pd1_b1_3_s7_af2pred` -> `pd1_b1_3`.

    Gate 2 decides per backbone using the best of its sequences, so the
    sequence suffix and the predictor suffix both have to come off. Done here
    once instead of in an awk `sub()` in each report script.
    """
    name = description
    for suffix in ("_af2pred", "_af2"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    parts = name.rsplit("_s", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return name
