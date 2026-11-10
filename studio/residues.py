#!/usr/bin/env python3
"""Canonical heavy-atom names per residue label.

Why this table is here. The stage-1 loader
(`rfdiffusion/inference/utils.py: parse_pdb_lines`) matches each ATOM record
**by name against the canonical tuple for the label on its own record**, and
only the first 14 heavy-atom slots are consulted. An atom whose name is absent
from that tuple falls out of the loop: its slot stays NaN, the mask stays False
and the coordinate is zeroed. There is no positional fallback, no warning and no
counter. So an uploaded structure can lose atoms on the way in and nothing says
so.

Three consequences this table lets the studio check for before spending GPU time:

1. An atom name not in the label's tuple is **dropped silently**.
2. Renaming a label without editing its atoms gives a **chimera**: only the
   names shared between the old and new tuple survive. TRP -> TYR keeps
   CG, CD1, CD2, CE2; drops NE1, CE3, CZ2, CZ3, CH2; and never supplies CE1,
   CZ, OH.
3. A residue with **no CA record** never enters the residue list at all and is
   absent from the loaded structure, side chain included.

Validated against the worked TRP -> TYR example in
`multi-target/PARSER-ATOM-NAME-LIMITATION.md`. Re-check against the cluster copy
before trusting it for a new residue type:

    grep -n "aa2long" ~/rfd/RFdiffusion/rfdiffusion/chemical.py
"""
from __future__ import annotations

BACKBONE = ("N", "CA", "C", "O")

CANONICAL: dict[str, tuple[str, ...]] = {
    "ALA": BACKBONE + ("CB",),
    "ARG": BACKBONE + ("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"),
    "ASN": BACKBONE + ("CB", "CG", "OD1", "ND2"),
    "ASP": BACKBONE + ("CB", "CG", "OD1", "OD2"),
    "CYS": BACKBONE + ("CB", "SG"),
    "GLN": BACKBONE + ("CB", "CG", "CD", "OE1", "NE2"),
    "GLU": BACKBONE + ("CB", "CG", "CD", "OE1", "OE2"),
    "GLY": BACKBONE,
    "HIS": BACKBONE + ("CB", "CG", "ND1", "CD2", "CE1", "NE2"),
    "ILE": BACKBONE + ("CB", "CG1", "CG2", "CD1"),
    "LEU": BACKBONE + ("CB", "CG", "CD1", "CD2"),
    "LYS": BACKBONE + ("CB", "CG", "CD", "CE", "NZ"),
    "MET": BACKBONE + ("CB", "CG", "SD", "CE"),
    "PHE": BACKBONE + ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "PRO": BACKBONE + ("CB", "CG", "CD"),
    "SER": BACKBONE + ("CB", "OG"),
    "THR": BACKBONE + ("CB", "OG1", "CG2"),
    "TRP": BACKBONE + ("CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "TYR": BACKBONE + ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"),
    "VAL": BACKBONE + ("CB", "CG1", "CG2"),
}

# The loader reads at most 14 heavy-atom slots. TRP uses all 14.
MAX_HEAVY_SLOTS = 14

ONE_LETTER = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}
THREE_LETTER = {v: k for k, v in ONE_LETTER.items()}


def survives_relabel(old: str, new: str) -> tuple[list[str], list[str], list[str]]:
    """What a bare label change would keep, drop and never supply.

    Returns (kept, dropped, never_supplied). A non-empty `dropped` or
    `never_supplied` means the result is a chimeric side chain, not a
    substitution: a label change must always be paired with an explicit atom set.
    """
    for label in (old, new):
        if label not in CANONICAL:
            raise KeyError(f"{label!r} is not a canonical residue label")
        if len(CANONICAL[label]) > MAX_HEAVY_SLOTS:
            raise AssertionError(f"{label} exceeds {MAX_HEAVY_SLOTS} heavy slots")
    a, b = set(CANONICAL[old]), set(CANONICAL[new])
    return sorted(a & b), sorted(a - b), sorted(b - a)
