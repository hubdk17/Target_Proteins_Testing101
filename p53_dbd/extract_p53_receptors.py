#!/usr/bin/env python3
"""
Extract clean, single-chain protein receptors for the p53 DBD BoltzGen runs,
from the raw RCSB structures in this directory. Stdlib only; run on the login node.

  cd ~/rfd/scripts && python3 p53_dbd/extract_p53_receptors.py

Writes (same dir):
  p53dbd_wt_receptor.pdb     2OCJ chain A, protein only (wild-type core, 96-289) -> whole-surface run
  p53dbd_y220c_receptor.pdb  2J1X chain A, protein only (Y220C core, 96-290)     -> Y220C pocket run

Keeps chain A only, protein residues only, one conformer (altloc blank/A), heavy atoms.
Drops the structural Zn and the extra crystal copies: BoltzGen holds the target rigid and the
fold is fixed by the coordinates, so we feed it the protein chain exactly like the MDM2 receptor.
Both chains are contiguous (no gaps), so 1-based index = crystal - 95 for 2J1X.
"""
import os

AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
       "MSE", "SEC"}
HERE = os.path.dirname(os.path.abspath(__file__))
JOBS = [("2OCJ.pdb", "p53dbd_wt_receptor.pdb"),
        ("2J1X.pdb", "p53dbd_y220c_receptor.pdb")]

for src, dst in JOBS:
    sp, dp = os.path.join(HERE, src), os.path.join(HERE, dst)
    keep, res = [], set()
    for l in open(sp):
        if l[:6] != "ATOM  ":
            continue
        if l[21] != "A":                       # chain A only
            continue
        if l[16] not in (" ", "A"):            # one conformer
            continue
        if l[12:16].strip().startswith("H") or l[76:78].strip() == "H":
            continue                           # heavy atoms only
        if l[17:20].strip() not in AA3:        # protein residues only
            continue
        keep.append(l.rstrip("\n"))
        res.add(int(l[22:26]))
    with open(dp, "w") as fh:
        fh.write("\n".join(keep) + "\nTER\nEND\n")
    print("%s -> %s : %d atoms, %d residues, span %d-%d"
          % (src, dst, len(keep), len(res), min(res), max(res)))
