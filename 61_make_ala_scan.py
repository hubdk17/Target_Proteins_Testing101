#!/usr/bin/env python3
"""
61_make_ala_scan.py -- build an exact alanine-scan series from a complex.

Why this is exact, and why that matters
---------------------------------------
Mutating any residue to alanine means deleting every side-chain heavy atom
beyond CB and renaming the residue. No atom is ever *added*. So unlike every
other mutation we have built in this project, an alanine scan needs no
PDBFixer, no tleap, no rotamer packing and no loop rebuilding -- which means
it cannot introduce a D-amino acid, and there is no rebuilt geometry to
distrust. (PDBFixer put five D-residues into 4ZQK PD-1 85-92 and every gate
was blind to them.)

It also sidesteps the parser trap from the IL-6 work: that experiment could
not build "label ALA with full atoms" because RFdiffusion's parser matches
atoms by name against the label's canonical tuple. Here we do the opposite
and fully supported thing -- label ALA with alanine's atoms.

The approximation that IS being made, stated plainly: the backbone and all
other side chains stay exactly where the crystal put them. No relaxation.
This is standard computational alanine scanning. It will systematically
over-read the loss from a large buried residue, because the cavity left
behind never closes.

Glycine and proline are skipped by default: GLY has nothing to delete, and
PRO -> ALA removes a ring bonded to the backbone N, which the fixed-backbone
assumption cannot honestly represent.

Usage
-----
  # scan every residue of chain B
  python3 61_make_ala_scan.py 1YCR.pdb --chain B --outdir arms/

  # scan only named positions
  python3 61_make_ala_scan.py 1YCR.pdb --chain B --positions 19,22,23,26 \
      --outdir arms/

  # also emit negative-control arms at positions with zero interface
  # contact, which SHOULD score identically -- a gate that moves on these
  # is reporting noise
  python3 61_make_ala_scan.py 1YCR.pdb --chain B --positions 21,24 \
      --outdir arms/ --tag control

Every output keeps the wild-type file alongside it, so the comparison is
always made within one batch rather than against a number from a previous
run on a different node.
"""

import argparse
import os
import sys

KEEP = {"N", "CA", "C", "O", "CB", "OXT"}
SKIP_RESNAMES = {"GLY", "PRO"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def read_lines(path):
    with open(path) as fh:
        return fh.readlines()


def residue_index(lines):
    """(chain, resnum) -> resname, for ATOM records only."""
    out = {}
    for ln in lines:
        if not ln.startswith("ATOM"):
            continue
        out[(ln[21], int(ln[22:26]))] = ln[17:20].strip()
    return out


def mutate_to_ala(lines, chain, resnum):
    """Return new lines with (chain, resnum) reduced to alanine.

    Deletes side-chain heavy atoms past CB and any hydrogens on them, and
    rewrites the residue name in every record that mentions the residue.
    Atom serial numbers are left alone deliberately: nothing downstream
    reads them, and renumbering is a silent way to break a reader.
    """
    out = []
    n_deleted = 0
    for ln in lines:
        if ln.startswith(("ATOM", "ANISOU", "HETATM")):
            try:
                num = int(ln[22:26])
            except ValueError:
                out.append(ln)
                continue
            if ln[21] == chain and num == resnum:
                name = ln[12:16].strip()
                if name not in KEEP:
                    n_deleted += 1
                    continue
                ln = ln[:17] + "ALA" + ln[20:]
        out.append(ln)
    return out, n_deleted


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build an exact fixed-backbone alanine scan.")
    ap.add_argument("pdb")
    ap.add_argument("--chain", required=True,
                    help="Chain to scan (the binder/peptide).")
    ap.add_argument("--positions", default=None,
                    help="Comma-separated residue numbers. Default: every "
                         "residue of the chain.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", default="ala",
                    help="Filename tag, e.g. 'control' for the "
                         "no-contact negative-control arms.")
    ap.add_argument("--include-pro", action="store_true",
                    help="Also mutate proline. Off by default: the "
                         "fixed-backbone assumption cannot represent it.")
    args = ap.parse_args(argv)

    lines = read_lines(args.pdb)
    index = residue_index(lines)
    chain_res = sorted(n for (c, n) in index if c == args.chain)
    if not chain_res:
        raise SystemExit("chain %r not found; file has %r"
                         % (args.chain,
                            sorted({c for (c, _) in index})))

    if args.positions:
        want = [int(x) for x in args.positions.split(",") if x.strip()]
        bad = [p for p in want if p not in chain_res]
        if bad:
            raise SystemExit("positions not in chain %s: %r"
                             % (args.chain, bad))
    else:
        want = chain_res

    os.makedirs(args.outdir, exist_ok=True)

    # the wild type always travels with its arms, in the same directory
    wt = os.path.join(args.outdir, "wt.pdb")
    with open(wt, "w") as fh:
        fh.writelines(lines)
    made = [("wt", "-", 0)]

    skipped = []
    for pos in want:
        resname = index[(args.chain, pos)]
        if resname == "ALA":
            skipped.append("%s%d (already ALA)" % (resname, pos))
            continue
        if resname in SKIP_RESNAMES and not (
                resname == "PRO" and args.include_pro):
            skipped.append("%s%d (skipped by default)" % (resname, pos))
            continue
        new, ndel = mutate_to_ala(lines, args.chain, pos)
        if ndel == 0:
            skipped.append("%s%d (nothing to delete)" % (resname, pos))
            continue
        one = THREE_TO_ONE.get(resname, "X")
        name = "%s_%s%dA.pdb" % (args.tag, one, pos)
        with open(os.path.join(args.outdir, name), "w") as fh:
            fh.writelines(new)
        made.append((name, "%s%d->ALA" % (resname, pos), ndel))

    print("wrote %d structures to %s" % (len(made), args.outdir))
    print("%-24s %-16s %s" % ("file", "mutation", "atoms deleted"))
    for name, mut, ndel in made:
        print("%-24s %-16s %d" % (name, mut, ndel))
    if skipped:
        print("\nskipped:")
        for s in skipped:
            print("  " + s)
    print("\nNext: score the whole directory in ONE batch, wild type "
          "included:")
    print("  python3 60_gate_battery.py %s/*.pdb --target-chain A "
          "--binder-chain %s --out gates_scan.tsv"
          % (args.outdir.rstrip("/"), args.chain))
    return 0


if __name__ == "__main__":
    sys.exit(main())
