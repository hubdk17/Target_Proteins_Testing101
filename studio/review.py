#!/usr/bin/env python3
"""The stage-0 review screen, as a command line tool.

This is what the studio will put in front of a user before any GPU time is
spent, printed instead of rendered. Everything here is a measurement or a
question; nothing is decided for you, because each of these four is a choice
where a wrong default runs to completion and produces plausible numbers:

  - which chain is the target        (picking by length called 4ZQK's PD-L1 "pd1")
  - which partner defines the site   (picking by contact count chose PCSK9's prodomain)
  - the hotspot set                  (1 -> 5 hotspots moved gate-1 pass 18.3% -> 31.3%)
  - whether to trust rebuilt atoms   (five D-residues in one rebuilt loop, invisible downstream)

Usage, on the login node, stdlib only, no GPU:

    python3 -m studio.review inputs/il6_clean.pdb --chain A
    python3 -m studio.review inputs/pd1_clean.pdb --chain A --partner B

Exit status is 1 if anything needs a human answer, so it can gate a submission.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import checks
from .config import load
from .gatekeeper import Gatekeeper, Verdict
from .residues import THREE_LETTER


def chains_in(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            counts[line[21]] = counts.get(line[21], 0) + 1
    return counts


def sequence_of(path: Path, chain: str) -> str:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
            out.append(THREE_LETTER.get(line[17:20].strip(), "X"))
    return "".join(out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb", type=Path)
    ap.add_argument("--chain", help="the chain you intend to design against")
    ap.add_argument("--partner", help="the chain whose contacts define the site")
    ap.add_argument("--hotspots", default="",
                    help="comma-separated residue numbers, e.g. 47,155,157,158")
    args = ap.parse_args(argv)

    cfg = load()
    gk = Gatekeeper(cfg)
    questions: list[str] = []

    print(f"\n=== {args.pdb}")
    print(f"judged by thresholds {cfg.version} ({cfg.hash[:19]}...)\n")

    # --- which chain -------------------------------------------------------
    counts = chains_in(args.pdb)
    if not counts:
        print("no CA atoms at all - is this a protein coordinate file?")
        return 1
    print("chains present (CA count):")
    for ch, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        mark = " <- selected" if ch == args.chain else ""
        print(f"  {ch}: {n:4d} residues{mark}")
    if not args.chain:
        longest = max(counts, key=counts.get)
        questions.append(
            f"which chain is the target? The longest is {longest}, but length is "
            "not identity: 4ZQK's longest chain is PD-L1, and selecting it by "
            "length named it 'pd1' through the whole pipeline."
        )
        print("\nno --chain given, nothing selected")
        _report(questions)
        return 1
    if args.chain not in counts:
        print(f"\nchain {args.chain} is not in this file")
        return 1

    seq = sequence_of(args.pdb, args.chain)
    print(f"\nchain {args.chain} sequence ({len(seq)} residues) - confirm this is "
          "the protein you mean:")
    for i in range(0, len(seq), 60):
        print(f"  {i + 1:4d}  {seq[i:i + 60]}")

    # --- breaks versus numbering ------------------------------------------
    walk = checks.chain_walk(args.pdb, args.chain,
                             break_cutoff=float(cfg.structure["break_cutoff_angstrom"]))
    print(f"\nresidues {walk.first}-{walk.last}, {walk.n_residues} with a CA")
    if walk.numbering_gaps:
        print("  numbering gaps (NOT breaks, the contig stays whole):")
        for a, b, d in walk.numbering_gaps:
            print(f"    {a}->{b}  CA-CA {d:.2f} A")
    if walk.real_breaks:
        print("  real chain breaks (the contig splits here):")
        for a, b, d in walk.real_breaks:
            print(f"    {a}->{b}  CA-CA {d:.2f} A")
    else:
        print("  no real chain breaks")
    print(f"\n  contig:  {walk.contig}")

    # --- what would actually load -----------------------------------------
    cov = checks.atom_coverage(args.pdb, args.chain)
    print(f"\natom names: {cov.summary()}")
    for label, items in (("would be DROPPED at load", cov.dropped),
                         ("no CA, would VANISH at load", cov.no_ca),
                         ("non-standard label, raises at load", cov.nonstandard)):
        if items:
            print(f"  {label}:")
            for item in items[:10]:
                print(f"    {item}")
            if len(items) > 10:
                print(f"    ... and {len(items) - 10} more")
            questions.append(f"{len(items)} residue(s) {label} - fix the file, "
                             "do not let the loader drop them silently")

    # --- chirality ---------------------------------------------------------
    chi = checks.chirality(args.pdb)
    j = gk.chirality(str(args.pdb), chi.n_d, chi.n_checked, chi.d_residues)
    print(f"\nchirality: {j.verdict.value.upper()} - {j.reason}")
    for r in chi.d_residues[:10]:
        print(f"    {r}")
    if j.verdict is not Verdict.PASS:
        questions.append(
            "the chirality check did not pass. Nothing downstream can see this, "
            "and rebuilding the loop is what causes it - do not 'fix' it by "
            "rebuilding again."
        )

    # --- partner and hotspots ---------------------------------------------
    others = {c: n for c, n in counts.items() if c != args.chain}
    print("\nsite definition:")
    if not others:
        print("  no other chain in this file, so hotspots cannot be derived from a partner")
        questions.append("no bound partner in the file - the hotspots have to come "
                         "from the literature, by hand, and be verified against the sequence")
    elif not args.partner:
        print(f"  candidate partner chains: {', '.join(sorted(others))}")
        questions.append(
            "which chain defines the site? Do not let contact count decide: on "
            "PCSK9 3BPS the highest-contact chain is the prodomain, while chain "
            "E (EGF-A of LDLR) gives the documented groove."
        )
    else:
        print(f"  partner chain {args.partner} "
              f"({'present' if args.partner in others else 'NOT IN FILE'})")
        if args.partner not in others:
            questions.append(f"partner chain {args.partner} is not in this file")

    if args.hotspots:
        wanted = [h.strip() for h in args.hotspots.split(",") if h.strip()]
        names = {}
        for line in args.pdb.read_text().splitlines():
            if line.startswith("ATOM") and line[21] == args.chain:
                names[line[22:26].strip()] = line[17:20].strip()
        print(f"\nhotspots ({len(wanted)}), read back out of the file:")
        for h in wanted:
            label = names.get(h)
            print(f"  {args.chain}{h}: {label or 'NOT PRESENT IN THIS CHAIN'}")
            if label is None:
                questions.append(f"hotspot {args.chain}{h} is not in chain {args.chain}")
        print("  confirm each residue name is the one you expect - this is the "
              "check that caught W157 being TRP, twice")
        print(f"\n  hotspot count: {len(wanted)}. Gate-1 pass rate moved 18.3% at "
              "one hotspot to 31.3% at five, so this number changes the "
              "statistics you will compare runs on. It aims a run; it does not "
              "make a bad site work.")
    else:
        questions.append("no hotspots given - set them deliberately, and record "
                         "how many, because the count changes the pass rate")

    _report(questions)
    return 1 if questions else 0


def _report(questions: list[str]) -> None:
    print()
    if not questions:
        print("nothing needs a human answer: this target is ready to submit.")
        return
    print(f"--- {len(questions)} thing(s) need your answer before any GPU time ---")
    for i, q in enumerate(questions, 1):
        print(f"{i}. {q}")
    print("\nNone of these can be defaulted safely: each one lets the run finish "
          "normally and produce numbers that look fine.")


if __name__ == "__main__":
    sys.exit(main())
