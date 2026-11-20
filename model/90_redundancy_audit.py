#!/usr/bin/env python3
"""
90_redundancy_audit.py -- how redundant is the corpus, and does that redundancy
leak between train and test?

Why this matters
----------------
Our recovery numbers (82/84/86) split BY STRUCTURE, which stops the SAME file
appearing on both sides. It does NOT stop two DIFFERENT files that are the same
protein (the PDB has lysozyme hundreds of times, every popular target many
times). If a near-identical chain sits in both train and test, the model can
score well by memorising a homolog -- a "check that cannot fail". This measures
how exposed we are to that, so we know how much to trust 0.267.

What it reports (stdlib only -- runs on the bare login node, no numpy):
  * number of structures and of UNIQUE longest-chain sequences
  * redundancy exposure: % of structures whose longest chain is byte-identical to
    some OTHER structure's (these straddle train/test under any random split)
  * the biggest exact-duplicate groups

This is EXACT-match redundancy = a LOWER bound. Homolog-level redundancy (a few
mutations apart) is higher and needs real clustering (mmseqs/cd-hit, not on this
cluster) -- if the exact number is already high, that is the next build.

Usage:
  python3 model/90_redundancy_audit.py ~/rfd/corpus_v1/raw
  python3 model/90_redundancy_audit.py ~/rfd/corpus_v2/raw
"""

import argparse
import glob
import importlib.util
import os
import sys
from collections import defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "fetch_corpus", os.path.join(_HERE, "80_fetch_corpus.py"))
fc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fc)


def longest_chain(text):
    chains = fc.chains_from_pdb(text)
    return max(chains.values(), key=len) if chains else ""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", help="corpus raw/ directory of .pdb files")
    ap.add_argument("--min-len", type=int, default=10,
                    help="ignore chains shorter than this")
    ap.add_argument("--emit-keep", default=None,
                    help="write one representative structure per unique "
                         "longest-chain sequence (a non-redundant keep-list for "
                         "86's --keep-list). Removes EXACT duplicates.")
    args = ap.parse_args(argv)

    files = sorted(glob.glob(os.path.join(args.raw, "*.pdb")))
    if not files:
        raise SystemExit("no .pdb files in %s" % args.raw)

    by_seq = defaultdict(list)
    n_ok = 0
    for p in files:
        try:
            s = longest_chain(open(p).read())
        except Exception:
            continue
        if len(s) < args.min_len:
            continue
        by_seq[s].append(os.path.basename(p))
        n_ok += 1

    n_unique = len(by_seq)
    # structures that share their longest-chain sequence with >=1 other structure
    in_dup = sum(len(v) for v in by_seq.values() if len(v) > 1)
    dup_groups = sorted((v for v in by_seq.values() if len(v) > 1),
                        key=len, reverse=True)

    print("=== redundancy audit: %s ===" % args.raw)
    print("structures with a usable protein chain : %d" % n_ok)
    print("UNIQUE longest-chain sequences         : %d" % n_unique)
    print("redundancy ratio (unique / structures) : %.3f"
          % (n_unique / n_ok if n_ok else 0))
    print("structures sharing a sequence w/ another: %d  (%.1f%% EXPOSED)"
          % (in_dup, 100.0 * in_dup / n_ok if n_ok else 0))
    print("exact-duplicate groups                 : %d" % len(dup_groups))
    print("\nbiggest exact-duplicate groups (copies : one example):")
    for v in dup_groups[:10]:
        print("  %4d copies   e.g. %s" % (len(v), v[0]))

    if args.emit_keep:
        # one representative per unique sequence = the exact-dedup keep set
        reps = sorted(v[0] for v in by_seq.values())
        with open(args.emit_keep, "w") as fh:
            fh.write("\n".join(reps) + "\n")
        print("\nemitted %d non-redundant representatives -> %s"
              % (len(reps), args.emit_keep))
        print("(feed to 86 with --keep-list to train/eval on unique sequences; "
              "collapses by LONGEST chain, so it is conservative.)")

    print("\nRead: 'EXPOSED %' is the share of the corpus that has an identical "
          "twin elsewhere; under a random by-structure split most of those twins "
          "land across train and test, so that much of the recovery number is "
          "suspect. This is EXACT match only (lower bound); homologs push it "
          "higher. If this is large, the fix is a sequence-CLUSTER split, not a "
          "by-structure split.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
