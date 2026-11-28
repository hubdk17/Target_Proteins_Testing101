#!/usr/bin/env python3
"""
83_holdout_recheck.py -- strengthen the blind holdout and re-audit an existing
corpus against it.

80_fetch_corpus.py only held out the MDM2 p53-binding domain and the p53 TAD
peptide. A p53 DNA-binding-domain or tetramerization complex, or a p63/p73
structure, would slip past that. This widens the held-out reference set and
re-scans a corpus that is already on disk, moving any newly-caught structure out
of raw/ so training never sees it.

Widened reference set (authoritative -- fetched from RCSB on the login node, not
hardcoded from memory):
  * MDM2 p53-binding domain      (embedded, from 1YCR; cut 0.30 also catches MDMX)
  * p53 transactivation peptide  (embedded, from 1YCR)
  * p53 DNA-binding domain       (fetched: 1TUP; cut 0.35 also catches p63/p73
                                  DBD, which are ~60% identical to p53 DBD)
  * p53 tetramerization domain   (fetched: 1OLG)

Runs on the LOGIN node (needs rcsb). Moves caught files raw/ -> held_out/.

Usage:
  cd ~/rfd/scripts
  python3 model/83_holdout_recheck.py ~/rfd/corpus_v1
"""

import argparse
import glob
import importlib.util
import os
import shutil
import sys
import urllib.request

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "fetch_corpus", os.path.join(_HERE, "80_fetch_corpus.py"))
fc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fc)

FILE_URL = "https://files.rcsb.org/download/%s.pdb"

# (label, pdb_id, min_overlap, cutoff). Embedded refs reuse 80's sequences.
FETCH_REFS = [
    ("p53_DBD", "1TUP", 60, 0.35),
    ("p53_TET", "1OLG", 20, 0.50),
]


def longest_protein_chain(text):
    chains = fc.chains_from_pdb(text)
    if not chains:
        return ""
    return max(chains.values(), key=len)


def fetch(pid):
    try:
        with urllib.request.urlopen(FILE_URL % pid, timeout=60) as r:
            return r.read().decode("utf-8", "replace")
    except Exception as e:
        sys.stderr.write("WARN: could not fetch %s (%s)\n" % (pid, e))
        return None


def build_refs():
    refs = [("MDM2", fc.MDM2_REF, 40, 0.30),
            ("p53_TAD", fc.P53_REF, 10, 0.60)]
    for label, pid, ov, cut in FETCH_REFS:
        txt = fetch(pid)
        if txt is None:
            sys.stderr.write("  (skipping %s -- coverage reduced)\n" % label)
            continue
        seq = longest_protein_chain(txt)
        if len(seq) < ov:
            sys.stderr.write("  (%s from %s too short: %d)\n"
                             % (label, pid, len(seq)))
            continue
        refs.append((label, seq, ov, cut))
        print("  ref %-8s from %s: %d residues (cut %.2f)"
              % (label, pid, len(seq), cut))
    return refs


def is_held(chain_seqs, refs):
    for ch, s in chain_seqs.items():
        for label, ref, ov, cut in refs:
            idn = fc.max_window_identity(s, ref, min_overlap=ov)
            if idn >= cut:
                return "chain %s %.0f%% id to %s" % (ch, 100 * idn, label)
    return None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("corpus", help="corpus dir (contains raw/)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report only; do not move files")
    args = ap.parse_args(argv)

    raw = os.path.join(args.corpus, "raw")
    held = os.path.join(args.corpus, "held_out")
    if not os.path.isdir(raw):
        raise SystemExit("no raw/ under %s" % args.corpus)

    print("building widened reference set ...")
    refs = build_refs()
    print("references in use: %s" % ", ".join(r[0] for r in refs))

    files = sorted(glob.glob(os.path.join(raw, "*.pdb")))
    print("re-scanning %d structures ...\n" % len(files))
    caught = []
    for p in files:
        why = is_held(fc.chains_from_pdb(open(p).read()), refs)
        if why:
            caught.append((os.path.basename(p), why))

    if not caught:
        print("CLEAN: 0 of %d structures match the widened holdout set." % len(files))
        print("(The original MDM2/p53-TAD guard already kept the corpus blind to "
              "the DBD/TET/paralog families too.)")
        return 0

    print("CAUGHT %d structures the original guard missed:" % len(caught))
    for name, why in caught:
        print("  %-10s %s" % (name, why))
    if args.dry_run:
        print("\n--dry-run: nothing moved.")
        return 0
    os.makedirs(held, exist_ok=True)
    for name, _ in caught:
        shutil.move(os.path.join(raw, name), os.path.join(held, name))
    print("\nmoved %d files raw/ -> held_out/. Re-run featurize on raw/ "
          "(fresh --out) so the table excludes them." % len(caught))
    return 0


if __name__ == "__main__":
    sys.exit(main())
