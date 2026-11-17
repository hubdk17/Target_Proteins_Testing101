#!/usr/bin/env python3
"""
80_fetch_corpus.py -- build a BLIND training corpus of protein-protein
interface complexes from the RCSB PDB.

Why this exists
---------------
We want to train our own model (OncoFlow-0: a learned, target-conditioned
model of interface chemistry) instead of relying on RFdiffusion / ProteinMPNN
/ AF2, each of which was trained on the whole PDB -- including 1YCR and many
other MDM2 complexes -- and so cannot be blind to our test target.

A model is only as blind as its training set. If a single MDM2 / MDMX / p53
complex leaks into the corpus, then "our model designed an MDM2 binder" is a
check that cannot fail: it memorised the answer. So the ONE job of this script
is to fetch a corpus and PROVE that the held-out structures are excluded.

What it does (runs on the cluster LOGIN node, which reaches rcsb.org;
the cloud container's egress blocks files.rcsb.org)
----------------------------------------------------------------------
1. Query the RCSB search API for entries with >= 2 protein entities and
   resolution <= --max-res (real protein-protein interfaces, decent quality).
2. Download each entry's .pdb from files.rcsb.org.
3. HOLDOUT GUARD: extract every chain's sequence and compare it, gaplessly,
   to the MDM2 p53-binding domain and the p53 transactivation peptide
   (reference sequences below, taken verbatim from 1YCR). Any structure with
   a chain above the identity threshold -- which also catches MDMX, ~50%
   identical to MDM2 over this domain -- is dropped and logged, never saved.
4. Write kept structures to <out>/raw/ and a manifest.tsv recording, for
   EVERY id, whether it was kept or held out and why. Nothing is silent.

Self-test (--selftest): run the guard on the two reference sequences
themselves and assert it flags them. A holdout guard that cannot be shown to
fire on 1YCR is not a guard. Exit non-zero if it fails to flag them.

Usage
-----
  module load slurm          # (not needed here; login node only)
  cd ~/rfd/scripts
  # validate first on a small pull, and confirm the self-test passes:
  python3 model/80_fetch_corpus.py --n 25 --out ~/rfd/corpus --selftest
  # then scale up once the manifest looks right:
  python3 model/80_fetch_corpus.py --n 2000 --out ~/rfd/corpus

This writes NO GitHub, needs no token, and keeps every download under --out.
"""

import argparse
import json
import os
import sys
import time
import urllib.request
import urllib.error

# -- reference sequences to hold out, verbatim from 1YCR (verified) ----------
# chain A, MDM2 25-109 (the p53-binding domain)
MDM2_REF = ("ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDL"
            "FGVPSFSVKEHRKIYTMIYRNLVV")
# chain B, p53 transactivation peptide 17-29
P53_REF = "ETFSDLWKLLPEN"

# belt-and-suspenders: a few well-known MDM2/MDMX/p53 entry ids. The sequence
# guard is the real safety; this just short-circuits the obvious ones.
BLOCK_IDS = {"1YCR", "1T4E", "1RV1", "3EQS", "4HG7", "4ZYI", "4ZYF",
             "2MPS", "3LBK", "1Z1M", "3DAC", "4OAS", "5AFG"}

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
FILE_URL = "https://files.rcsb.org/download/%s.pdb"


def max_window_identity(seq, ref, min_overlap):
    """Best gapless identity of `ref` against `seq`, over all offsets.

    Slides the shorter sequence across the longer one with no gaps and returns
    the highest fraction of matching residues over any overlap of length
    >= min_overlap. Gapless on purpose: it is crude but conservative and needs
    no alignment library on the login node. It can miss a homolog that differs
    by an insertion/deletion, so the MDM2 threshold is set generously low to
    absorb that -- better to drop a borderline structure than leak one.
    """
    if not seq or not ref:
        return 0.0
    a, b = (seq, ref) if len(seq) >= len(ref) else (ref, seq)
    best = 0.0
    for off in range(-len(b) + 1, len(a)):
        match = tot = 0
        for i in range(len(b)):
            j = off + i
            if 0 <= j < len(a):
                tot += 1
                if b[i] == a[j]:
                    match += 1
        if tot >= min_overlap:
            best = max(best, match / tot)
    return best


def is_holdout(chain_seqs, mdm2_cut, p53_cut):
    """Return (True, reason) if any chain is too close to a held-out target."""
    for ch, s in chain_seqs.items():
        im = max_window_identity(s, MDM2_REF, min_overlap=40)
        if im >= mdm2_cut:
            return True, "chain %s %.0f%% id to MDM2/MDMX domain" % (
                ch, 100 * im)
        ip = max_window_identity(s, P53_REF, min_overlap=10)
        if ip >= p53_cut:
            return True, "chain %s %.0f%% id to p53 peptide" % (ch, 100 * ip)
    return False, ""


def chains_from_pdb(text):
    """Map chain id -> one-letter sequence, first model, protein atoms only."""
    seq = {}
    seen = set()
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"):
            break
        if not ln.startswith("ATOM"):
            continue
        rn = ln[17:20].strip()
        if rn not in THREE2ONE:
            continue
        ch = ln[21]
        key = (ch, ln[22:27])            # resSeq + iCode, de-dupe by residue
        if key in seen:
            continue
        seen.add(key)
        seq.setdefault(ch, []).append(THREE2ONE[rn])
    return {c: "".join(v) for c, v in seq.items()}


def rcsb_search(n, max_res, page=10000):
    """Return up to n PDB ids. RCSB caps a single request at 10,000 rows and
    returns HTTP 400 above that, so page through the result set when n is larger.
    """
    def _query(start, rows):
        return {
            "query": {"type": "group", "logical_operator": "and", "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                    "operator": "greater_or_equal", "value": 2}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal", "value": max_res}},
            ]},
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": start, "rows": rows},
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "score", "direction": "desc"}]},
        }

    ids, start = [], 0
    while len(ids) < n:
        rows = min(page, n - len(ids))
        data = json.dumps(_query(start, rows)).encode()
        req = urllib.request.Request(
            SEARCH_URL, data=data,
            headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            res = json.loads(r.read().decode())
        batch = [row["identifier"] for row in res.get("result_set", [])]
        if not batch:
            break                       # ran out of matching entries
        ids.extend(batch)
        print("  search: %d ids so far" % len(ids))
        start += len(batch)
        if len(batch) < rows:
            break                       # last (partial) page
    return ids[:n]


def fetch(pdb_id, tries=3):
    for k in range(tries):
        try:
            with urllib.request.urlopen(FILE_URL % pdb_id, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None            # no .pdb (likely cif-only); skip
            time.sleep(1.5 * (k + 1))
        except Exception:
            time.sleep(1.5 * (k + 1))
    return None


def selftest(mdm2_cut, p53_cut):
    ok = True
    hold, why = is_holdout({"A": MDM2_REF}, mdm2_cut, p53_cut)
    print("[selftest] MDM2 ref -> holdout=%s (%s)" % (hold, why))
    ok = ok and hold
    hold, why = is_holdout({"B": P53_REF}, mdm2_cut, p53_cut)
    print("[selftest] p53  ref -> holdout=%s (%s)" % (hold, why))
    ok = ok and hold
    # a clearly unrelated sequence must NOT be held out (guard can be wrong
    # in the other direction too)
    decoy = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKA"
    hold, why = is_holdout({"A": decoy}, mdm2_cut, p53_cut)
    print("[selftest] decoy   -> holdout=%s (expect False)" % hold)
    ok = ok and (not hold)
    print("[selftest] %s" % ("PASS" if ok else "FAIL"))
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=25,
                    help="number of entries to pull (start small to validate)")
    ap.add_argument("--max-res", type=float, default=2.5)
    ap.add_argument("--out", default="corpus")
    ap.add_argument("--mdm2-cut", type=float, default=0.30,
                    help="identity to MDM2 domain at/above which to hold out "
                         "(0.30 also catches MDMX)")
    ap.add_argument("--p53-cut", type=float, default=0.60)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)

    if not selftest(args.mdm2_cut, args.p53_cut) and args.selftest:
        raise SystemExit("HOLDOUT GUARD SELF-TEST FAILED -- not fetching.")
    if args.selftest:
        print("[selftest] guard verified; proceeding to fetch.\n")

    raw = os.path.join(args.out, "raw")
    os.makedirs(raw, exist_ok=True)
    man = os.path.join(args.out, "manifest.tsv")
    if os.path.exists(man):
        raise SystemExit("%s exists; pass a fresh --out (no silent append)."
                         % man)

    print("searching RCSB for <= %g A, >=2 protein chains ..." % args.max_res)
    ids = rcsb_search(args.n, args.max_res)
    print("got %d entry ids" % len(ids))

    cols = ["pdb_id", "status", "n_chains", "reason"]
    kept = held = skipped = 0
    with open(man, "w") as mh:
        mh.write("\t".join(cols) + "\n")
        for i, pid in enumerate(ids, 1):
            if pid.upper() in BLOCK_IDS:
                mh.write("%s\tholdout\tNA\tblocklisted id\n" % pid)
                held += 1
                continue
            text = fetch(pid)
            if text is None:
                mh.write("%s\tskip\tNA\tno .pdb (cif-only or 404)\n" % pid)
                skipped += 1
                continue
            cseq = chains_from_pdb(text)
            if len(cseq) < 2:
                mh.write("%s\tskip\t%d\t<2 protein chains in file\n"
                         % (pid, len(cseq)))
                skipped += 1
                continue
            hold, why = is_holdout(cseq, args.mdm2_cut, args.p53_cut)
            if hold:
                mh.write("%s\tholdout\t%d\t%s\n" % (pid, len(cseq), why))
                held += 1
                continue
            with open(os.path.join(raw, "%s.pdb" % pid), "w") as fh:
                fh.write(text)
            mh.write("%s\tkept\t%d\t-\n" % (pid, len(cseq)))
            kept += 1
            if i % 10 == 0:
                print("  %d/%d  kept=%d held=%d skip=%d"
                      % (i, len(ids), kept, held, skipped))
            time.sleep(0.2)

    print("\nDONE  kept=%d  held_out=%d  skipped=%d" % (kept, held, skipped))
    print("corpus structures : %s/" % raw)
    print("manifest          : %s" % man)
    print("Sanity check the held-out rows before training -- every MDM2/MDMX/"
          "p53 hit must appear there, never in raw/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
