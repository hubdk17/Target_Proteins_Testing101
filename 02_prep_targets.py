#!/usr/bin/env python3
"""
Prepare RFdiffusion inputs for batch 1 of the multi-target run.

For each target, this does the three things that cannot be guessed:
  1. picks the target chain and says how it chose,
  2. finds REAL chain breaks (CA-CA > 4.5 A) and writes the contig from them,
  3. derives hotspots from the binding partner that is already in the crystal,
     rather than from residue numbers recited from memory.

Every target in this batch was chosen because its PDB entry contains the actual
partner (a protein, a peptide or an inhibitor) sitting in the site we want to
block. So the epitope is measured here, not asserted.

Run on the LOGIN NODE. No GPU, no venv needed, standard library only.
    cd ~/rfd && python3 02_prep_targets.py

Writes:  inputs/<name>_clean.pdb  (target chain, protein only)
         inputs/targets.tsv       (the table the array job reads)

Deliberately noisy: it prints what it chose and why for every target, and
flags anything that smells wrong. Read the output before submitting anything.
"""

import os
import sys
import urllib.request

# One entry per target.
#   want_len       expected residue count of the chain we design against; used to
#                  pick the chain, and a mismatch is reported, not silently taken
#   chain          force the target chain instead of picking by length
#   partner_chain  force which chain defines the epitope. REQUIRED whenever the
#                  entry holds more than one thing touching the target, because
#                  "most contacts wins" picks the biggest neighbour, not the
#                  therapeutically interesting one
TARGETS = [
    {"name": "pd1", "pdb": "4ZQK", "want_len": 115,
     "partner_desc": "the other half of the PD-1/PD-L1 pair",
     "site": "PD-1/PD-L1 interface"},

    # Chain P is PCSK9's own prodomain and chain E is the EGF-A domain of LDLR.
    # The prodomain contacts more of the catalytic domain, so picking by contact
    # count aims the run at an intramolecular interface instead of the drug site.
    {"name": "pcsk9", "pdb": "3BPS", "want_len": 435, "partner_chain": "E",
     "partner_desc": "EGF-A domain of LDLR", "site": "EGF-A groove"},

    # Chain C is the BAX BH3 peptide; B and D are a second copy in the cell.
    {"name": "bcl2", "pdb": "2XA0", "want_len": 137, "partner_chain": "C",
     "partner_desc": "BAX BH3 peptide", "site": "BH3 groove"},

    {"name": "brca1", "pdb": "1T29", "want_len": 211, "partner_chain": "B",
     "partner_desc": "phosphorylated BACH1 peptide",
     "site": "BRCT phosphopeptide site"},

    {"name": "glp1r", "pdb": "3IOL", "want_len": 100, "partner_chain": "B",
     "partner_desc": "GLP-1 / exendin-4 peptide",
     "site": "peptide-binding face of the ECD"},

    # 2J30 is PROcaspase-3: one uncleaved 245-residue chain, not the mature
    # p17/p12 pair, and the 174->180 break is the activation loop.
    {"name": "casp3", "pdb": "2J30", "want_len": 245, "partner_chain": "B",
     "partner_desc": "active-site inhibitor", "site": "substrate groove"},
]

AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL",
}
# things that are solvent or buffer, not a binding partner
JUNK = {
    "HOH","DOD","SO4","PO4","CL","NA","K","MG","CA","ZN","MN","GOL","EDO",
    "PEG","PGE","ACT","ACY","DMS","TRS","MPD","NO3","IOD","BR","FMT","CIT",
}

CONTACT_CUTOFF = 5.0   # A, heavy atom to heavy atom
BREAK_CUTOFF   = 4.5   # A, consecutive CA-CA; a real chain break
N_HOTSPOTS     = 5
LEN_TOLERANCE  = 0.30  # warn if the chosen chain is this far off the expectation


def fetch(pdb_id, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return "cached"
    url = "https://files.rcsb.org/download/%s.pdb" % pdb_id
    urllib.request.urlretrieve(url, path)
    return "downloaded"


def parse(path):
    """Return (title, protein_atoms, hetatms). Atoms are dicts."""
    title, prot, het = [], [], []
    with open(path) as fh:
        for line in fh:
            rec = line[:6]
            if rec == "TITLE ":
                title.append(line[10:80].rstrip())
                continue
            if rec not in ("ATOM  ", "HETATM"):
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):          # keep one conformer only
                continue
            element = line[76:78].strip()
            if element == "H" or line[12:16].strip().startswith("H"):
                continue                          # heavy atoms only
            resname = line[17:20].strip()
            if resname in JUNK:
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            atom = {
                "line": line,
                "chain": line[21],
                "resnum": line[22:26].strip(),
                "icode": line[26],
                "resname": resname,
                "name": line[12:16].strip(),
                "xyz": xyz,
            }
            if rec == "ATOM  " and resname in AA3:
                prot.append(atom)
            else:
                het.append(atom)
    return " ".join(title), prot, het


def chain_molecules(path):
    """{chain letter: molecule name} from the COMPND records.

    This is what says which protein a chain actually is. Two chains of similar
    length (PD-1 and PD-L1 are both ~110 residue Ig domains) cannot be told
    apart by size, and picking the wrong one aims the whole run at the wrong
    protein while every other number still looks reasonable.
    """
    out, molecule = {}, None
    pending = ""
    with open(path) as fh:
        for line in fh:
            if line[:6] != "COMPND":
                if out or molecule:
                    break
                continue
            text = line[10:80].strip()
            pending = (pending + " " + text).strip() if pending else text
            if not pending.endswith(";"):
                continue
            field = pending.rstrip(";")
            pending = ""
            if field.startswith("MOLECULE:"):
                molecule = field[len("MOLECULE:"):].strip()
            elif field.startswith("CHAIN:") and molecule:
                for c in field[len("CHAIN:"):].split(","):
                    c = c.strip()
                    if c:
                        out[c] = molecule
    return out


def chain_lengths(prot):
    seen, order = {}, {}
    for a in prot:
        key = (a["chain"], a["resnum"], a["icode"])
        if key in seen:
            continue
        seen[key] = True
        order[a["chain"]] = order.get(a["chain"], 0) + 1
    return order


def dist2(p, q):
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    dz = p[2] - q[2]
    return dx * dx + dy * dy + dz * dz


def contig_and_breaks(cas):
    """cas: list of (resnum_int, xyz) in file order. Returns (contig, breaks)."""
    segs, breaks = [], []
    start = prev = cas[0][0]
    prev_xyz = cas[0][1]
    for num, xyz in cas[1:]:
        if dist2(xyz, prev_xyz) > BREAK_CUTOFF ** 2:
            segs.append((start, prev))
            breaks.append((prev, num))
            start = num
        prev, prev_xyz = num, xyz
    segs.append((start, prev))
    return segs, breaks


def main():
    root = os.path.expanduser("~/rfd")
    indir = os.path.join(root, "inputs")
    os.makedirs(indir, exist_ok=True)
    rows, problems = [], []

    for spec in TARGETS:
        name = spec["name"]
        pdb_id = spec["pdb"]
        want_len = spec["want_len"]
        partner_desc = spec["partner_desc"]
        site_desc = spec["site"]
        force_chain = spec.get("chain")
        force_partner = spec.get("partner_chain")
        raw = os.path.join(indir, "%s.pdb" % pdb_id)
        print("=" * 72)
        try:
            how = fetch(pdb_id, raw)
        except Exception as exc:
            print("%-7s %s  FETCH FAILED: %s" % (name, pdb_id, exc))
            problems.append("%s: could not download %s" % (name, pdb_id))
            continue

        title, prot, het = parse(raw)
        print("%-7s %s (%s)" % (name, pdb_id, how))
        print("  TITLE: %s" % title[:200])
        if not prot:
            problems.append("%s: no protein atoms parsed" % name)
            print("  NO PROTEIN ATOMS -- wrong file?")
            continue

        lens = chain_lengths(prot)
        mols = chain_molecules(raw)
        print("  chains: %s" % ", ".join(
            "%s=%d" % (c, n) for c, n in sorted(lens.items())))
        for c in sorted(lens):
            print("    chain %s  %d res  %s" % (c, lens[c], mols.get(c, "(unnamed)")))

        # the named chain if the entry forces one, else the closest to the
        # expected domain length
        if force_chain:
            if force_chain not in lens:
                problems.append("%s: forced chain %s absent (have %s)"
                                % (name, force_chain, sorted(lens)))
                print("  *** FORCED CHAIN %s ABSENT" % force_chain)
                continue
            chain = force_chain
        else:
            chain = min(lens, key=lambda c: abs(lens[c] - want_len))
        got = lens[chain]
        print("  chose chain %s = %s" % (chain, mols.get(chain, "(unnamed)")))
        print("    %d residues; expected ~%d%s -- site: %s"
              % (got, want_len,
                 " [forced]" if force_chain else "", site_desc))
        off = abs(got - want_len) / float(want_len)
        if off > LEN_TOLERANCE:
            msg = ("%s: chain %s has %d residues, expected ~%d (%.0f%% off)"
                   % (name, chain, got, want_len, 100 * off))
            problems.append(msg)
            print("  *** LENGTH MISMATCH -- check this is the right chain/entry")

        tgt = [a for a in prot if a["chain"] == chain]

        icodes = sorted({a["icode"] for a in tgt} - {" "})
        if icodes:
            problems.append("%s: insertion codes %s in chain %s"
                            % (name, icodes, chain))
            print("  *** INSERTION CODES %s -- a numeric contig will be wrong"
                  % icodes)

        # contig from real chain breaks
        cas, seen = [], set()
        for a in tgt:
            if a["name"] != "CA" or a["resnum"] in seen:
                continue
            seen.add(a["resnum"])
            cas.append((int(a["resnum"]), a["xyz"]))
        segs, breaks = contig_and_breaks(cas)
        contig = "/".join("%s%d-%d" % (chain, s, e) for s, e in segs)
        print("  CAs: %d   span: %d-%d" % (len(cas), cas[0][0], cas[-1][0]))
        if breaks:
            print("  REAL chain breaks (CA-CA > %.1f A): %s"
                  % (BREAK_CUTOFF, ", ".join("%d->%d" % b for b in breaks)))
        else:
            print("  no chain breaks")
        print("  contig: %s" % contig)

        # hotspots, measured against the partner in the crystal.
        # Exclude anything on the target's OWN chain: a ligand or modified
        # residue carries the target's chain letter and would otherwise be
        # counted as a binding partner against itself.
        partner = [a for a in prot + het if a["chain"] != chain]
        if not partner:
            problems.append("%s: no partner in the file, cannot derive hotspots"
                            % name)
            print("  *** NO PARTNER FOUND -- expected %s" % partner_desc)
            continue
        # An entry can hold several partners (a second copy of the complex, a
        # crystallisation additive). Hotspots have to mark ONE patch, so count
        # against each partner chain separately and keep the best-contacting one.
        cut2 = CONTACT_CUTOFF ** 2
        by_chain = {}
        for a in partner:
            by_chain.setdefault(a["chain"], []).append(a["xyz"])

        per_partner = {}
        for pch, pts in by_chain.items():
            counts_p = {}
            for a in tgt:
                n = sum(1 for q in pts if dist2(a["xyz"], q) <= cut2)
                if n:
                    counts_p[a["resnum"]] = counts_p.get(a["resnum"], 0) + n
            if counts_p:
                per_partner[pch] = counts_p

        if len(per_partner) > 1:
            print("  partners in contact: %s" % ", ".join(
                "%s(%d res)" % (c, len(v)) for c, v in sorted(per_partner.items())))

        counts = {}
        if per_partner:
            if force_partner:
                # Named explicitly, because contact count picks the biggest
                # neighbour and that is not always the site we care about. On
                # 3BPS the prodomain outranks the EGF-A domain we are aiming at.
                if force_partner not in per_partner:
                    problems.append(
                        "%s: forced partner chain %s has no contacts (in contact: %s)"
                        % (name, force_partner, sorted(per_partner)))
                    print("  *** FORCED PARTNER %s NOT IN CONTACT" % force_partner)
                    continue
                best = force_partner
            else:
                best = max(per_partner, key=lambda c: sum(per_partner[c].values()))
                if len(per_partner) > 1:
                    problems.append(
                        "%s: %d things touch the target and no partner_chain is set; "
                        "picked %s by contact count - confirm it is the right site"
                        % (name, len(per_partner), best))
            counts = per_partner[best]
            partner_desc = "%s, chain %s (%d atoms)%s" % (
                partner_desc, best, len(by_chain[best]),
                " [forced]" if force_partner else "")

        if not counts:
            problems.append("%s: partner present but zero contacts" % name)
            print("  *** ZERO CONTACTS -- partner may be a crystal neighbour")
            continue

        ranked = sorted(counts.items(), key=lambda kv: -kv[1])[:N_HOTSPOTS]
        resnames = {a["resnum"]: a["resname"] for a in tgt}
        hotspots = ",".join("%s%s" % (chain, r) for r, _ in ranked)
        print("  partner: %s, %d residues in contact"
              % (partner_desc, len(counts)))
        print("  hotspots: %s" % hotspots)
        for r, n in ranked:
            print("      %s%-5s %s  %d contact atoms"
                  % (chain, r, resnames.get(r, "???"), n))

        # Hotspots must be one patch. If they are spread out, RFdiffusion is
        # being asked to reach two places at once and will do neither well.
        hot_ca = {}
        for a in tgt:
            if a["name"] == "CA":
                hot_ca[a["resnum"]] = a["xyz"]
        pts = [hot_ca[r] for r, _ in ranked if r in hot_ca]
        if len(pts) > 1:
            cx = sum(p[0] for p in pts) / len(pts)
            cy = sum(p[1] for p in pts) / len(pts)
            cz = sum(p[2] for p in pts) / len(pts)
            spread = max(dist2(p, (cx, cy, cz)) for p in pts) ** 0.5
            print("  hotspot patch radius: %.1f A" % spread)
            if spread > 15.0:
                problems.append("%s: hotspots span %.1f A, not one patch"
                                % (name, spread))
                print("  *** HOTSPOTS NOT ONE PATCH -- narrow the site by hand")

        # cleaned input: target chain, protein only
        clean = os.path.join(indir, "%s_clean.pdb" % name)
        with open(clean, "w") as out:
            for a in tgt:
                out.write(a["line"])
            out.write("TER\nEND\n")
        rows.append((name, pdb_id, chain, str(len(cas)), contig, hotspots))

    tsv = os.path.join(indir, "targets.tsv")
    with open(tsv, "w") as out:
        for r in rows:
            out.write("\t".join(r) + "\n")

    print("=" * 72)
    print("wrote %s with %d of %d targets" % (tsv, len(rows), len(TARGETS)))
    for r in rows:
        print("  %-7s %s chain %s  n=%-4s %s  hotspots %s"
              % (r[0], r[1], r[2], r[3], r[4], r[5]))
    if problems:
        print("\n%d THING(S) TO LOOK AT BEFORE SUBMITTING:" % len(problems))
        for p in problems:
            print("  - %s" % p)
    else:
        print("\nno complaints from any target")
    return 0


if __name__ == "__main__":
    sys.exit(main())
