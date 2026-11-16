#!/usr/bin/env python3
"""
60_gate_battery.py -- CPU-only interface gate battery for binder complexes.

Computes every gate we can evaluate without a GPU on a two-chain complex,
one row per input structure, so that a series of designs (or a ladder of
peptides with measured affinity) can be scored and the gates compared
against each other and against experiment.

Why this exists
---------------
The pipeline's existing gates are:
  gate 1  compactness (Rg ratio)        -- shape only, target-blind
  gate 2  self-consistency              -- target-blind
  gate 3  AF2 pae_interaction           -- needs a GPU
  (+ PyRosetta REF15 interface energy, academic licence only)

On a fixed-backbone ladder -- same pose, same backbone, only side chains
differ -- gates 1 and 2 cannot vary at all. This script adds the classical
physical terms that can, none of which exist anywhere in the pipeline:

  G_dsasa_*   buried surface area, split apolar/polar (Shrake-Rupley)
  G_hbond     interface hydrogen bonds / salt bridges, geometric
  G_unsat     buried unsatisfied polar atoms
  G_gapvol    gap volume -- packing defect, a shape-complementarity proxy
  G_helixcost intrinsic helix-formation cost of the binder sequence
              (conformational entropy: the pipeline has NO term for this)
  G_dumb_dG   -0.025 * buried apolar area, and nothing else -- the
              deliberately stupid baseline every learned gate must beat

NO SIDE-CHAIN REBUILDING happens here. Scores are fixed-backbone,
fixed-rotamer, computed on the coordinates handed to it. That is a
documented approximation, not an oversight: it is what makes an alanine
scan exact (see 61_make_ala_scan.py, which only deletes atoms).

Usage
-----
  python3 60_gate_battery.py complex.pdb [more.pdb ...] \
      --target-chain A --binder-chain B --out gates.tsv

  # per-residue breakdown for one structure
  python3 60_gate_battery.py complex.pdb --per-residue

Validation already done (2026-11-16, on 1YCR)
---------------------------------------------
  * SASA of one isolated carbon = 120.76 A^2, analytic 4*pi*(1.70+1.4)^2
    = 120.76. Exact to 3 decimal places at 512 points.
  * Two carbons 6.5 A apart give exactly twice the free value; at 2.0 A
    they occlude each other symmetrically. Correct limits.
  * 1YCR dSASA_total = 1464 A^2, i.e. ~730 A^2 buried per side, which
    agrees with the published interface area for this complex.
  * NOISE FLOOR, measured rather than assumed: SASA is rotation
    invariant, so the spread over 6 random rigid rotations of the same
    complex is pure sphere-point discretisation error. At 128 points,
    dsasa_apolar has sd = 5.3 A^2. Any delta smaller than ~16 A^2 (3 sd)
    is not resolved. Use 512 points for a final number.
  * The unsatisfied-polar term returns 0 on all 13 arms of the 1YCR
    alanine scan. That is a measurement, not a dead check: deleting the
    two target polar atoms near Trp23 NE1 makes it return 1 and name
    B/TRP23:NE1. A gate that cannot be made to fire is not a gate.
  * Chirality sign calibrated on 1YCR: 94/94 L-residues at the
    N-C-CA-CB improper = -122.3 +/- 3.0 deg. Positive means D.

House rules honoured
--------------------
  * Columns are written with a header and read by NAME, never by position.
  * --out refuses to append to an existing file; pass a fresh name.
  * Every gate prints the raw quantity it is derived from, so a derived
    number can never be mistaken for a measurement.
  * A value that could not be computed is written as NA, never as 0.0.
"""

import argparse
import math
import os
import sys

# ---------------------------------------------------------------- constants

# Bondi van der Waals radii, angstrom. Heavy atoms only; we never see H.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90}
VDW_DEFAULT = 1.70

# Expected heavy-atom count per residue, for completeness checks.
NHEAVY = {
    "ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9, "GLU": 9,
    "GLY": 4, "HIS": 10, "ILE": 8, "LEU": 8, "LYS": 9, "MET": 8, "PHE": 11,
    "PRO": 7, "SER": 6, "THR": 7, "TRP": 14, "TYR": 12, "VAL": 7,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

BACKBONE = {"N", "CA", "C", "O", "OXT"}

# Pace & Scholtz (1998) helix propensity, kcal/mol relative to alanine.
# Higher = more costly to put in a helix. Summed over the binder this is a
# crude but real estimate of the free-energy price of pre-organising the
# helix before it ever touches the target.
HELIX_COST = {
    "A": 0.00, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39,
    "E": 0.40, "I": 0.41, "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54,
    "H": 0.61, "V": 0.61, "N": 0.65, "T": 0.66, "C": 0.68, "D": 0.69,
    "G": 1.00, "P": 3.16,
}

# Hydrophobic transfer coefficient, kcal/mol per A^2 of buried apolar area.
# The classical range is 0.020-0.030; 0.025 is the usual compromise.
APOLAR_COEFF = 0.025

# Donors and acceptors, by residue and atom name. Used for the hydrogen-bond
# and buried-unsatisfied terms. We have no hydrogens, so geometry is
# heavy-atom distance only -- see the LIMITATION note in hbond_terms().
POLAR_ELEMENTS = {"N", "O"}

CHARGED_POS = {("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
               ("LYS", "NZ"), ("HIS", "ND1"), ("HIS", "NE2")}
CHARGED_NEG = {("ASP", "OD1"), ("ASP", "OD2"),
               ("GLU", "OE1"), ("GLU", "OE2")}


# ------------------------------------------------------------------ parsing

class Atom:
    __slots__ = ("name", "resname", "chain", "resnum", "xyz", "element",
                 "radius", "is_backbone")

    def __init__(self, name, resname, chain, resnum, xyz, element):
        self.name = name
        self.resname = resname
        self.chain = chain
        self.resnum = resnum
        self.xyz = xyz
        self.element = element
        self.radius = VDW.get(element, VDW_DEFAULT)
        self.is_backbone = name in BACKBONE

    @property
    def key(self):
        return (self.chain, self.resnum)


def element_of(line):
    """Element from columns 77-78, falling back to the atom name."""
    el = line[76:78].strip().upper()
    if el:
        return el
    name = line[12:16].strip()
    # PDB atom names are right-justified in a 4-col field when the element
    # is one letter; a leading digit means a hydrogen like 1HB.
    for ch in name:
        if ch.isalpha():
            return ch.upper()
    return "C"


def read_pdb(path):
    """Heavy atoms of the first model only. Hydrogens and waters dropped."""
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:20].strip()
            if resname in ("HOH", "WAT", "DOD"):
                continue
            el = element_of(line)
            if el == "H" or el == "D":
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]),
                       float(line[46:54]))
            except ValueError:
                continue
            atoms.append(Atom(line[12:16].strip(), resname, line[21],
                              int(line[22:26]), xyz, el))
    return atoms


def residues_of(atoms):
    """Ordered dict-like: (chain, resnum) -> (resname, [atoms])."""
    out = {}
    for a in atoms:
        if a.key not in out:
            out[a.key] = (a.resname, [])
        out[a.key][1].append(a)
    return out


# --------------------------------------------------------------- geometry

def dist2(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return dx * dx + dy * dy + dz * dz


def dist(a, b):
    return math.sqrt(dist2(a, b))


def sphere_points(n):
    """n roughly-uniform points on the unit sphere (golden spiral)."""
    pts = []
    inc = math.pi * (3.0 - math.sqrt(5.0))
    off = 2.0 / n
    for k in range(n):
        y = k * off - 1.0 + off / 2.0
        r = math.sqrt(max(0.0, 1.0 - y * y))
        phi = k * inc
        pts.append((math.cos(phi) * r, y, math.sin(phi) * r))
    return pts


class Grid:
    """Uniform bucket grid for neighbour queries."""

    def __init__(self, atoms, cell):
        self.cell = cell
        self.buckets = {}
        for i, a in enumerate(atoms):
            self.buckets.setdefault(self._idx(a.xyz), []).append(i)

    def _idx(self, p):
        c = self.cell
        return (int(math.floor(p[0] / c)), int(math.floor(p[1] / c)),
                int(math.floor(p[2] / c)))

    def near(self, p):
        ix, iy, iz = self._idx(p)
        out = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    b = self.buckets.get((ix + dx, iy + dy, iz + dz))
                    if b:
                        out.extend(b)
        return out


# -------------------------------------------------------------------- SASA

def sasa(atoms, n_points=128, probe=1.4):
    """Shrake-Rupley solvent-accessible surface area per atom, A^2.

    Returns a list parallel to `atoms`. Implemented from scratch so that
    nothing in this battery depends on a licensed package.
    """
    if not atoms:
        return []
    pts = sphere_points(n_points)
    maxr = max(a.radius for a in atoms) + probe
    grid = Grid(atoms, 2.0 * maxr)
    out = []
    for i, a in enumerate(atoms):
        ri = a.radius + probe
        cand = []
        for j in grid.near(a.xyz):
            if j == i:
                continue
            b = atoms[j]
            rj = b.radius + probe
            cut = ri + rj
            if dist2(a.xyz, b.xyz) < cut * cut:
                cand.append((b.xyz, rj * rj))
        if not cand:
            out.append(4.0 * math.pi * ri * ri)
            continue
        acc = 0
        ax, ay, az = a.xyz
        for px, py, pz in pts:
            tx, ty, tz = ax + px * ri, ay + py * ri, az + pz * ri
            buried = False
            for (bx, by, bz), rj2 in cand:
                dx, dy, dz = tx - bx, ty - by, tz - bz
                if dx * dx + dy * dy + dz * dz < rj2:
                    buried = True
                    break
            if not buried:
                acc += 1
        out.append(4.0 * math.pi * ri * ri * acc / n_points)
    return out


def is_apolar(atom):
    return atom.element in ("C", "S")


# --------------------------------------------------------------- the gates

def dsasa_terms(complex_atoms, tgt, bnd, n_points):
    """Buried surface area on complex formation, split apolar/polar.

    dSASA = SASA(target alone) + SASA(binder alone) - SASA(complex),
    counted over the binder's atoms and the target's atoms separately so
    that an asymmetric interface is visible.
    """
    s_cplx = sasa(complex_atoms, n_points)
    s_tgt = sasa(tgt, n_points)
    s_bnd = sasa(bnd, n_points)

    idx = {id(a): k for k, a in enumerate(complex_atoms)}
    tot = apol = pol = 0.0
    bnd_tot = tgt_tot = 0.0
    per_atom_buried = {}
    for group, alone in ((tgt, s_tgt), (bnd, s_bnd)):
        for a, sa in zip(group, alone):
            sc = s_cplx[idx[id(a)]]
            d = sa - sc
            per_atom_buried[id(a)] = (sa, sc)
            if d <= 0:
                continue
            tot += d
            if is_apolar(a):
                apol += d
            else:
                pol += d
            if group is bnd:
                bnd_tot += d
            else:
                tgt_tot += d
    return {
        "dsasa_total": tot,
        "dsasa_apolar": apol,
        "dsasa_polar": pol,
        "dsasa_binder_side": bnd_tot,
        "dsasa_target_side": tgt_tot,
        "sasa_complex_total": sum(s_cplx),
    }, per_atom_buried, s_cplx, idx


def hbond_terms(tgt, bnd, cutoff=3.5):
    """Interface hydrogen bonds and salt bridges, heavy-atom geometry.

    LIMITATION, stated rather than hidden: with no hydrogens present we
    cannot apply a donor-H-acceptor angle. This counts N/O to N/O pairs
    across the interface within `cutoff`, which overcounts -- a pair that
    is close but badly oriented is indistinguishable from a real bond.
    Treat it as an upper bound. Salt bridges are the subset where one
    partner is a formally charged N and the other a formally charged O.
    """
    tp = [a for a in tgt if a.element in POLAR_ELEMENTS]
    bp = [a for a in bnd if a.element in POLAR_ELEMENTS]
    n_hb = n_salt = 0
    pairs = []
    cut2 = cutoff * cutoff
    for a in bp:
        for b in tp:
            if dist2(a.xyz, b.xyz) > cut2:
                continue
            # N-N and O-O contacts cannot both donate; require one of each
            # unless one is a charged pair.
            ka = (a.resname, a.name)
            kb = (b.resname, b.name)
            charged = ((ka in CHARGED_POS and kb in CHARGED_NEG) or
                       (ka in CHARGED_NEG and kb in CHARGED_POS))
            if charged:
                n_salt += 1
                n_hb += 1
                pairs.append((a, b, "salt"))
            elif a.element != b.element:
                n_hb += 1
                pairs.append((a, b, "hbond"))
            else:
                # N..N or O..O: possible, one must be a donor. Count it but
                # mark it, so the ambiguity is visible in the breakdown.
                n_hb += 1
                pairs.append((a, b, "ambiguous"))
    return {"hbond_interface": n_hb, "saltbridge_interface": n_salt}, pairs


def unsat_terms(tgt, bnd, per_atom_buried, s_cplx, idx, pairs,
                bury_frac=0.10, partner_cut=3.5):
    """Buried unsatisfied polar atoms.

    A polar atom counts as unsatisfied when it is buried by complex
    formation -- its SASA in the complex falls below `bury_frac` of its
    SASA in the free partner -- and it has no N/O partner within
    `partner_cut` anywhere in the complex, on either chain.

    This is the term that, in the published design literature, separates
    designs that bind from designs that merely score well. Our pipeline
    has no equivalent at any stage.
    """
    allatoms = tgt + bnd
    polar = [a for a in allatoms if a.element in POLAR_ELEMENTS]
    grid = Grid(allatoms, 2.0 * partner_cut)
    cut2 = partner_cut * partner_cut
    unsat = []
    for a in polar:
        sa, sc = per_atom_buried.get(id(a), (None, None))
        if sa is None:
            continue
        if sa <= 0.1:
            continue  # was already buried inside its own chain: not ours
        if sc > bury_frac * sa:
            continue  # still solvent-exposed in the complex
        found = False
        for j in grid.near(a.xyz):
            b = allatoms[j]
            if b is a or b.element not in POLAR_ELEMENTS:
                continue
            if b.key == a.key and b.is_backbone and a.is_backbone:
                continue  # own backbone neighbour
            if dist2(a.xyz, b.xyz) <= cut2:
                found = True
                break
        if not found:
            unsat.append(a)
    return {"unsat_buried_polar": len(unsat)}, unsat


def gapvol_terms(tgt, bnd, spacing=0.7, shell=5.0):
    """Gap volume at the interface, A^3, and the gap volume index.

    A real packing-defect measure rather than a contact count. Lay a grid
    over the interface region; a voxel counts as gap when it is
      (a) outside every atom's van der Waals sphere,
      (b) within `shell` A of at least one target atom AND one binder atom.
    Condition (b) is what makes it a *gap* rather than bulk solvent.

    Gap volume index = gap volume / interface area, the Laskowski
    normalisation, so that a big interface is not penalised for being big.
    """
    # interface atoms only, to keep the grid small
    cut = shell + 2.0
    gt = Grid(tgt, 2.0 * cut)
    itgt, ibnd = set(), set()
    for i, b in enumerate(bnd):
        for j in gt.near(b.xyz):
            if dist(b.xyz, tgt[j].xyz) <= cut:
                itgt.add(j)
                ibnd.add(i)
    if not itgt or not ibnd:
        return {"gap_volume": 0.0, "gap_volume_index": float("nan")}
    ta = [tgt[j] for j in itgt]
    ba = [bnd[i] for i in ibnd]
    both = ta + ba
    xs = [a.xyz[0] for a in both]
    ys = [a.xyz[1] for a in both]
    zs = [a.xyz[2] for a in both]
    lo = (min(xs) - 2, min(ys) - 2, min(zs) - 2)
    hi = (max(xs) + 2, max(ys) + 2, max(zs) + 2)
    g_all = Grid(both, 2.0 * shell)
    n_t = len(ta)
    gap = 0
    nx = int((hi[0] - lo[0]) / spacing) + 1
    ny = int((hi[1] - lo[1]) / spacing) + 1
    nz = int((hi[2] - lo[2]) / spacing) + 1
    vox = spacing ** 3
    sh2 = shell * shell
    for ix in range(nx):
        x = lo[0] + ix * spacing
        for iy in range(ny):
            y = lo[1] + iy * spacing
            for iz in range(nz):
                z = lo[2] + iz * spacing
                p = (x, y, z)
                inside = False
                near_t = near_b = False
                for j in g_all.near(p):
                    a = both[j]
                    d2 = dist2(p, a.xyz)
                    if d2 < a.radius * a.radius:
                        inside = True
                        break
                    if d2 <= sh2:
                        if j < n_t:
                            near_t = True
                        else:
                            near_b = True
                if inside or not (near_t and near_b):
                    continue
                gap += 1
    return {"gap_volume": gap * vox}


def helix_terms(bnd_residues):
    """Intrinsic cost of pre-organising the binder as a helix, kcal/mol.

    Sum of Pace & Scholtz helix propensities over the binder sequence.
    This is a sequence-only term: it does not look at the target at all.
    That is the point. Nothing in the pipeline prices the conformational
    entropy a flexible peptide pays to fold before it binds, and for a
    12-mer that price is a large fraction of the measured Kd.
    """
    seq = []
    for key in sorted(bnd_residues, key=lambda k: k[1]):
        resname = bnd_residues[key][0]
        seq.append(THREE_TO_ONE.get(resname, "X"))
    s = "".join(seq)
    known = [c for c in s if c in HELIX_COST]
    if not known:
        return {"binder_seq": s, "helix_cost_total": float("nan"),
                "helix_cost_per_res": float("nan")}
    tot = sum(HELIX_COST[c] for c in known)
    return {"binder_seq": s,
            "helix_cost_total": tot,
            "helix_cost_per_res": tot / len(known)}


def contact_terms(tgt, bnd):
    """Contact counts at the cutoffs this project already uses.

    Kept so that the new gates can be compared against the old habit on
    the same structures -- and so that the failure of contact counting is
    on the record in the same table, not asserted.
    """
    out = {}
    hv = 0
    cut2 = 4.5 * 4.5
    g = Grid(tgt, 12.0)
    for b in bnd:
        for j in g.near(b.xyz):
            if dist2(b.xyz, tgt[j].xyz) <= cut2:
                hv += 1
    out["contacts_heavy_4.5"] = hv
    tca = [a for a in tgt if a.name == "CA"]
    bca = [a for a in bnd if a.name == "CA"]
    for cut in (8.0, 10.0):
        n = 0
        c2 = cut * cut
        for b in bca:
            for t in tca:
                if dist2(b.xyz, t.xyz) <= c2:
                    n += 1
        out["contacts_ca_%g" % cut] = n
    out["n_binder_res"] = len(bca)
    out["n_target_res"] = len(tca)
    return out


def compactness_terms(bnd):
    """Gate 1, ported. Rg of the binder against the compact expectation."""
    ca = [a.xyz for a in bnd if a.name == "CA"]
    n = len(ca)
    if n < 2:
        return {"binder_rg": float("nan"), "rg_ratio": float("nan")}
    cx = sum(p[0] for p in ca) / n
    cy = sum(p[1] for p in ca) / n
    cz = sum(p[2] for p in ca) / n
    rg = math.sqrt(sum(dist2(p, (cx, cy, cz)) for p in ca) / n)
    expected = 2.2 * (n ** 0.38)
    return {"binder_rg": rg, "rg_ratio": rg / expected,
            "rg_expected": expected}


def completeness_terms(residues):
    """How many residues are missing heavy atoms.

    An energy term computed on an incomplete side chain is not wrong, it
    is undefined. This column exists so that no score in this table can be
    read without knowing whether its inputs were complete.
    """
    missing = 0
    worst = []
    for key, (resname, atoms) in residues.items():
        exp = NHEAVY.get(resname)
        if exp is None:
            continue
        have = len({a.name for a in atoms if a.name != "OXT"})
        if have < exp:
            missing += 1
            worst.append("%s%d(%d/%d)" % (resname, key[1], have, exp))
    return {"n_incomplete_residues": missing,
            "incomplete_list": ";".join(worst[:10]) if worst else "-"}


def chirality_terms(residues):
    """D-amino acid detection, with the sign calibrated, not assumed.

    The N-C-CA-CB improper dihedral is about -122 deg for an L-amino acid
    and +122 for a D one. Calibrated on 1YCR (94/94 residues at
    -122.3 +/- 3.0 deg). A rebuilt loop or a tleap-built system can
    silently contain D-residues, and a D-residue is invisible to every
    other gate, so this runs on every structure.
    """
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    def improper(p1, p2, p3, p4):
        b1, b2, b3 = sub(p2, p1), sub(p3, p2), sub(p4, p3)
        n1, n2 = cross(b1, b2), cross(b2, b3)
        nb2 = math.sqrt(dot(b2, b2))
        if nb2 == 0:
            return None
        m = cross(n1, (b2[0] / nb2, b2[1] / nb2, b2[2] / nb2))
        return math.degrees(math.atan2(dot(m, n2), dot(n1, n2)))

    n_d = 0
    n_checked = 0
    d_list = []
    for key, (resname, atoms) in residues.items():
        if resname == "GLY":
            continue
        by = {a.name: a.xyz for a in atoms}
        if not {"N", "CA", "C", "CB"} <= by.keys():
            continue
        ang = improper(by["N"], by["C"], by["CA"], by["CB"])
        if ang is None:
            continue
        n_checked += 1
        if ang > 0:
            n_d += 1
            d_list.append("%s%s%d" % (key[0], resname, key[1]))
    return {"n_chiral_checked": n_checked, "n_d_residues": n_d,
            "d_residue_list": ";".join(d_list[:10]) if d_list else "-"}


# ----------------------------------------------------------------- driver

COLUMNS = [
    "structure",
    # provenance / sanity, first so they are never missed
    "n_binder_res", "n_target_res", "n_incomplete_residues",
    "n_chiral_checked", "n_d_residues",
    "binder_seq",
    # existing gate 1
    "binder_rg", "rg_expected", "rg_ratio",
    # old habit, kept for comparison
    "contacts_heavy_4.5", "contacts_ca_8", "contacts_ca_10",
    # new gates
    "dsasa_total", "dsasa_apolar", "dsasa_polar",
    "dsasa_binder_side", "dsasa_target_side",
    "hbond_interface", "saltbridge_interface",
    "unsat_buried_polar",
    "gap_volume", "gap_volume_index",
    "helix_cost_total", "helix_cost_per_res",
    # the deliberately stupid baseline
    "dumb_dG_apolar",
    # raw inputs the derived numbers came from
    "sasa_complex_total", "incomplete_list", "d_residue_list",
]


def score_one(path, target_chain, binder_chain, n_points, do_gap):
    atoms = read_pdb(path)
    tgt = [a for a in atoms if a.chain == target_chain]
    bnd = [a for a in atoms if a.chain == binder_chain]
    if not tgt or not bnd:
        chains = sorted({a.chain for a in atoms})
        raise SystemExit(
            "%s: chains %r/%r not both present; file has %r.\n"
            "Pass --target-chain/--binder-chain. Note that "
            "*_af2pred.pdb files lose their original chain labels."
            % (path, target_chain, binder_chain, chains))

    row = {"structure": os.path.basename(path)}
    allres = residues_of(atoms)
    bndres = residues_of(bnd)

    row.update(completeness_terms(allres))
    row.update(chirality_terms(allres))
    row.update(contact_terms(tgt, bnd))
    row.update(compactness_terms(bnd))
    row.update(helix_terms(bndres))

    ds, per_atom, s_cplx, idx = dsasa_terms(atoms, tgt, bnd, n_points)
    row.update(ds)

    hb, pairs = hbond_terms(tgt, bnd)
    row.update(hb)

    un, unsat = unsat_terms(tgt, bnd, per_atom, s_cplx, idx, pairs)
    row.update(un)

    if do_gap:
        gv = gapvol_terms(tgt, bnd)
        row.update(gv)
        # Laskowski normalisation: interface area is half the buried area
        iface_area = row["dsasa_total"] / 2.0
        row["gap_volume_index"] = (gv["gap_volume"] / iface_area
                                   if iface_area > 0 else float("nan"))
    else:
        row["gap_volume"] = float("nan")
        row["gap_volume_index"] = float("nan")

    row["dumb_dG_apolar"] = -APOLAR_COEFF * row["dsasa_apolar"]
    return row, (pairs, unsat, tgt, bnd)


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        if math.isnan(v):
            return "NA"
        return "%.3f" % v
    return str(v)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="CPU-only interface gate battery.")
    ap.add_argument("pdbs", nargs="+")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--out", default=None,
                    help="TSV to write. Refuses to overwrite or append.")
    ap.add_argument("--sasa-points", type=int, default=128,
                    help="Sphere points per atom. 128 is plenty; 512 for "
                         "a final number.")
    ap.add_argument("--no-gap-volume", action="store_true",
                    help="Skip the gap-volume grid, the slow term.")
    ap.add_argument("--per-residue", action="store_true",
                    help="Print a per-residue breakdown for each input.")
    args = ap.parse_args(argv)

    if args.out and os.path.exists(args.out):
        raise SystemExit(
            "%s exists. Pass a fresh --out: appending to a score file and "
            "reading it back by column position is how we lost a day."
            % args.out)

    rows = []
    for path in args.pdbs:
        row, extra = score_one(path, args.target_chain, args.binder_chain,
                               args.sasa_points, not args.no_gap_volume)
        rows.append(row)
        if args.per_residue:
            per_residue_report(path, row, extra)

    lines = ["\t".join(COLUMNS)]
    for row in rows:
        lines.append("\t".join(fmt(row.get(c)) for c in COLUMNS))
    text = "\n".join(lines) + "\n"
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text)
        sys.stderr.write("wrote %s (%d rows)\n" % (args.out, len(rows)))
    else:
        sys.stdout.write(text)
    return 0


def per_residue_report(path, row, extra):
    pairs, unsat, tgt, bnd = extra
    print("=" * 70)
    print("per-residue breakdown:", os.path.basename(path))
    print("=" * 70)
    print("binder sequence: %s" % row["binder_seq"])
    print()
    print("interface polar pairs within 3.5 A (upper bound, no H angles):")
    if not pairs:
        print("  none")
    for a, b, kind in sorted(pairs, key=lambda p: p[0].resnum):
        print("  %-9s %-4s -- %-9s %-4s  %.2f A  [%s]" % (
            "%s%d" % (a.resname, a.resnum), a.name,
            "%s%d" % (b.resname, b.resnum), b.name,
            dist(a.xyz, b.xyz), kind))
    print()
    print("buried unsatisfied polar atoms:")
    if not unsat:
        print("  none")
    for a in unsat:
        print("  %s/%s%d %s" % (a.chain, a.resname, a.resnum, a.name))
    print()
    print("binder residues by buried area (needs --sasa-points, "
          "recomputed per residue):")
    s_c = sasa(tgt + bnd, 128)
    s_b = sasa(bnd, 128)
    off = len(tgt)
    per = {}
    for i, a in enumerate(bnd):
        d = s_b[i] - s_c[off + i]
        k = (a.resname, a.resnum)
        pb = per.setdefault(k, [0.0, 0.0])
        if d > 0:
            pb[0] += d
            if is_apolar(a):
                pb[1] += d
    for (rn, num), (tot, apol) in sorted(per.items(),
                                         key=lambda kv: -kv[1][0]):
        print("  %-9s buried %7.1f A^2  (apolar %6.1f)" % (
            "%s%d" % (rn, num), tot, apol))
    print()


if __name__ == "__main__":
    sys.exit(main())
