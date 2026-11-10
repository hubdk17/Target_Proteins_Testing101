#!/usr/bin/env python3
"""Tests for the studio layer.

House style on this project: design checks that can actually fail. So these
tests are mostly about the refusals - the cases where the right behaviour is an
exception or a REFUSED verdict rather than a number. Every one of them encodes a
mistake this pipeline has already made once.

    python3 -m pytest tests/ -q
    python3 tests/test_studio.py        # no pytest needed
"""
from __future__ import annotations

import math
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from studio import checks, provenance, residues, scores, stats  # noqa: E402
from studio.config import load  # noqa: E402
from studio.gatekeeper import ForbiddenRanking, Gatekeeper, Verdict  # noqa: E402
from studio.store import Store  # noqa: E402

CFG = load()
GK = Gatekeeper(CFG)


# --------------------------------------------------------------------- helpers

def atom(serial, name, resname, chain, resnum, x, y, z):
    return (
        f"ATOM  {serial:5d} {name:<4s}{resname:>4s}{chain:>2s}{resnum:4d}"
        f"    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]}"
    )


def write_pdb(lines, suffix=".pdb"):
    fh = tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False)
    fh.write("\n".join(lines) + "\nEND\n")
    fh.close()
    return Path(fh.name)


def one_residue(resnum, chain="A", resname="ALA", flip_cb=False, origin=(0.0, 0.0, 0.0)):
    """An L residue, or its mirror image if flip_cb.

    With N, CA and C laid out in the z=0 plane as below, the signed volume
    (N-CA) x (C-CA) . (CB-CA) points along -z, so CB at negative z is the
    positive (L) case. Checked by hand rather than guessed: getting this
    backwards is how a mirror-image fixture silently tests nothing.
    """
    ox, oy, oz = origin
    cb_z = 1.2 if flip_cb else -1.2
    return [
        atom(1, "N", resname, chain, resnum, ox - 1.2, oy + 0.6, oz),
        atom(2, "CA", resname, chain, resnum, ox, oy, oz),
        atom(3, "C", resname, chain, resnum, ox + 1.2, oy + 0.7, oz),
        atom(4, "O", resname, chain, resnum, ox + 1.3, oy + 1.9, oz),
        atom(5, "CB", resname, chain, resnum, ox, oy - 0.8, oz + cb_z),
    ]


def ca_chain(numbers, spacing=3.8, chain="A", resname="ALA"):
    """CA-only backbone with a given residue numbering and CA-CA spacing."""
    out = []
    for i, num in enumerate(numbers):
        out.append(atom(i + 1, "CA", resname, chain, num, i * spacing, 0.0, 0.0))
    return out


# ----------------------------------------------------------------------- config

def test_config_has_one_hash_and_all_sections():
    cfg = load()
    assert cfg.hash.startswith("sha256:")
    assert cfg.version
    assert cfg.gate1["rho_max"] == 1.30
    assert cfg.gate3["min_replicates"] >= 3
    assert cfg.stamp()["config_hash"] == cfg.hash


def test_config_hash_changes_with_content():
    """Two different threshold files must not share a hash, or a verdict could
    claim to have been judged by numbers it was not judged by."""
    body = CFG.path.read_text().replace("rho_max = 1.30", "rho_max = 1.40")
    fh = tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False)
    fh.write(body)
    fh.close()
    tweaked = Path(fh.name)
    try:
        other = load(tweaked)
        assert other.hash != CFG.hash
        assert other.gate1["rho_max"] == 1.40
    finally:
        tweaked.unlink()


# ------------------------------------------------------------------- scores.py

SC_REORDERED = """\
SCORE:     pae_interaction  binder_aligned_rmsd  plddt_binder description
SCORE:               8.100                1.400        78.00 d1_s1_af2pred
SCORE:              23.000                5.900        94.00 d1_s2_af2pred
SCORE:     binder_aligned_rmsd  pae_interaction  plddt_binder description
SCORE:                   0.800            4.200        76.00 d2_s1_af2pred
"""


def test_score_columns_are_read_by_name_not_position():
    """The bug this module exists for: `41_gate2_report.sh` uses $2 as the RMSD
    regardless of the header, so a reordered score file hands it a real number
    from a different column. Here the two blocks below have the RMSD in
    different positions and both must still read correctly."""
    p = write_pdb([SC_REORDERED.rstrip()], suffix=".sc")
    try:
        rows = scores.read_sc(p)
        assert len(rows) == 3
        rmsd = scores.column(rows, scores.GATE2_RMSD, float)
        assert rmsd == [1.4, 5.9, 0.8]
        pae = scores.column(rows, scores.GATE3_PAE, float)
        assert pae == [8.1, 23.0, 4.2]
    finally:
        p.unlink()


def test_missing_column_raises_rather_than_returning_zero():
    p = write_pdb([SC_REORDERED.rstrip()], suffix=".sc")
    try:
        rows = scores.read_sc(p)
        try:
            scores.column(rows, "iptm", float)
        except scores.ScoreFileError as e:
            assert "iptm" in str(e)
        else:
            raise AssertionError("a missing column must raise, not default to 0.0")
    finally:
        p.unlink()


def test_backbone_of_strips_sequence_and_predictor_suffixes():
    assert scores.backbone_of("pd1_b1_3_s7_af2pred") == "pd1_b1_3"
    assert scores.backbone_of("pd1_b1_3_s12") == "pd1_b1_3"
    assert scores.backbone_of("pd1_b1_3") == "pd1_b1_3"


# -------------------------------------------------------------------- stats.py

def test_wilson_is_sane_at_the_edges():
    """0/50 must not come out as a zero-width interval; brca1 was 0/50 and the
    honest reading is 0 - 7.1%, not 'zero'."""
    r = stats.wilson(0, 50)
    assert r.low == 0.0
    assert 0.06 < r.high < 0.08
    r = stats.wilson(20, 50)
    assert 0.27 < r.low < 0.29 and 0.53 < r.high < 0.55


def test_overlapping_cells_are_not_ranked():
    """The six-target ladder at n=50 per cell.

    pcsk9 40% and glp1r 38% are one tier, not two. On Wilson intervals pd1 20%
    is not separable from pcsk9 either, so the honest grouping is two tiers
    rather than the three that pairwise z-tests suggest. The low group does
    separate from the top, which is the one statement the data supports.
    """
    tiers = GK.rank_cells({
        "pcsk9": (20, 50), "glp1r": (19, 50), "pd1": (10, 50),
        "bcl2": (3, 50), "casp3": (3, 50), "brca1": (0, 50),
    })
    flat = {name: i for i, tier in enumerate(tiers) for name in tier}
    assert flat["pcsk9"] == flat["glp1r"], "pcsk9 must not be ranked above glp1r"
    assert flat["bcl2"] == flat["casp3"] == flat["brca1"]
    assert flat["pcsk9"] < flat["bcl2"], "the top and bottom groups do separate"
    assert len(tiers) == 2, f"expected two defensible tiers, got {tiers}"


def test_small_cells_are_flagged_underpowered_not_null():
    r = GK.pass_rate(0, 10)
    assert r.underpowered and "UNDERPOWERED" in str(r)


def test_unequal_arms_refuse_to_be_compared():
    try:
        stats.require_equal_n({"wt": [1.0] * 150, "mut": [1.0] * 150, "ctl": [1.0] * 97})
    except stats.UnequalArms as e:
        assert "ctl=97" in str(e)
    else:
        raise AssertionError("unequal n must refuse, not compare")
    assert stats.require_equal_n({"wt": [1.0] * 150, "mut": [1.0] * 150}) == 150


# ---------------------------------------------------------------- gate verdicts

def test_gate1_separates_compact_borderline_and_rod():
    assert GK.gate1("d1", 1.05).verdict is Verdict.PASS
    assert GK.gate1("d2", 1.31).verdict is Verdict.BORDERLINE
    rod = GK.gate1("d3", 2.2)
    assert rod.verdict is Verdict.FAIL and "rod" in rod.reason


def test_gate2_refuses_a_backbone_sampled_once():
    """Within-backbone spread is 11-20 units, so one sequence is not a verdict.
    Backbone 87 scored 23.0 on its single sample and is the pool's best."""
    j = GK.gate2_backbone("bb87", {"s1": 1.2})
    assert j.verdict is Verdict.REFUSED and not j.verdict.carries_forward
    many = {f"s{i}": 3.0 for i in range(1, 9)}
    many["s4"] = 1.1
    ok = GK.gate2_backbone("bb87", many)
    assert ok.verdict is Verdict.PASS and ok.detail["best_sequence"] == "s4"


def test_gate3_refuses_a_single_replicate():
    j = GK.gate3("d1", [8.0])
    assert j.verdict is Verdict.REFUSED
    assert "1 in 10" in j.reason


def test_gate3_borderline_covers_the_noise_band():
    """A design at 9.4 with a wide spread and one at 4.1 with a tight spread must
    not look the same."""
    near = GK.gate3("near", [9.4, 12.2, 6.6])
    far = GK.gate3("far", [4.1, 4.6, 3.7])
    assert near.verdict is Verdict.BORDERLINE
    assert far.verdict is Verdict.PASS
    assert near.verdict.carries_forward and far.verdict.carries_forward
    assert GK.gate3("bad", [18.0, 19.0, 20.0]).verdict is Verdict.FAIL


def test_gate3_reports_saturation_instead_of_a_number():
    j = GK.gate3("d1", [31.8, 32.0, 31.9])
    assert j.verdict is Verdict.SATURATED and not j.verdict.carries_forward


def test_gate3_marks_whether_the_tight_noise_floor_was_reached():
    assert GK.gate3("d1", [4.0] * 3).detail["tight_floor"] is False
    assert GK.gate3("d2", [4.0] * 9).detail["tight_floor"] is True


def test_every_verdict_carries_the_threshold_hash():
    for j in (GK.gate1("d", 1.0), GK.gate3("d", [4.0] * 3),
              GK.chirality("d", 0, 100)):
        assert j.config_hash == CFG.hash
        assert j.thresholds_version == CFG.version


# -------------------------------------------------------------------- ranking

def test_confidence_sorting_is_refused_with_no_override():
    rows = [{"design": "a", "plddt_binder": 94.0}, {"design": "b", "plddt_binder": 78.0}]
    for key in ("plddt_binder", "PLDDT", "iptm", "mpnn_score", "-confidence"):
        try:
            GK.rank(rows, key)
        except ForbiddenRanking as e:
            assert "2.6%" in str(e)
        else:
            raise AssertionError(f"ranking by {key} must be refused")


def test_ranking_refuses_to_fall_back_to_filename_order():
    """`43_build_complex.py --top N` slices alphabetically while its docstring
    promises grip order. Here a missing key is an error."""
    rows = [{"design": "a", "grip": 31}, {"design": "b"}]
    try:
        GK.rank(rows, "grip", descending=True)
    except KeyError as e:
        assert "grip" in str(e)
    else:
        raise AssertionError("a missing rank key must raise")
    ranked = GK.rank([{"design": "a", "grip": 31}, {"design": "b", "grip": 44}],
                     "grip", descending=True)
    assert [r["design"] for r in ranked] == ["b", "a"]


# ------------------------------------------------------------------- residues

def test_relabelling_trp_to_tyr_is_a_chimera():
    """The worked example from PARSER-ATOM-NAME-LIMITATION.md."""
    kept, dropped, never = residues.survives_relabel("TRP", "TYR")
    # The doc names the side-chain atoms; the backbone is shared by every label.
    assert set(kept) == {"CG", "CD1", "CD2", "CE2"} | {"N", "CA", "C", "O", "CB"}
    assert dropped == sorted(["NE1", "CE3", "CZ2", "CZ3", "CH2"])
    assert never == sorted(["CE1", "CZ", "OH"])


def test_ala_can_only_ever_load_five_heavy_atoms():
    assert len(residues.CANONICAL["ALA"]) == 5
    assert len(residues.CANONICAL["TRP"]) == residues.MAX_HEAVY_SLOTS


# --------------------------------------------------------------------- checks

def test_atom_coverage_finds_dropped_atoms_and_missing_ca():
    lines = one_residue(1)
    lines += [atom(6, "ZZ1", "ALA", "A", 1, 9.0, 9.0, 9.0)]        # dropped at load
    lines += [atom(7, "N", "ALA", "A", 2, 0.0, 0.0, 0.0),           # no CA: vanishes
              atom(8, "C", "ALA", "A", 2, 1.0, 0.0, 0.0)]
    lines += [atom(9, "CA", "XYZ", "A", 3, 5.0, 0.0, 0.0)]          # non-standard label
    p = write_pdb(lines)
    try:
        cov = checks.atom_coverage(p)
        assert any("ZZ1" in d for d in cov.dropped)
        assert any("ALA2" in m for m in cov.no_ca)
        assert any("XYZ3" in n for n in cov.nonstandard)
        assert not cov.clean
    finally:
        p.unlink()


def test_chain_walk_tells_a_numbering_gap_from_a_real_break():
    """IL-6's shape: residue numbers jump 51->61 at a real break, while other
    number jumps are only convention. Measured, never assumed."""
    # 19..21 contiguous, then a number jump WITH normal spacing (numbering gap),
    # then a number jump WITH a long distance (real break).
    lines = [atom(1, "CA", "ALA", "A", 19, 0.0, 0.0, 0.0),
             atom(2, "CA", "ALA", "A", 20, 3.8, 0.0, 0.0),
             atom(3, "CA", "ALA", "A", 25, 7.6, 0.0, 0.0),     # gap, not a break
             atom(4, "CA", "ALA", "A", 61, 30.0, 0.0, 0.0)]    # real break
    p = write_pdb(lines)
    try:
        walk = checks.chain_walk(p, "A", break_cutoff=CFG.structure["break_cutoff_angstrom"])
        assert [(a, b) for a, b, _ in walk.numbering_gaps] == [(20, 25)]
        assert [(a, b) for a, b, _ in walk.real_breaks] == [(25, 61)]
        assert walk.contig == "A19-25/A61-61"
    finally:
        p.unlink()


def test_chain_walk_keeps_a_contig_whole_across_a_pure_numbering_gap():
    p = write_pdb(ca_chain([1, 2, 3, 10, 11, 12]))
    try:
        walk = checks.chain_walk(p, "A")
        assert walk.real_breaks == []
        assert walk.contig == "A1-12"
    finally:
        p.unlink()


def test_chirality_wrapper_reads_the_scripts_own_report():
    good = write_pdb(one_residue(1) + one_residue(2, origin=(5.0, 0.0, 0.0)))
    bad = write_pdb(one_residue(1) + one_residue(2, flip_cb=True, origin=(5.0, 0.0, 0.0)))
    try:
        g = checks.chirality(good)
        assert g.n_checked == 2 and g.n_d == 0
        assert GK.chirality("ok", g.n_d, g.n_checked).verdict is Verdict.PASS

        b = checks.chirality(bad)
        assert b.n_checked == 2 and b.n_d == 1, f"expected one D residue, got {b.n_d}"
        j = GK.chirality("bad", b.n_d, b.n_checked, b.d_residues)
        assert j.verdict is Verdict.FAIL and not j.verdict.carries_forward
    finally:
        good.unlink()
        bad.unlink()


def test_a_mirrored_file_is_refused_not_reported_as_all_d():
    lines = []
    for i in range(1, 5):
        lines += one_residue(i, flip_cb=True, origin=(5.0 * i, 0.0, 0.0))
    p = write_pdb(lines)
    try:
        c = checks.chirality(p)
        assert c.mirrored_warning
        j = GK.chirality("mirrored", c.n_d, c.n_checked)
        assert j.verdict is Verdict.REFUSED
    finally:
        p.unlink()


def test_a_chirality_parse_failure_raises_instead_of_passing():
    """A check that silently reports 'no problems' because its parse failed is
    the exact failure shape this project catalogues."""
    empty = write_pdb(["REMARK nothing here"])
    try:
        try:
            checks.chirality(empty)
        except checks.CheckError:
            pass
        else:
            raise AssertionError("an unparseable report must raise")
    finally:
        empty.unlink()


def test_stage_success_is_artifact_count_not_exit_code():
    ok = checks.verify_output(expected=100, produced=100, elapsed_s=4_600,
                              exit_code=0, min_seconds_per_item=1.0)
    assert ok.ok and ok.kind == "finished"

    silent = checks.verify_output(expected=100, produced=12, elapsed_s=4_600, exit_code=0)
    assert not silent.ok and silent.kind == "silent_failure"
    assert "despite exit 0" in silent.message

    never_started = checks.verify_output(expected=100, produced=0, elapsed_s=45, exit_code=0)
    assert never_started.kind == "hard_failure"

    too_fast = checks.verify_output(expected=100, produced=100, elapsed_s=3,
                                    exit_code=0, min_seconds_per_item=1.0)
    assert not too_fast.ok and too_fast.kind == "silent_failure"


# ---------------------------------------------------------------------- store

def test_two_runs_on_one_target_cannot_mix_their_pools():
    """Code-review findings 2, 3 and 4 in one test: `TAG` never reached the
    derived directories, so two tags on one target overwrote each other's gate
    lists; and `gate2/<name>.pass` was never truncated, so a run where nothing
    passed printed the previous run's survivors as its own."""
    with tempfile.TemporaryDirectory() as tmp:
        st = Store(tmp)
        a = st.create_run("il6", {"tag": "b1"}, CFG.stamp(), 50.0)
        b = st.create_run("il6", {"tag": "b2"}, CFG.stamp(), 50.0)
        assert a.run_id != b.run_id and a.root != b.root

        st.record(a.run_id, [GK.gate1("il6_1", 1.0), GK.gate1("il6_2", 1.1)])
        st.record(b.run_id, [GK.gate1("il6_1", 2.2)])      # nothing passes in run b

        assert st.survivors(a.run_id, "gate1") == ["il6_1", "il6_2"]
        assert st.survivors(b.run_id, "gate1") == [], "run b must not inherit run a's survivors"


def test_refused_designs_do_not_carry_forward_but_borderline_does():
    with tempfile.TemporaryDirectory() as tmp:
        st = Store(tmp)
        r = st.create_run("il6", {}, CFG.stamp(), 50.0)
        st.record(r.run_id, [
            GK.gate3("pass", [4.0, 4.2, 3.9]),
            GK.gate3("borderline", [9.4, 12.2, 6.6]),
            GK.gate3("refused", [8.0]),
            GK.gate3("fail", [18.0, 19.0, 20.0]),
        ])
        assert st.survivors(r.run_id, "gate3") == ["borderline", "pass"]
        funnel = st.funnel(r.run_id)["gate3"]
        assert funnel == {"pass": 1, "borderline": 1, "refused": 1, "fail": 1}


def test_spend_accumulates_and_unknown_event_kinds_raise():
    with tempfile.TemporaryDirectory() as tmp:
        st = Store(tmp)
        r = st.create_run("il6", {}, CFG.stamp(), 50.0)
        st.charge(r.run_id, "stage1", 4.6, 8.23, "rfdiffusion")
        total = st.charge(r.run_id, "gate2", 9.3, 16.65, "self-consistency")
        assert math.isclose(total, 24.88)
        st.event(r.run_id, "silent_failure", "12 of 100 structures", stage="stage1")
        assert st.events(r.run_id)[0]["kind"] == "silent_failure"
        try:
            st.event(r.run_id, "looks_fine", "typo kind")
        except ValueError:
            pass
        else:
            raise AssertionError("an unknown event kind must raise")


# ----------------------------------------------------------------- provenance

def test_incomplete_provenance_refuses_to_be_written():
    with tempfile.TemporaryDirectory() as tmp:
        p = provenance.Provenance(run_id="r1", target="il6", **CFG.stamp())
        assert "image_digest" in p.missing()
        try:
            p.write(Path(tmp) / "prov.json")
        except ValueError as e:
            assert "image_digest" in str(e)
        else:
            raise AssertionError("an incomplete provenance record must refuse")

        weight = Path(tmp) / "params.bin"
        weight.write_bytes(b"not really 5 GB")
        p.image_digest = "sha256:deadbeef"
        p.add_weight("af2_monomer", weight)
        p.seeds = {"stage1": [1, 2, 3]}
        p.artifact_counts = {"stage1": 100}
        p.refuse("stage0", "two partner chains in contact and none chosen")
        out = p.write(Path(tmp) / "prov.json")
        assert out.exists()
        assert p.weights["af2_monomer"].startswith("sha256:")
        assert p.refusals and p.refusals[0]["where"] == "stage0"


# --------------------------------------------------------------------- runner

def _run_all():
    tests = [(k, v) for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    failed = []
    for name, fn in tests:
        try:
            fn()
            print(f"  ok    {name}")
        except Exception as e:  # noqa: BLE001
            failed.append((name, e))
            print(f"  FAIL  {name}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - len(failed)} passed, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(_run_all())
