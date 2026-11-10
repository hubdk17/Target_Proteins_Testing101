#!/usr/bin/env python3
"""
Hotspot coverage must not be a by-product of the CA-CA grip test.

The fixture is built so that NO CA-CA pair falls inside the 8 A grip shell
(grip == 0), while one hotspot's side chain reaches to 6 A of a binder CA.
Before the fix, `hot_hit` was populated inside `if dd <= g2`, so it was a
subset of the grip pairs and this hotspot was reported as missed. If this test
ever runs with grip > 0 it would pass for the wrong reason and prove nothing,
hence the explicit assertion on grip.

    python3 tests/test_grip_hotspot.py     # or: python3 -m pytest tests/
"""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, os.pardir, "42_grip_rank.py")
FIXTURE = os.path.join(HERE, "fixture_hotspot_outside_ca_shell.pdb")

_spec = importlib.util.spec_from_file_location("grip", SCRIPT)
grip = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(grip)


def test_hotspot_coverage_independent_of_grip():
    m = grip.measure(FIXTURE, [157, 120])
    assert m is not None, "fixture needs both chain A and chain B"

    assert m["grip"] == 0, (
        "fixture is meant to have no CA-CA contact inside the grip shell; "
        "got grip=%d, so this test would pass for the wrong reason"
        % m["grip"])

    assert m["hot_hit"] == [157], (
        "hotspot 157 reaches 6 A via its side chain and must be reported as "
        "engaged even with grip==0; got %r" % (m["hot_hit"],))

    assert abs(m["hot_near"][157] - 6.0) < 1e-6, m["hot_near"][157]
    assert abs(m["hot_near"][120] - 40.0) < 1e-6, m["hot_near"][120]
    assert 120 not in m["hot_hit"], "the distal control must not be engaged"


if __name__ == "__main__":
    test_hotspot_coverage_independent_of_grip()
    print("ok: hotspot coverage is independent of the CA-CA grip count")
