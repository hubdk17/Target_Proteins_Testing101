#!/usr/bin/env python3
"""The statistics the studio is allowed to show, and the comparisons it refuses.

Three rules this project paid for:

- A pass rate is a Wilson interval, not a percentage. At 50 designs per cell the
  six-target ladder resolves about three tiers and no finer, so "pcsk9 40% vs
  glp1r 38%" is not a difference.
- A cell with too few items is reported as UNDERPOWERED, never as a null.
- Arms with unequal n are not compared at all. Unequal n between arms is how
  this stack misleads.

Stdlib only; the normal quantile is a closed-form approximation so nothing here
needs scipy.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Rate:
    passed: int
    total: int
    low: float
    high: float
    underpowered: bool

    @property
    def point(self) -> float:
        return self.passed / self.total if self.total else float("nan")

    def __str__(self) -> str:
        if not self.total:
            return "no data"
        s = f"{self.point:.1%} ({self.passed}/{self.total}, 95% CI {self.low:.1%}-{self.high:.1%})"
        return s + " UNDERPOWERED" if self.underpowered else s


@dataclass(frozen=True)
class Spread:
    n: int
    mean: float
    sd: float
    sem: float
    lo: float
    hi: float

    def __str__(self) -> str:
        return f"{self.mean:.2f} +/- {self.sd:.2f} (n={self.n}, range {self.lo:.2f}-{self.hi:.2f})"


def z_for(confidence: float = 0.95) -> float:
    """Two-sided normal quantile. Abramowitz & Stegun 26.2.23, |error| < 4.5e-4."""
    p = 1.0 - (1.0 - confidence) / 2.0
    if not 0.5 < p < 1.0:
        raise ValueError(f"confidence {confidence} out of range")
    t = math.sqrt(-2.0 * math.log(1.0 - p))
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    return t - ((c2 * t + c1) * t + c0) / (((d3 * t + d2) * t + d1) * t + 1.0)


def wilson(passed: int, total: int, confidence: float = 0.95,
           underpowered_below: int = 30) -> Rate:
    """Wilson score interval. Correct at 0/n and n/n, where the normal
    approximation gives a zero-width interval and invites a false claim."""
    if passed < 0 or total < 0 or passed > total:
        raise ValueError(f"{passed} of {total} is not a rate")
    if total == 0:
        return Rate(0, 0, float("nan"), float("nan"), True)

    z = z_for(confidence)
    p = passed / total
    denom = 1.0 + z * z / total
    centre = (p + z * z / (2.0 * total)) / denom
    half = (z / denom) * math.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total))
    return Rate(
        passed=passed,
        total=total,
        low=max(0.0, centre - half),
        high=min(1.0, centre + half),
        underpowered=total < underpowered_below,
    )


def comparable(a: Rate, b: Rate) -> bool:
    """False when the intervals overlap, i.e. when the two must not be ranked."""
    if not a.total or not b.total:
        return False
    return a.low > b.high or b.low > a.high


def spread(values: list[float]) -> Spread:
    n = len(values)
    if n == 0:
        raise ValueError("no values")
    mean = sum(values) / n
    if n == 1:
        return Spread(1, mean, float("nan"), float("nan"), mean, mean)
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    sd = math.sqrt(var)
    return Spread(n, mean, sd, sd / math.sqrt(n), min(values), max(values))


class UnequalArms(Exception):
    """Raised instead of returning a comparison across arms of different n."""


def require_equal_n(arms: dict[str, list]) -> int:
    """Refuse an arm comparison unless every arm has the same n.

    `52_mutation_report.py` prints each arm's n and then compares regardless.
    This raises, because an unequal-n comparison looks exactly like a real
    finding. Read only complete groups.
    """
    if not arms:
        raise UnequalArms("no arms given")
    sizes = {name: len(values) for name, values in arms.items()}
    if len(set(sizes.values())) != 1:
        detail = ", ".join(f"{k}={v}" for k, v in sorted(sizes.items()))
        raise UnequalArms(
            f"arms have unequal n ({detail}); refusing to compare. "
            "Drop to the smallest complete group or finish the short arms."
        )
    n = next(iter(sizes.values()))
    if n == 0:
        raise UnequalArms("every arm is empty")
    return n
