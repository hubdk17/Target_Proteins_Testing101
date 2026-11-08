#!/usr/bin/env python3
"""The only thing in the studio allowed to write a verdict.

Everything else measures. This module decides, from one versioned threshold
file, and stamps that file's hash on every verdict it writes. If a number ever
needs explaining, the row says which thresholds judged it.

Four deliberate design choices, each one a trap this project already fell into:

- A verdict is never a bare boolean. `BORDERLINE` exists because gate 3 is not
  deterministic where it decides: roughly 1 in 10 repeats flips pass/fail, and a
  single design needs ~3 units to be resolvable. A design at 9.4 +/- 2.8 and one
  at 4.1 +/- 0.6 must not look the same on screen.
- Gate 3 REFUSES below `min_replicates`. One draw is not a verdict.
- Ranking by any confidence-like key raises. `plddt_binder` explains ~2.6% of
  fold-back variance and inverts on 17 of 47 backbones, so a "best designs" list
  sorted by confidence is near random. There is no flag to turn this off.
- A saturated gate-3 score is reported as saturated, not as a number, because
  above ~31.75 the score no longer separates anything.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any, Sequence

from . import stats
from .config import Config, load


class Verdict(str, Enum):
    PASS = "pass"
    FAIL = "fail"
    BORDERLINE = "borderline"      # measured, but inside the noise of the cut
    REFUSED = "refused"            # we will not judge this: not enough evidence
    SATURATED = "saturated"        # the score stopped discriminating

    @property
    def carries_forward(self) -> bool:
        """Does a design with this verdict go on to the next stage?

        BORDERLINE does, flagged. REFUSED does not: the answer is more
        replicates, not a guess.
        """
        return self in (Verdict.PASS, Verdict.BORDERLINE)


@dataclass(frozen=True)
class Judgement:
    gate: str
    design: str
    verdict: Verdict
    value: float | None
    threshold: float | None
    reason: str
    detail: dict[str, Any] = field(default_factory=dict)
    thresholds_version: str = ""
    config_hash: str = ""

    def as_row(self) -> dict[str, Any]:
        row = asdict(self)
        row["verdict"] = self.verdict.value
        return row

    def __str__(self) -> str:
        v = "n/a" if self.value is None else f"{self.value:.3f}"
        return f"[{self.gate}] {self.design}: {self.verdict.value.upper()} ({v}) - {self.reason}"


class ForbiddenRanking(Exception):
    """Raised when something tries to sort designs by a confidence metric."""


class Gatekeeper:
    def __init__(self, cfg: Config | None = None):
        self.cfg = cfg or load()

    # ---------------------------------------------------------------- helpers

    def _judge(self, gate, design, verdict, value, threshold, reason, **detail):
        return Judgement(
            gate=gate,
            design=design,
            verdict=verdict,
            value=value,
            threshold=threshold,
            reason=reason,
            detail=detail,
            **self.cfg.stamp(),
        )

    # ----------------------------------------------------------------- gate 1

    def gate1(self, design: str, rho: float) -> Judgement:
        """Compactness. rho = Rg(CA) / (2.2 * N^0.38)."""
        cut = float(self.cfg.gate1["rho_max"])
        band = float(self.cfg.gate1["rho_borderline"])
        if rho <= cut - band:
            return self._judge("gate1", design, Verdict.PASS, rho, cut, "compact")
        if rho <= cut + band:
            return self._judge(
                "gate1", design, Verdict.BORDERLINE, rho, cut,
                f"within {band} of the compactness cut",
            )
        reason = "extended" if rho < 1.8 else "rod-like; a straight helix scores ~2.2"
        return self._judge("gate1", design, Verdict.FAIL, rho, cut, reason)

    # ----------------------------------------------------------------- gate 2

    def gate2_sequence(self, design: str, foldback_rmsd: float) -> Judgement:
        """Does this sequence fold back into the backbone it was written for?"""
        cut = float(self.cfg.gate2["foldback_rmsd_max"])
        band = float(self.cfg.gate2["foldback_rmsd_borderline"])
        if foldback_rmsd <= cut - band:
            return self._judge("gate2", design, Verdict.PASS, foldback_rmsd, cut,
                               "folds back onto its design")
        if foldback_rmsd <= cut + band:
            return self._judge("gate2", design, Verdict.BORDERLINE, foldback_rmsd, cut,
                               f"within {band} A of the fold-back cut")
        return self._judge("gate2", design, Verdict.FAIL, foldback_rmsd, cut,
                           "does not adopt the designed shape")

    def gate2_backbone(self, backbone: str, per_sequence: dict[str, float]) -> Judgement:
        """A backbone's verdict is the best of its sequences - and it is only a
        verdict at all if enough sequences were sampled.

        Within-backbone pae spread is 11-20 units, so a pass belongs to the
        sequence. Backbone 87 scored 23.0 on its single sample and is the pool's
        best design: screening one sequence per backbone is how that gets thrown
        away silently.
        """
        need = int(self.cfg.gate2["min_sequences_per_backbone"])
        n = len(per_sequence)
        if n < need:
            return self._judge(
                "gate2", backbone, Verdict.REFUSED, None, None,
                f"only {n} sequence(s) sampled, need {need}; a pass belongs to "
                "the sequence, not the backbone",
                n_sequences=n,
            )
        best_seq = min(per_sequence, key=per_sequence.get)
        j = self.gate2_sequence(backbone, per_sequence[best_seq])
        return self._judge(
            "gate2", backbone, j.verdict, j.value, j.threshold,
            f"{j.reason}; best of {n} sequences ({best_seq})",
            n_sequences=n, best_sequence=best_seq,
        )

    # ----------------------------------------------------------------- gate 3

    def gate3(self, design: str, pae_values: Sequence[float]) -> Judgement:
        """Interface, over replicates. Never from a single draw.

        The score is not deterministic even on one node and the noise is
        heavy-tailed. What is reported is a mean with its spread, and a
        borderline band wide enough to cover that noise.
        """
        need = int(self.cfg.gate3["min_replicates"])
        n = len(pae_values)
        if n < need:
            return self._judge(
                "gate3", design, Verdict.REFUSED, None, None,
                f"{n} replicate(s), need {need}; one draw of this score is not "
                "a verdict (~1 in 10 flips pass/fail)",
                n_replicates=n,
            )

        s = stats.spread(list(pae_values))
        cut = float(self.cfg.gate3["pae_interaction_max"])
        band = float(self.cfg.gate3["pae_borderline"])
        sat = float(self.cfg.gate3["pae_saturation"])
        tight = int(self.cfg.gate3["replicates_for_tight_floor"])
        detail = {
            "n_replicates": n,
            "sd": s.sd,
            "sem": s.sem,
            "lo": s.lo,
            "hi": s.hi,
            "tight_floor": n >= tight,
        }

        if s.mean >= sat:
            return self._judge(
                "gate3", design, Verdict.SATURATED, s.mean, cut,
                f"at or above the saturation point ({sat}); the score no longer "
                "separates designs here",
                **detail,
            )
        if s.mean <= cut - band:
            return self._judge("gate3", design, Verdict.PASS, s.mean, cut,
                               f"below the gate by more than the noise band ({band})",
                               **detail)
        if s.mean <= cut + band:
            return self._judge(
                "gate3", design, Verdict.BORDERLINE, s.mean, cut,
                f"within the noise band ({band}) of the cut; "
                f"{'' if n >= tight else f'{tight} replicates would tighten this'}".strip()
                or "inside the noise of the cut",
                **detail,
            )
        return self._judge("gate3", design, Verdict.FAIL, s.mean, cut,
                           "above the gate", **detail)

    # ------------------------------------------------------------- chirality

    def chirality(self, design: str, n_d_residues: int, n_checked: int,
                  d_residues: Sequence[str] = ()) -> Judgement:
        """Refuse, never repair.

        A D-residue is invisible to Rg, backbone RMSD, pae_interaction and every
        contact metric in this stack, and the force field carries no term that
        would restore it. PDBFixer built five into one loop deterministically;
        addHydrogens inverted a centre in one build in five.
        """
        frac = float(self.cfg.chirality["majority_warn_fraction"])
        if n_checked == 0:
            return self._judge("chirality", design, Verdict.REFUSED, None, None,
                               "no residue had all of N, CA, C, CB to measure")
        if n_d_residues == 0:
            return self._judge("chirality", design, Verdict.PASS, 0.0, 0.0,
                               f"all {n_checked} measurable residues are L",
                               n_checked=n_checked)
        if n_d_residues / n_checked > frac:
            return self._judge(
                "chirality", design, Verdict.REFUSED, float(n_d_residues), 0.0,
                f"{n_d_residues} of {n_checked} residues read as D - the file is "
                "probably mirrored or the convention is backwards; stopping "
                "rather than reporting every residue as D",
                n_checked=n_checked,
            )
        return self._judge(
            "chirality", design, Verdict.FAIL, float(n_d_residues), 0.0,
            f"{n_d_residues} D-residue(s) found; no other check in this pipeline "
            "can see them, and rebuilding is not a fix",
            n_checked=n_checked, d_residues=list(d_residues),
        )

    # --------------------------------------------------------------- ranking

    def check_sort_key(self, key: str) -> None:
        """Raise if `key` is a confidence metric. There is no override."""
        k = key.lower().lstrip("-")
        for forbidden in self.cfg.ranking["forbidden_sort_keys"]:
            if k == forbidden or k.startswith(forbidden + "_") or k.endswith("_" + forbidden):
                raise ForbiddenRanking(
                    f"cannot rank designs by {key!r}: confidence explains ~2.6% of "
                    "fold-back variance pool-wide and inverts on 17 of 47 "
                    "backbones. Rank on fold-back RMSD or grip."
                )

    def rank(self, designs: Sequence[dict], key: str, descending: bool = False) -> list[dict]:
        """Sort gate survivors. Refuses confidence keys; refuses to invent order.

        `43_build_complex.py --top N` slices the alphabetically first N while its
        docstring promises grip order. Here, order is explicit or it is an error.
        """
        self.check_sort_key(key)
        missing = [d.get("design", "?") for d in designs if key not in d]
        if missing:
            raise KeyError(
                f"{len(missing)} design(s) have no {key!r} to rank on "
                f"(e.g. {missing[0]}); refusing to fall back to filename order"
            )
        return sorted(designs, key=lambda d: d[key], reverse=descending)

    # ------------------------------------------------------------ pass rates

    def pass_rate(self, passed: int, total: int) -> stats.Rate:
        return stats.wilson(
            passed, total,
            confidence=float(self.cfg.statistics["confidence_level"]),
            underpowered_below=int(self.cfg.statistics["underpowered_below"]),
        )

    def rank_cells(self, cells: dict[str, tuple[int, int]]) -> list[list[str]]:
        """Group cells into tiers, refusing to order cells whose intervals overlap.

        A cell joins the current tier unless it is separable from that tier's
        LEADER. Comparing each cell only to its immediate neighbour instead
        would chain: adjacent intervals overlap all the way down the six-target
        ladder, so every cell would land in one tier and the gate would say
        nothing at all.

        Note this is deliberately stricter than reading the same data with
        pairwise z-tests. Those suggest three tiers, with pd1 ~20% sitting
        between the ~40% pair and the <=6% group; on Wilson intervals pd1 is not
        separable from pcsk9 and so joins the top tier. The uncorrected ~2 sigma
        that separated them does not survive five adjacent comparisons, so the
        conservative grouping is the one to put in front of a user.

        Returns tiers, best first. Names inside one tier are not ranked against
        each other, ever.
        """
        rates = {name: self.pass_rate(p, t) for name, (p, t) in cells.items()}
        order = sorted(rates, key=lambda n: rates[n].point, reverse=True)
        tiers: list[list[str]] = []
        for name in order:
            if tiers and not stats.comparable(rates[tiers[-1][0]], rates[name]):
                tiers[-1].append(name)
            else:
                tiers.append([name])
        return tiers
