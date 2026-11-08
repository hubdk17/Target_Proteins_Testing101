#!/usr/bin/env python3
"""The one-page record a user can put in a report, and we can audit a result from.

What it has to pin, and why each item is here rather than nice-to-have:

- the threshold file's version and hash, so a verdict can be traced to the
  numbers that judged it;
- the container image digest and every weight file's checksum, because a
  silently different parameter set is exactly the kind of fault that produces
  plausible numbers;
- every seed, and the node or pod the work ran on, because the scoring gate is
  not deterministic even on one node;
- what ARRIVED, not what was requested: the hotspots actually used, the contig
  actually used, counts of files actually produced. Stage 1 never logs its
  hotspots or potentials, so a run's own record is the only place this exists;
- every refusal that fired, so a quiet run and a run with nine suppressed
  warnings do not look the same.
"""
from __future__ import annotations

import hashlib
import json
import platform
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any


def file_digest(path: Path | str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return "sha256:" + h.hexdigest()


@dataclass
class Provenance:
    run_id: str
    target: str
    created_at: float = field(default_factory=time.time)

    # judged by
    thresholds_version: str = ""
    config_hash: str = ""

    # ran on
    image_digest: str = ""
    host: str = field(default_factory=platform.node)
    gpu: str = ""
    provider: str = ""
    weights: dict[str, str] = field(default_factory=dict)
    versions: dict[str, str] = field(default_factory=dict)

    # what actually arrived
    inputs: dict[str, Any] = field(default_factory=dict)
    stage_params: dict[str, Any] = field(default_factory=dict)
    seeds: dict[str, Any] = field(default_factory=dict)
    artifact_counts: dict[str, int] = field(default_factory=dict)

    # what we refused to do
    refusals: list[dict[str, Any]] = field(default_factory=list)
    gpu_hours: float = 0.0
    usd: float = 0.0

    def add_weight(self, name: str, path: Path | str) -> None:
        self.weights[name] = file_digest(path)

    def refuse(self, where: str, why: str, **detail) -> None:
        self.refusals.append({"where": where, "why": why, "at": time.time(), **detail})

    def missing(self) -> list[str]:
        """Fields that must be filled before a result is shown to a user."""
        gaps = []
        if not self.config_hash:
            gaps.append("config_hash")
        if not self.image_digest:
            gaps.append("image_digest")
        if not self.weights:
            gaps.append("weights")
        if not self.seeds:
            gaps.append("seeds")
        if not self.artifact_counts:
            gaps.append("artifact_counts")
        return gaps

    def write(self, path: Path | str, require_complete: bool = True) -> Path:
        """Write the record. Refuses an incomplete one by default.

        A result with no provenance is a number with no origin, and this project
        has spent real time working out after the fact which inputs a run used.
        """
        gaps = self.missing()
        if require_complete and gaps:
            raise ValueError(
                "provenance is incomplete: " + ", ".join(gaps) +
                ". Fill these or pass require_complete=False to write a draft."
            )
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        record = asdict(self)
        record["complete"] = not gaps
        record["missing"] = gaps
        path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        return path
