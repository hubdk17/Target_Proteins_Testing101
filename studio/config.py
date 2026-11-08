#!/usr/bin/env python3
"""Load the one threshold file and hash it.

Thresholds used to be awk literals scattered across six scripts. Here they live
in `thresholds.toml` and nothing else is allowed to hold one. Every verdict the
gatekeeper writes carries `config_hash`, so a row can always be traced to the
numbers that judged it.

    from studio.config import load
    cfg = load()
    cfg.gate1["rho_max"]        # 1.3
    cfg.hash                    # 'sha256:9f2c...'
"""
from __future__ import annotations

import hashlib
import tomllib
from dataclasses import dataclass
from pathlib import Path

DEFAULT_PATH = Path(__file__).with_name("thresholds.toml")

# Sections that must exist. A typo in a section name would otherwise surface as
# a KeyError deep inside a verdict, after the GPU time was already spent.
REQUIRED = (
    "gate1",
    "gate2",
    "gate3",
    "ranking",
    "chirality",
    "structure",
    "statistics",
    "cost",
)


@dataclass(frozen=True)
class Config:
    version: str
    hash: str
    path: Path
    gate1: dict
    gate2: dict
    gate3: dict
    ranking: dict
    chirality: dict
    structure: dict
    statistics: dict
    cost: dict

    def stamp(self) -> dict:
        """What goes on every verdict row and into the provenance record."""
        return {"thresholds_version": self.version, "config_hash": self.hash}


def load(path: Path | str = DEFAULT_PATH) -> Config:
    path = Path(path)
    raw = path.read_bytes()
    data = tomllib.loads(raw.decode("utf-8"))

    missing = [s for s in REQUIRED if s not in data]
    if missing:
        raise ValueError(f"{path}: missing threshold sections: {', '.join(missing)}")
    if "version" not in data:
        raise ValueError(f"{path}: no version; bump it on every edit")

    digest = "sha256:" + hashlib.sha256(raw).hexdigest()
    return Config(
        version=data["version"],
        hash=digest,
        path=path,
        **{s: data[s] for s in REQUIRED},
    )
