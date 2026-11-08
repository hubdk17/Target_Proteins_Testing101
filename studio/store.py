#!/usr/bin/env python3
"""Run state in SQLite, keyed by a run UUID.

This replaces append-only `.sc` text files and target-named directories, and it
is what structurally removes three findings from the 2026-11-07 code review
rather than patching six scripts:

- `TAG` never reached `gate1/`, `mpnn_out_`, `sc_in/`, `gate2/`, `cplx_in/`, so
  two tags on one target overwrote each other's gate lists and mixed their
  pools. Here the run id is the directory, so two runs cannot collide.
- `41_gate2_report.sh` never truncated `gate2/<name>.pass`, so a run where
  nothing passed printed the PREVIOUS run's survivors as its own. Here
  survivors are rows belonging to a run id; there is nothing to go stale.
- `thread_multi.py` did not clear `sc_in/<name>`, so old threaded designs were
  re-predicted into the new score file. Same fix.

One more rule from the project, enforced here: a stage's success is its output
count and elapsed time, never its exit code or the provider's status.
"""
from __future__ import annotations

import json
import sqlite3
import time
import uuid
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    run_id        TEXT PRIMARY KEY,
    target        TEXT NOT NULL,
    created_at    REAL NOT NULL,
    state         TEXT NOT NULL,
    thresholds_version TEXT NOT NULL,
    config_hash   TEXT NOT NULL,
    cost_ceiling_usd REAL NOT NULL,
    spec_json     TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS designs (
    run_id   TEXT NOT NULL,
    design   TEXT NOT NULL,
    stage    TEXT NOT NULL,
    backbone TEXT,
    path     TEXT,
    created_at REAL NOT NULL,
    PRIMARY KEY (run_id, design, stage)
);
CREATE TABLE IF NOT EXISTS judgements (
    run_id    TEXT NOT NULL,
    gate      TEXT NOT NULL,
    design    TEXT NOT NULL,
    verdict   TEXT NOT NULL,
    value     REAL,
    threshold REAL,
    reason    TEXT NOT NULL,
    detail_json TEXT NOT NULL,
    thresholds_version TEXT NOT NULL,
    config_hash TEXT NOT NULL,
    created_at REAL NOT NULL,
    PRIMARY KEY (run_id, gate, design)
);
CREATE TABLE IF NOT EXISTS events (
    run_id  TEXT NOT NULL,
    at      REAL NOT NULL,
    stage   TEXT,
    kind    TEXT NOT NULL,
    message TEXT NOT NULL,
    detail_json TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS spend (
    run_id   TEXT NOT NULL,
    at       REAL NOT NULL,
    stage    TEXT NOT NULL,
    gpu_hours REAL NOT NULL,
    usd      REAL NOT NULL,
    note     TEXT NOT NULL
);
"""

# Event kinds. `silent_failure` is the important one: a stage that reported
# success and produced too few artifacts. It stops the run and is never
# reported as a result.
EVENT_KINDS = (
    "queued", "started", "progress", "finished",
    "hard_failure", "silent_failure", "stall",
    "refusal", "ceiling_reached", "cancelled",
)

TERMINAL_STATES = ("finished", "failed", "stopped", "cancelled")


@dataclass
class Run:
    run_id: str
    target: str
    root: Path
    cost_ceiling_usd: float


class Store:
    def __init__(self, root: Path | str):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.db_path = self.root / "studio.sqlite"
        with closing(self._conn()) as c:
            c.executescript(SCHEMA)
            c.commit()

    def _conn(self) -> sqlite3.Connection:
        c = sqlite3.connect(self.db_path)
        c.row_factory = sqlite3.Row
        return c

    # ------------------------------------------------------------------- runs

    def create_run(self, target: str, spec: dict, config_stamp: dict,
                   cost_ceiling_usd: float) -> Run:
        """A new run. Its id is its directory; nothing is keyed by target name."""
        run_id = uuid.uuid4().hex
        run_root = self.root / "runs" / run_id
        for sub in ("inputs", "stage1", "stage2", "gate2", "stage3", "artifacts"):
            (run_root / sub).mkdir(parents=True, exist_ok=True)
        with closing(self._conn()) as c:
            c.execute(
                "INSERT INTO runs (run_id,target,created_at,state,"
                "thresholds_version,config_hash,cost_ceiling_usd,spec_json) "
                "VALUES (?,?,?,?,?,?,?,?)",
                (run_id, target, time.time(), "created",
                 config_stamp["thresholds_version"], config_stamp["config_hash"],
                 cost_ceiling_usd, json.dumps(spec, sort_keys=True)),
            )
            c.commit()
        return Run(run_id, target, run_root, cost_ceiling_usd)

    def set_state(self, run_id: str, state: str) -> None:
        with closing(self._conn()) as c:
            c.execute("UPDATE runs SET state=? WHERE run_id=?", (state, run_id))
            c.commit()

    def get_run(self, run_id: str) -> dict[str, Any]:
        with closing(self._conn()) as c:
            row = c.execute("SELECT * FROM runs WHERE run_id=?", (run_id,)).fetchone()
        if row is None:
            raise KeyError(f"no run {run_id}")
        d = dict(row)
        d["spec"] = json.loads(d.pop("spec_json"))
        return d

    # ---------------------------------------------------------------- designs

    def add_designs(self, run_id: str, stage: str,
                    designs: Iterable[tuple[str, str | None, str | None]]) -> int:
        rows = [(run_id, name, stage, backbone, path, time.time())
                for name, backbone, path in designs]
        with closing(self._conn()) as c:
            c.executemany(
                "INSERT OR REPLACE INTO designs "
                "(run_id,design,stage,backbone,path,created_at) VALUES (?,?,?,?,?,?)",
                rows,
            )
            c.commit()
        return len(rows)

    def count_designs(self, run_id: str, stage: str) -> int:
        with closing(self._conn()) as c:
            return c.execute(
                "SELECT COUNT(*) FROM designs WHERE run_id=? AND stage=?",
                (run_id, stage),
            ).fetchone()[0]

    # ------------------------------------------------------------- judgements

    def record(self, run_id: str, judgements: Iterable) -> int:
        rows = []
        for j in judgements:
            r = j.as_row()
            rows.append((
                run_id, r["gate"], r["design"], r["verdict"], r["value"],
                r["threshold"], r["reason"], json.dumps(r["detail"], sort_keys=True),
                r["thresholds_version"], r["config_hash"], time.time(),
            ))
        with closing(self._conn()) as c:
            c.executemany(
                "INSERT OR REPLACE INTO judgements (run_id,gate,design,verdict,"
                "value,threshold,reason,detail_json,thresholds_version,"
                "config_hash,created_at) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                rows,
            )
            c.commit()
        return len(rows)

    def survivors(self, run_id: str, gate: str) -> list[str]:
        """Designs that carry forward from `gate`, for this run only.

        PASS and BORDERLINE carry forward; REFUSED does not, because the answer
        to "not enough evidence" is more replicates, not a guess.
        """
        with closing(self._conn()) as c:
            rows = c.execute(
                "SELECT design FROM judgements WHERE run_id=? AND gate=? "
                "AND verdict IN ('pass','borderline') ORDER BY design",
                (run_id, gate),
            ).fetchall()
        return [r["design"] for r in rows]

    def funnel(self, run_id: str) -> dict[str, dict[str, int]]:
        """Per gate, a count of each verdict. This is the user's funnel view."""
        with closing(self._conn()) as c:
            rows = c.execute(
                "SELECT gate, verdict, COUNT(*) AS n FROM judgements "
                "WHERE run_id=? GROUP BY gate, verdict", (run_id,),
            ).fetchall()
        out: dict[str, dict[str, int]] = {}
        for r in rows:
            out.setdefault(r["gate"], {})[r["verdict"]] = r["n"]
        return out

    # ------------------------------------------------------------ events/spend

    def event(self, run_id: str, kind: str, message: str,
              stage: str | None = None, **detail) -> None:
        if kind not in EVENT_KINDS:
            raise ValueError(f"unknown event kind {kind!r}; known: {EVENT_KINDS}")
        with closing(self._conn()) as c:
            c.execute(
                "INSERT INTO events (run_id,at,stage,kind,message,detail_json) "
                "VALUES (?,?,?,?,?,?)",
                (run_id, time.time(), stage, kind, message,
                 json.dumps(detail, sort_keys=True)),
            )
            c.commit()

    def events(self, run_id: str) -> list[dict[str, Any]]:
        with closing(self._conn()) as c:
            rows = c.execute(
                "SELECT * FROM events WHERE run_id=? ORDER BY at", (run_id,)
            ).fetchall()
        return [dict(r) for r in rows]

    def charge(self, run_id: str, stage: str, gpu_hours: float, usd: float,
               note: str = "") -> float:
        """Record spend and return the run's total. The caller compares it to
        the ceiling: on rented GPU a stalled job keeps billing."""
        with closing(self._conn()) as c:
            c.execute(
                "INSERT INTO spend (run_id,at,stage,gpu_hours,usd,note) "
                "VALUES (?,?,?,?,?,?)",
                (run_id, time.time(), stage, gpu_hours, usd, note),
            )
            c.commit()
        return self.spent(run_id)

    def spent(self, run_id: str) -> float:
        with closing(self._conn()) as c:
            v = c.execute(
                "SELECT COALESCE(SUM(usd),0) FROM spend WHERE run_id=?", (run_id,)
            ).fetchone()[0]
        return float(v)
