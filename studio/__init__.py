"""studio - the hosted layer over the numbered pipeline scripts.

The scripts in the repository root stay the source of truth and are invoked as
they are, so a hosted result and a cluster result remain comparable. This
package adds only the things a hosted service needs and a hand-run pipeline does
without: run isolation, one place for thresholds, one writer of verdicts,
provenance, and honest failure reporting.

Plan: /mnt/project-files/studio/DENOVO-STUDIO-PLAN-2026-11-08.md
"""

__all__ = [
    "checks",
    "config",
    "gatekeeper",
    "provenance",
    "residues",
    "scores",
    "stats",
    "store",
]
