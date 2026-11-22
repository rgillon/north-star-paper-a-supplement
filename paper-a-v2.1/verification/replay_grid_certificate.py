#!/usr/bin/env python3
"""Replay the complete 24-grid terminal certificate with a JSON receipt.

Run from the directory containing ``experiments`` (or supply its path):

    python replay_grid_certificate.py --experiments experiments \
        --checkpoint grid_replay_checkpoint.json > grid_replay.json

Standard output is one JSON document in the same schema as the recorded
Paper A replay. Progress is printed as JSON lines on standard error.
The optional checkpoint is replaced atomically after each completed grid.
No numerical prime tuples are sampled; every grid uses symbolic exclusion.
"""
from __future__ import annotations

import argparse
import json
from itertools import permutations
from pathlib import Path
import platform
import sys
import time

SOURCE_COMMIT = "20dd528d5fdd4bf7a16cc5ce05f9cd25ae155c25"
BASE_SHAPES = ((2,2,1,1), (3,1,1,1), (2,1,1,1),
               (4,1,1), (3,2,1), (2,2,2))
EXCLUSION_KEYS = ("cube", "incidence", "root_factor",
                  "row_factor", "largest_prime")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--experiments", type=Path, default=Path("experiments"),
                        help="directory containing the pinned verifier and dependencies; default: ./experiments")
    parser.add_argument("--checkpoint", type=Path,
                        help="optional JSON file updated after every completed grid")
    args = parser.parse_args()
    experiments = args.experiments.resolve()
    if not (experiments / "level1_grid_gap_certificate.py").is_file():
        parser.error(f"verifier not found in {experiments}")
    sys.path.insert(0, str(experiments))

    import numpy as np
    import sympy as sp
    from level1_grid_gap_certificate import certify

    shapes = sorted({ordered for shape in BASE_SHAPES
                     for ordered in permutations(shape)})
    assert len(shapes) == 24
    report = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "sympy": sp.__version__,
        "source_commit": SOURCE_COMMIT,
        "rows": [],
        "complete": False,
    }
    start = time.perf_counter()

    def checkpoint() -> None:
        report["elapsed_seconds"] = time.perf_counter() - start
        if args.checkpoint is not None:
            args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
            temporary = args.checkpoint.with_name(args.checkpoint.name + ".tmp")
            temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            temporary.replace(args.checkpoint)

    try:
        for shape in shapes:
            row_start = time.perf_counter()
            counts, residual = certify(shape)
            row = {"shape": shape, "counts": dict(counts),
                   "residual": residual,
                   "seconds": time.perf_counter() - row_start}
            report["rows"].append(row)
            checkpoint()
            print(json.dumps(row), file=sys.stderr, flush=True)
            if residual:
                raise AssertionError(f"unresolved patterns in {shape}: {residual}")
            exclusions = sum(counts.get(key, 0) for key in EXCLUSION_KEYS)
            if exclusions != counts.get("root", 0):
                raise AssertionError(f"exclusion counts do not balance in {shape}")
        report["complete"] = True
    finally:
        checkpoint()
        print(json.dumps(report, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
