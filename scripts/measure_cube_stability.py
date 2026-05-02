#!/usr/bin/env python3
"""Measure structural stability of an elastic cube after floor collision.

Reads a Sibernetic position_buffer.txt produced by demo1 (or any config
with elastic particles) and computes the vertical extent of the elastic
particles in the first vs last logged frame. The ratio (final / initial)
tells you whether the cube maintained its structural integrity after
hitting the floor or collapsed into a pancake.

Healthy cube  : extent retention > 0.5  (cube settled but kept its shape)
Pancake bug   : extent retention < 0.2  (cube flattened against the floor)

Background: the README claims Taichi's elastic body flattens to mean Y
~0.24 after a 3s floor collision (vs PyTorch's 1.45). The proximate cause
is the well-documented coordinate-space bug — elastic forces are
effectively ~287x weaker than they should be relative to SPH pressure
because they're computed in scaled coordinates while pressure is in world
coordinates. This script gives that visual symptom an objective metric so
fixes can be measured without watching gifs.

Position buffer format (single header at top, then concatenated frames):

    line  0..5   simulation box bounds (xmin, xmax, ymin, ymax, zmin, zmax)
    line  6      numOfElasticP
    line  7      numOfLiquidP
    line  8      numOfBoundaryP
    line  9      time_step
    line 10      log_step
    line 11+     particles, all frames concatenated, in elastic→liquid→
                 boundary order. Each line: "x y z type", with type ~= 2
                 for elastic, ~= 1 for liquid, == 3 for boundary.

Usage:

    python3 scripts/measure_cube_stability.py path/to/position_buffer.txt
    python3 scripts/measure_cube_stability.py path/to/position_buffer.txt --json
    python3 scripts/measure_cube_stability.py path/to/position_buffer.txt --threshold 0.5
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def parse_header(lines):
    """Return (n_elastic, n_liquid, n_boundary, time_step, log_step)."""
    return (
        int(lines[6].split()[0]),
        int(lines[7].split()[0]),
        int(lines[8].split()[0]),
        float(lines[9].split()[0]),
        int(lines[10].split()[0]),
    )


def measure(position_file: Path):
    """Measure first- and last-frame elastic Y stats from a position buffer.

    Streams the file twice — once to count frames and grab frame 0, once
    to grab the final frame — so it runs in O(frame_size) memory rather
    than loading the entire buffer (which can be 100s of MB).
    """
    with position_file.open() as f:
        head = [next(f) for _ in range(11)]
    n_e, n_l, n_b, dt, log_step = parse_header(head)
    frame_size = n_e + n_l + n_b
    if frame_size == 0:
        raise ValueError("buffer header reports 0 particles")

    body_size = position_file.stat().st_size
    # We don't know the exact byte length per particle line, so count lines
    # rather than seeking. One pass is unavoidable to find the last frame.
    total_lines = 0
    with position_file.open() as f:
        for _ in f:
            total_lines += 1
    n_frames = (total_lines - 11) // frame_size
    if n_frames < 1:
        raise ValueError(
            f"buffer has fewer than 1 complete frame "
            f"({total_lines - 11} particle lines / {frame_size} per frame)"
        )

    # Pass 1: capture first frame's elastic Y values (lines 11 .. 11+n_e).
    first_y = []
    with position_file.open() as f:
        for _ in range(11):
            next(f)
        for _ in range(n_e):
            parts = next(f).split()
            first_y.append(float(parts[1]))

    # Pass 2: skip to the last frame's elastic block.
    last_frame_start = 11 + (n_frames - 1) * frame_size
    last_y = []
    with position_file.open() as f:
        for _ in range(last_frame_start):
            next(f)
        for _ in range(n_e):
            parts = next(f).split()
            last_y.append(float(parts[1]))

    initial_extent = max(first_y) - min(first_y)
    final_extent = max(last_y) - min(last_y)
    initial_mean_y = sum(first_y) / len(first_y)
    final_mean_y = sum(last_y) / len(last_y)
    initial_min_y = min(first_y)
    final_min_y = min(last_y)
    retention = final_extent / initial_extent if initial_extent > 0 else 0.0

    return {
        "position_file": str(position_file),
        "buffer_size_bytes": body_size,
        "n_elastic": n_e,
        "n_liquid": n_l,
        "n_boundary": n_b,
        "n_frames": n_frames,
        "log_step": log_step,
        "buffer_header_dt": dt,
        "initial_elastic_extent_y": initial_extent,
        "final_elastic_extent_y": final_extent,
        "extent_retention": retention,
        "initial_elastic_mean_y": initial_mean_y,
        "final_elastic_mean_y": final_mean_y,
        "initial_elastic_min_y": initial_min_y,
        "final_elastic_min_y": final_min_y,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("position_file", type=Path)
    p.add_argument(
        "--threshold", type=float, default=0.5,
        help="extent_retention below this is treated as a pancake failure (default 0.5)",
    )
    p.add_argument(
        "--json", action="store_true",
        help="emit metrics as JSON (in addition to the human summary)",
    )
    args = p.parse_args(argv)

    if not args.position_file.exists():
        print(f"error: {args.position_file} does not exist", file=sys.stderr)
        return 2

    m = measure(args.position_file)
    passed = m["extent_retention"] >= args.threshold
    status = "PASS" if passed else "FAIL"

    print(f"=== Cube stability: {status} ===")
    print(f"  file              : {m['position_file']} ({m['buffer_size_bytes']/1024/1024:.1f} MB)")
    print(f"  particles         : {m['n_elastic']} elastic + {m['n_liquid']} liquid + {m['n_boundary']} boundary")
    print(f"  frames            : {m['n_frames']} (logstep={m['log_step']}, buffer-header dt={m['buffer_header_dt']:g})")
    print()
    print(f"  elastic Y extent  : {m['initial_elastic_extent_y']:7.3f} → {m['final_elastic_extent_y']:7.3f}  (retention {m['extent_retention']:.1%})")
    print(f"  elastic mean Y    : {m['initial_elastic_mean_y']:7.3f} → {m['final_elastic_mean_y']:7.3f}")
    print(f"  elastic lowest Y  : {m['initial_elastic_min_y']:7.3f} → {m['final_elastic_min_y']:7.3f}")
    print()
    print(f"  threshold         : retention >= {args.threshold:.0%}")
    print(f"  result            : {status}")
    if not passed:
        print(f"  → cube collapsed; structural integrity not preserved.")

    if args.json:
        print()
        print(json.dumps({**m, "threshold": args.threshold, "passed": passed}, indent=2))

    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
