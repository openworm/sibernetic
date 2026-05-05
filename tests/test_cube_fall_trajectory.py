"""Trajectory-equivalence regression for the demo1_centered cube fall.

Codifies the "looks correct" visual inspection as automated checks.
A cube starting at the boundary geometric center, falling under gravity
in a symmetric box, should:

  1. Fall straight down — minimal horizontal drift
  2. Land on the floor — final cube min_y close to floor (y ≈ 0)
  3. Stay intact — cube y-extent remains a healthy fraction of initial
     (no pancaking, no explosion)
  4. Settle — final cube center y stable (not still falling, not flying up)

Reads a Sibernetic-format position_buffer.txt (or our normalized variant)
and checks all four criteria. Both OpenCL gold-standard output and our
Metal substrate output must pass. If either fails, that's a regression
of physical correctness.

Usage:
    python3 tests/test_cube_fall_trajectory.py path/to/position_buffer.txt
    python3 tests/test_cube_fall_trajectory.py --json path/to/position_buffer.txt

Targets (from demo1_centered config — cube at boundary geo center):
  Initial cube center: (43.59, 44.42, 43.59)
  Initial cube extent: (10.02, 10.02, 10.02)
  Floor at y ≈ 0; expected settled min_y < 5
  Expected horizontal drift: < 5 units (small, since no horizontal force)
  Expected y-extent retention: > 50% (intact, no pancake)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


# Demo1_centered initial cube geometry (from configuration/demo1_centered).
INIT_CENTER = np.array([43.59, 44.42, 43.59])
INIT_EXTENT = np.array([10.02, 10.02, 10.02])

# Pass thresholds — calibrated from healthy reference runs.
THRESHOLDS = {
    'horiz_drift_max': 5.0,        # max |Δxz| of cube center, units
    'final_min_y_max': 5.0,        # cube bottom should be near floor
    'extent_y_ratio_min': 0.5,     # >50% retention = no pancake
    'extent_y_ratio_max': 2.0,     # <200% = no explosion
    'final_settle_window': 0.5,    # last N seconds — cube center y must be stable
    'settle_y_drift_max': 1.0,     # cube center y can wobble by this in settle window
}


def parse_position_buffer(path: Path):
    """Return (n_active, n_boundary, dt, log_step, frames) where each
    frame is (active_xyz, types) numpy arrays."""
    with open(path) as f:
        header = [f.readline().rstrip('\n') for _ in range(11)]
        n_elastic = int(header[6])
        n_liquid = int(header[7])
        n_boundary = int(header[8])
        dt = float(header[9])
        log_step = int(header[10])
        n_active = n_elastic + n_liquid

        # First frame: all particles
        first_active_xyz, first_active_t = [], []
        for _ in range(n_active + n_boundary):
            ws = f.readline().split()
            x, y, z, t = float(ws[0]), float(ws[1]), float(ws[2]), float(ws[3])
            if t < 3.0:
                first_active_xyz.append([x, y, z])
                first_active_t.append(t)

        frames = [(np.array(first_active_xyz), np.array(first_active_t))]

        # Subsequent frames — auto-detect format:
        #   Sibernetic: only active per frame (n_active lines)
        #   Our normalized dump: all particles per frame (n_total lines)
        rest = [ln for ln in f if ln.strip()]
        # Count particles in the next "frame" by looking for a boundary line
        # (type >= 3.0) — if first frame has all-active, format is Sibernetic.
        if rest:
            chunk_active = n_active
            chunk_full = n_active + n_boundary
            # Try active-only first
            if len(rest) % chunk_active == 0:
                chunk = chunk_active
            elif len(rest) % chunk_full == 0:
                chunk = chunk_full
            else:
                chunk = chunk_active

            for i in range(0, len(rest), chunk):
                xyz, ts = [], []
                for ln in rest[i:i + chunk]:
                    ws = ln.split()
                    t = float(ws[3])
                    if t < 3.0:
                        xyz.append([float(ws[0]), float(ws[1]), float(ws[2])])
                        ts.append(t)
                if len(xyz) == n_active:
                    frames.append((np.array(xyz), np.array(ts)))

    return n_active, n_boundary, dt, log_step, frames


def cube_metrics_per_frame(frames):
    """Per-frame cube center and extent of elastic particles (type 2.x)."""
    centers, extents = [], []
    for xyz, ts in frames:
        elastic = xyz[(ts >= 2.0) & (ts < 3.0)]
        centers.append(elastic.mean(axis=0))
        extents.append(elastic.max(axis=0) - elastic.min(axis=0))
    return np.array(centers), np.array(extents)


def check_trajectory(path: Path, *, init_center=INIT_CENTER,
                     init_extent=INIT_EXTENT, thresholds=THRESHOLDS):
    n_active, n_boundary, dt, log_step, frames = parse_position_buffer(path)
    centers, extents = cube_metrics_per_frame(frames)

    n_frames = len(frames)
    sim_t = np.arange(n_frames) * log_step * dt
    total_time = sim_t[-1] if n_frames > 0 else 0.0

    # === Check 1: initial cube position matches expected (within 1 unit) ===
    init_err = float(np.linalg.norm(centers[0] - init_center))
    init_ok = init_err < 1.0

    # === Check 2: horizontal drift over the whole trajectory ===
    drift_xz = np.sqrt((centers[:, 0] - init_center[0]) ** 2
                     + (centers[:, 2] - init_center[2]) ** 2)
    max_drift = float(drift_xz.max())
    drift_ok = max_drift < thresholds['horiz_drift_max']

    # === Check 3: cube reached the floor (final min_y near 0) ===
    final_min_y = float(centers[-1, 1] - extents[-1, 1] / 2)
    floor_ok = final_min_y < thresholds['final_min_y_max']

    # === Check 4: cube didn't pancake or explode (y-extent retention) ===
    final_ext_y = float(extents[-1, 1])
    ext_ratio = final_ext_y / init_extent[1]
    ext_ok = (thresholds['extent_y_ratio_min'] < ext_ratio
              < thresholds['extent_y_ratio_max'])

    # === Check 5: cube has settled (center y stable in last N seconds) ===
    settle_window_t = thresholds['final_settle_window']
    settle_idx = sim_t >= max(0, total_time - settle_window_t)
    if settle_idx.sum() >= 2:
        y_centers_settle = centers[settle_idx, 1]
        settle_drift = float(y_centers_settle.max() - y_centers_settle.min())
        settle_ok = settle_drift < thresholds['settle_y_drift_max']
    else:
        settle_drift = 0.0
        settle_ok = True

    all_ok = init_ok and drift_ok and floor_ok and ext_ok and settle_ok

    return {
        'path': str(path),
        'n_frames': n_frames,
        'sim_time_s': float(total_time),
        'init_center': centers[0].tolist(),
        'init_err': init_err,
        'init_ok': init_ok,
        'horiz_drift_max': max_drift,
        'drift_ok': drift_ok,
        'final_min_y': final_min_y,
        'floor_ok': floor_ok,
        'extent_y_ratio': ext_ratio,
        'ext_ok': ext_ok,
        'settle_drift': settle_drift,
        'settle_ok': settle_ok,
        'all_ok': all_ok,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('path', type=Path)
    ap.add_argument('--json', action='store_true')
    args = ap.parse_args()

    r = check_trajectory(args.path)

    if args.json:
        print(json.dumps(r, indent=2))
        return 0 if r['all_ok'] else 1

    print(f"Cube-fall trajectory check: {r['path']}")
    print(f"  Frames: {r['n_frames']}  sim_t: 0 → {r['sim_time_s']:.2f}s")
    print()
    chk = lambda ok: '[PASS]' if ok else '[FAIL]'
    print(f"  {chk(r['init_ok'])} initial center {r['init_center']}  "
          f"(error vs target = {r['init_err']:.3f})")
    print(f"  {chk(r['drift_ok'])} horiz drift max = {r['horiz_drift_max']:.3f} "
          f"(threshold {THRESHOLDS['horiz_drift_max']})")
    print(f"  {chk(r['floor_ok'])} final min_y = {r['final_min_y']:.3f} "
          f"(threshold < {THRESHOLDS['final_min_y_max']})")
    print(f"  {chk(r['ext_ok'])} y-extent retention = {r['extent_y_ratio']:.3f} "
          f"(must be in [{THRESHOLDS['extent_y_ratio_min']}, "
          f"{THRESHOLDS['extent_y_ratio_max']}])")
    print(f"  {chk(r['settle_ok'])} settle drift = {r['settle_drift']:.3f} "
          f"(threshold < {THRESHOLDS['settle_y_drift_max']})")
    print()
    if r['all_ok']:
        print("  [OVERALL PASS] cube falls, lands intact, settles, no slide")
    else:
        print("  [OVERALL FAIL] one or more trajectory checks failed")
    return 0 if r['all_ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
