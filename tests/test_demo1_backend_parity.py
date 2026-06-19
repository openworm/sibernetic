"""Backend equivalence regression for the unmodified demo1 cube fall.

Two trajectories of the SAME scenario (`configuration/demo1`) — one from
the OpenCL gold-standard, one from the Metal native substrate — should
produce equivalent particle motion within tolerance. "Equivalent" here
means the substrates agree, not that the trajectory looks textbook.
demo1's boundary is asymmetric (3459 -x face vs 3671 +x face particles),
so a cube placed at the box center sits ~0.83 units off the boundary's
inner geometric center; both backends slide as a result. They should
slide the *same way* — that's the equivalence claim.

Reads two Sibernetic-format position_buffer.txt files (variable or
fixed-frame), samples matched sim times, computes:

  1. Initial-position agreement (sanity — same starting state)
  2. Per-particle L2 trajectory error vs OpenCL reference
  3. Cube-center trajectory diff (xyz, over time)
  4. Cube-extent trajectory diff (xyz, over time)
  5. Final-state agreement

Tolerances are loose at first run — we calibrate from the first OpenCL +
Metal pair we capture, then tighten as the substrates converge. See
the printed metrics on first run; bake them into THRESHOLDS once stable.

Usage:
    python3 tests/test_demo1_backend_parity.py REFERENCE.txt CANDIDATE.txt
    python3 tests/test_demo1_backend_parity.py --json REF.txt CAND.txt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


# Loose initial tolerances. Tighten after first calibration run.
THRESHOLDS = {
    "init_l2_max": 0.5,  # initial state should match closely
    "mean_per_particle_l2_max": 8.0,  # average over time of mean per-particle L2
    "cube_center_diff_max": 8.0,  # max ||Δcube_center|| over trajectory
    "cube_extent_diff_max": 4.0,  # max ||Δextent|| (componentwise sum)
    "final_l2_max": 12.0,  # final-state per-particle L2
}


def parse_position_buffer(path: Path):
    """Return (n_active, n_boundary, dt, log_step, frames, frame_times).

    frames: list of (active_xyz [n_active, 3], types [n_active]) — boundary
    discarded; only non-boundary particles compared. Auto-detects whether
    subsequent frames carry boundary or not.

    frame_times: list of explicit per-frame sim times (seconds) if a
    `<path>.times.txt` sidecar exists (handles variable chunk sizes).
    Otherwise None — caller derives times from frame index * dt * log_step.
    """
    with open(path) as f:
        header = [f.readline().rstrip("\n") for _ in range(11)]
        n_elastic = int(header[6])
        n_liquid = int(header[7])
        n_boundary = int(header[8])
        dt = float(header[9])
        log_step = int(header[10])
        n_active = n_elastic + n_liquid

        # Frame 0 always has all particles in Sibernetic format.
        first_active_xyz, first_active_t = [], []
        for _ in range(n_active + n_boundary):
            ws = f.readline().split()
            x, y, z, t = float(ws[0]), float(ws[1]), float(ws[2]), float(ws[3])
            if t < 3.0:
                first_active_xyz.append([x, y, z])
                first_active_t.append(t)

        frames = [(np.array(first_active_xyz), np.array(first_active_t))]

        rest = [ln for ln in f if ln.strip()]
        if rest:
            chunk_active = n_active
            chunk_full = n_active + n_boundary
            if len(rest) % chunk_active == 0:
                chunk = chunk_active
            elif len(rest) % chunk_full == 0:
                chunk = chunk_full
            else:
                chunk = chunk_active
            for i in range(0, len(rest), chunk):
                xyz, ts = [], []
                for ln in rest[i : i + chunk]:
                    ws = ln.split()
                    t = float(ws[3])
                    if t < 3.0:
                        xyz.append([float(ws[0]), float(ws[1]), float(ws[2])])
                        ts.append(t)
                if len(xyz) == n_active:
                    frames.append((np.array(xyz), np.array(ts)))

    sidecar = path.with_suffix(path.suffix + ".times.txt")
    frame_times = None
    if sidecar.exists():
        with open(sidecar) as f:
            frame_times = [float(ln.strip()) for ln in f if ln.strip()]
        if len(frame_times) != len(frames):
            print(
                f"WARN: sidecar has {len(frame_times)} times but {len(frames)} "
                f"frames in {path.name} — falling back to uniform timing"
            )
            frame_times = None

    return n_active, n_boundary, dt, log_step, frames, frame_times


def cube_metrics(xyz, types):
    """Cube center and extent for elastic particles in this frame."""
    elastic = xyz[(types >= 2.0) & (types < 3.0)]
    if len(elastic) == 0:
        return np.zeros(3), np.zeros(3)
    return elastic.mean(axis=0), elastic.max(axis=0) - elastic.min(axis=0)


def sample_at_times(frames, dt, log_step, target_times, frame_times=None):
    """Pick frames closest to each target sim time. Uses explicit
    `frame_times` (one float per frame) if provided, else falls back to
    uniform `frame_dt = dt * log_step`."""
    if frame_times is not None:
        ft = np.array(frame_times)
        out = []
        for t in target_times:
            idx = int(np.argmin(np.abs(ft - t)))
            out.append(frames[idx])
        return out
    frame_dt = dt * log_step
    out = []
    for t in target_times:
        idx = min(int(round(t / frame_dt)), len(frames) - 1)
        out.append(frames[idx])
    return out


def check_parity(
    ref_path: Path, cand_path: Path, *, n_samples: int = 30, thresholds=THRESHOLDS
):
    r_n, _, r_dt, r_ls, r_frames, r_times = parse_position_buffer(ref_path)
    c_n, _, c_dt, c_ls, c_frames, c_times = parse_position_buffer(cand_path)

    if r_n != c_n:
        return {
            "error": f"particle-count mismatch: ref={r_n} cand={c_n}",
            "all_ok": False,
        }

    r_t_max = r_times[-1] if r_times else (len(r_frames) - 1) * r_dt * r_ls
    c_t_max = c_times[-1] if c_times else (len(c_frames) - 1) * c_dt * c_ls
    t_max = min(r_t_max, c_t_max)
    target_times = np.linspace(0, t_max, n_samples)

    r_samples = sample_at_times(r_frames, r_dt, r_ls, target_times, r_times)
    c_samples = sample_at_times(c_frames, c_dt, c_ls, target_times, c_times)

    per_particle_l2 = []  # mean L2 per frame
    cube_center_diffs = []  # ||Δcube_center||
    cube_extent_diffs = []  # ||Δextent||
    for (r_xyz, r_t), (c_xyz, c_t) in zip(r_samples, c_samples):
        diff = r_xyz - c_xyz
        l2 = np.sqrt((diff**2).sum(axis=1))
        per_particle_l2.append(float(l2.mean()))

        r_cc, r_ext = cube_metrics(r_xyz, r_t)
        c_cc, c_ext = cube_metrics(c_xyz, c_t)
        cube_center_diffs.append(float(np.linalg.norm(r_cc - c_cc)))
        cube_extent_diffs.append(float(np.linalg.norm(r_ext - c_ext)))

    init_l2 = per_particle_l2[0]
    final_l2 = per_particle_l2[-1]
    mean_l2 = float(np.mean(per_particle_l2))
    max_cube_center_diff = float(np.max(cube_center_diffs))
    max_cube_extent_diff = float(np.max(cube_extent_diffs))

    init_ok = init_l2 < thresholds["init_l2_max"]
    mean_ok = mean_l2 < thresholds["mean_per_particle_l2_max"]
    center_ok = max_cube_center_diff < thresholds["cube_center_diff_max"]
    extent_ok = max_cube_extent_diff < thresholds["cube_extent_diff_max"]
    final_ok = final_l2 < thresholds["final_l2_max"]
    all_ok = init_ok and mean_ok and center_ok and extent_ok and final_ok

    return {
        "ref_path": str(ref_path),
        "cand_path": str(cand_path),
        "n_active_particles": r_n,
        "sim_time_s": float(t_max),
        "n_samples": n_samples,
        "init_l2": init_l2,
        "init_ok": init_ok,
        "mean_per_particle_l2": mean_l2,
        "mean_ok": mean_ok,
        "max_cube_center_diff": max_cube_center_diff,
        "center_ok": center_ok,
        "max_cube_extent_diff": max_cube_extent_diff,
        "extent_ok": extent_ok,
        "final_l2": final_l2,
        "final_ok": final_ok,
        "all_ok": all_ok,
        "per_sample_per_particle_l2": per_particle_l2,
        "per_sample_cube_center_diff": cube_center_diffs,
        "per_sample_cube_extent_diff": cube_extent_diffs,
        "sample_times_s": target_times.tolist(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "reference", type=Path, help="OpenCL gold-standard position_buffer.txt"
    )
    ap.add_argument("candidate", type=Path, help="Metal substrate position_buffer.txt")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--samples", type=int, default=30)
    args = ap.parse_args()

    r = check_parity(args.reference, args.candidate, n_samples=args.samples)

    if args.json:
        print(json.dumps(r, indent=2))
        return 0 if r["all_ok"] else 1

    if "error" in r:
        print(f"ERROR: {r['error']}")
        return 1

    def chk(ok):
        return "[PASS]" if ok else "[FAIL]"

    print("Backend parity check on demo1")
    print(f"  reference: {r['ref_path']}")
    print(f"  candidate: {r['cand_path']}")
    print(
        f"  particles: {r['n_active_particles']}  sim_t: 0 → {r['sim_time_s']:.2f}s  samples: {r['n_samples']}"
    )
    print()
    print(
        f"  {chk(r['init_ok'])} initial per-particle L2 = {r['init_l2']:.3f} "
        f"(threshold < {THRESHOLDS['init_l2_max']})"
    )
    print(
        f"  {chk(r['mean_ok'])} mean per-particle L2 over time = {r['mean_per_particle_l2']:.3f} "
        f"(threshold < {THRESHOLDS['mean_per_particle_l2_max']})"
    )
    print(
        f"  {chk(r['center_ok'])} max cube-center diff = {r['max_cube_center_diff']:.3f} "
        f"(threshold < {THRESHOLDS['cube_center_diff_max']})"
    )
    print(
        f"  {chk(r['extent_ok'])} max cube-extent diff = {r['max_cube_extent_diff']:.3f} "
        f"(threshold < {THRESHOLDS['cube_extent_diff_max']})"
    )
    print(
        f"  {chk(r['final_ok'])} final per-particle L2 = {r['final_l2']:.3f} "
        f"(threshold < {THRESHOLDS['final_l2_max']})"
    )
    print()
    if r["all_ok"]:
        print("  [OVERALL PASS] backends agree on demo1 within tolerance")
    else:
        print("  [OVERALL FAIL] substrates have diverged on demo1")
    return 0 if r["all_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
