#!/usr/bin/env python3
"""Run the same demo1 cube-drop scenario through one or more local
Sibernetic-compatible binaries and compare cube-stability metrics.

The "gold standard" is whichever local OpenCL Sibernetic binary you
point at first via --binary. Every subsequent binary's output is
compared against it. Each --binary takes a NAME=PATH form so the table
labels are readable.

Usage:

    # Single backend — just print metrics, no comparison
    python3 scripts/cross_backend_regression.py \\
        --binary "OpenCL=./Release/Sibernetic"

    # Compare two backends — first one is the baseline
    python3 scripts/cross_backend_regression.py \\
        --binary "OpenCL=./Release/Sibernetic" \\
        --binary "Metal-native=src/metal_diff/sib_metal"

    # Override scenario / sim length
    python3 scripts/cross_backend_regression.py --config demo1 --timelimit 1.0 \\
        --binary "OpenCL=./Release/Sibernetic"

The native-Metal substrate (src/metal_diff/sib_metal) takes a different
CLI than the main Sibernetic binary; for that path use
src/metal_diff/dump_metal_trajectory.py to generate a position_buffer.txt
and feed it to scripts/measure_cube_stability.py directly.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

# Tolerance bands for the cube-stability metric. Calibrated against
# OpenCL on demo1 / 1-sec — accept anything in [80%, 200%] extent
# retention. Outside that means the cube either pancaked or grew
# unphysically.
EXTENT_RETENTION_OK = (0.80, 2.00)
# Mean Y must drop by at least 10× — confirms the cube actually fell.
MEAN_Y_FELL_FACTOR = 0.5


def run_local(
    *, binary: Path, config: str, timelimit: float, logstep: int, out_dir: Path
) -> dict:
    """Run a Sibernetic binary locally, return a result-shaped dict."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(binary),
        "-no_g",
        "-f", config,
        f"timelimit={timelimit}",
        f"logstep={logstep}",
        "-l_to",
        f"lpath={out_dir}",
    ]
    start = time.monotonic()
    proc = subprocess.run(
        cmd,
        cwd=str(binary.parent.parent),  # binary is at <root>/Release/Sibernetic
        capture_output=True,
        text=True,
        timeout=1800,
    )
    wall = time.monotonic() - start

    pos_file = out_dir / "position_buffer.txt"
    metrics = None
    if pos_file.exists():
        sys.path.insert(0, str(SCRIPT_DIR))
        from measure_cube_stability import measure  # noqa: E402
        try:
            metrics = measure(pos_file)
        except Exception as e:
            metrics = {"error": repr(e)}

    return {
        "status": "succeeded" if proc.returncode == 0 else "failed",
        "result": {
            "wall_clock_seconds": wall,
            "sibernetic_exit_code": proc.returncode,
            "stability_metrics": metrics,
            "sibernetic_stdout_tail": proc.stdout[-2000:],
            "sibernetic_stderr_tail": proc.stderr[-2000:],
        },
    }


def evaluate(metrics: dict | None) -> tuple[str, str]:
    """Return (verdict, why) for a cube-stability metrics dict."""
    if not metrics or "extent_retention" not in metrics:
        return "FAIL", "no metrics"
    ret = metrics["extent_retention"]
    init_y = metrics["initial_elastic_mean_y"]
    final_y = metrics["final_elastic_mean_y"]
    fell = final_y / init_y if init_y > 0 else 1.0
    notes = []
    if not (EXTENT_RETENTION_OK[0] <= ret <= EXTENT_RETENTION_OK[1]):
        notes.append(f"extent retention {ret:.1%} outside [{EXTENT_RETENTION_OK[0]:.0%}, {EXTENT_RETENTION_OK[1]:.0%}]")
    if fell > MEAN_Y_FELL_FACTOR:
        notes.append(f"cube didn't fall (mean_y {init_y:.1f} → {final_y:.1f})")
    if notes:
        return "FAIL", "; ".join(notes)
    return "PASS", "intact + fell as expected"


def parse_binary_arg(s: str) -> tuple[str, Path]:
    """Parse 'NAME=PATH' or just 'PATH' (uses basename as name)."""
    if "=" in s:
        name, path = s.split("=", 1)
        return name, Path(path)
    p = Path(s)
    return p.name, p


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--config", default="demo1")
    p.add_argument("--timelimit", type=float, default=1.0)
    p.add_argument("--logstep", type=int, default=500)
    p.add_argument(
        "--binary", action="append", default=[],
        help="A `NAME=PATH` (or bare PATH) of a Sibernetic-compatible "
             "binary to run. Pass multiple times to compare. The first "
             "successful run with PASS verdict becomes the baseline.",
    )
    p.add_argument(
        "--out-dir", type=Path, default=Path(f"/tmp/cross_backend_{int(time.time())}"),
        help="Where to put per-backend output buffers.",
    )
    args = p.parse_args(argv)

    if not args.binary:
        p.error("at least one --binary is required")

    runs: list[dict] = []

    for spec in args.binary:
        name, binary = parse_binary_arg(spec)
        if not binary.exists():
            print(f"binary not found: {binary}", file=sys.stderr)
            return 2
        print(f"==> {name}: running {binary} ...", flush=True)
        t0 = time.monotonic()
        resp = run_local(
            binary=binary,
            config=args.config,
            timelimit=args.timelimit,
            logstep=args.logstep,
            out_dir=args.out_dir / name,
        )
        wall = time.monotonic() - t0
        runs.append({"name": name, "wall_clock_outer": wall, "resp": resp})
        inner = resp["result"].get("wall_clock_seconds")
        print(f"    done in {wall:.1f}s outer (inner sim {inner:.1f}s)")

    # Summary
    print()
    print(f"=== Cross-backend regression: {args.config} / timelimit={args.timelimit}s ===")
    print(f"{'backend':28} {'sim wall':>10} {'extent':>10} {'mean_y init→final':>22} {'verdict':>8}")
    print("-" * 80)

    baseline_metrics = None
    baseline_name = None
    for run in runs:
        name = run["name"]
        res = run["resp"].get("result") or {}
        m = res.get("stability_metrics") or {}
        inner = res.get("wall_clock_seconds", 0.0) or 0.0
        ret = m.get("extent_retention")
        i_y = m.get("initial_elastic_mean_y", 0)
        f_y = m.get("final_elastic_mean_y", 0)
        verdict, why = evaluate(m)
        # First PASS becomes the baseline.
        if verdict == "PASS" and baseline_metrics is None:
            baseline_metrics = m
            baseline_name = name
        ret_str = f"{ret:.1%}" if ret is not None else "—"
        y_str = f"{i_y:.1f}→{f_y:.1f}" if m else "—"
        print(f"{name:28} {inner:>9.1f}s {ret_str:>10} {y_str:>22} {verdict:>8}  {why}")

    if baseline_metrics and len(runs) > 1:
        print()
        print(f"=== diff vs baseline ({baseline_name}) ===")
        b_ret = baseline_metrics["extent_retention"]
        b_y = baseline_metrics["final_elastic_mean_y"]
        for run in runs:
            if run["name"] == baseline_name:
                continue
            m = (run["resp"].get("result") or {}).get("stability_metrics") or {}
            if not m or "extent_retention" not in m:
                print(f"  {run['name']}: no metrics to compare")
                continue
            d_ret = m["extent_retention"] - b_ret
            d_y = m["final_elastic_mean_y"] - b_y
            print(f"  {run['name']}: Δretention {d_ret:+.1%}   Δfinal_mean_y {d_y:+.2f}")

    fails = [
        run["name"] for run in runs
        if evaluate((run["resp"].get("result") or {}).get("stability_metrics"))[0] == "FAIL"
    ]
    if fails:
        print(f"\n{len(fails)} backend(s) FAILED: {', '.join(fails)}")
        return 1
    print(f"\nAll {len(runs)} backend(s) PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
