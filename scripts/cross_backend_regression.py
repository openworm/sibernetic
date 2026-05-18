#!/usr/bin/env python3
"""Run the same demo1 cube-drop scenario through one or more local
Sibernetic-compatible substrates and compare cube-stability metrics.

The "gold standard" is whichever local OpenCL Sibernetic binary you
point at first via --binary. Every subsequent binary's output is
compared against it. Each --binary takes a NAME=PATH form so the table
labels are readable.

Substrates
----------

Three substrate kinds are dispatched by `--substrate` (paired
index-aligned with `--binary`):

* `opencl` (default)
    The main Sibernetic binary CLI:
        sibernetic -no_g -f <config> timelimit=<T> logstep=<L> -l_to lpath=<dir>
    Produces `position_buffer.txt` in `lpath`.

* `cuda`
    Native CUDA differentiable substrate (`src/cuda_diff/sib_cuda`).
    Its CLI is op-positional (`xpbd_full_fwd ...`), not Sibernetic CLI,
    so we shell out to `src/cuda_diff/dump_cuda_full_trajectory.py`
    which locates `sib_cuda` HERE-relative and emits a Sibernetic-format
    `position_buffer.txt`. The `--binary` value is used as a display
    label here; the actual binary is found by the helper.

* `metal-native`
    Native Metal differentiable substrate (`src/metal_diff/sib_metal`).
    Same shape as cuda: we invoke `src/metal_diff/dump_metal_trajectory.py`
    and let it find `sib_metal` HERE-relative. `--binary` is a label.

Downstream `measure(pos_file)` is unchanged across all three substrates,
so the regression table and PASS/FAIL bands are identical.

Usage
-----

    # Single OpenCL backend - just print metrics, no comparison.
    python scripts/cross_backend_regression.py \\
        --binary "OpenCL=./Release/Sibernetic"

    # Compare OpenCL against CUDA on demo1.
    python scripts/cross_backend_regression.py \\
        --binary "OpenCL=./Release/Sibernetic" \\
        --binary "CUDA=src/cuda_diff/sib_cuda.exe" --substrate opencl --substrate cuda

    # Three-way OpenCL vs Metal vs CUDA on demo1 / 1s.
    python scripts/cross_backend_regression.py \\
        --binary "OpenCL=./Release/Sibernetic" --substrate opencl \\
        --binary "Metal=src/metal_diff/sib_metal" --substrate metal-native \\
        --binary "CUDA=src/cuda_diff/sib_cuda.exe" --substrate cuda

If --substrate is omitted entirely, every binary defaults to `opencl`
(matches the legacy behavior of this script). If --substrate is given
at all, the number of --substrate values must equal the number of
--binary values.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

SUBSTRATE_CHOICES = ("opencl", "cuda", "metal-native")

# Tolerance bands for the cube-stability metric. Calibrated against
# OpenCL on demo1 / 1-sec - accept anything in [80%, 200%] extent
# retention. Outside that means the cube either pancaked or grew
# unphysically.
EXTENT_RETENTION_OK = (0.80, 2.00)
# Mean Y must drop by at least 10x - confirms the cube actually fell.
MEAN_Y_FELL_FACTOR = 0.5


def run_opencl(
    *, binary: Path, config: str, timelimit: float, logstep: int, out_dir: Path
) -> dict:
    """Run the main Sibernetic OpenCL binary, return a result-shaped dict."""
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
    metrics = _try_measure(pos_file)

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


def run_cuda(
    *, config: str, timelimit: float, logstep: int, out_dir: Path,
    dt: float,
) -> dict:
    """Run the CUDA differentiable substrate via dump_cuda_full_trajectory.py.

    The helper script locates `sib_cuda` HERE-relative; we pass it the
    same `config`, derived `K = round(timelimit/dt)`, and `logstep` the
    OpenCL path uses, then point `--out` at a per-run file inside
    `out_dir` named `position_buffer.txt` so the downstream `measure()`
    flow stays identical.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    pos_file = out_dir / "position_buffer.txt"
    K = max(1, int(round(timelimit / dt)))
    dump_script = REPO_ROOT / "src" / "cuda_diff" / "dump_cuda_full_trajectory.py"
    cmd = [
        sys.executable, str(dump_script),
        "--scenario", config,
        "--K", str(K),
        "--logstep", str(logstep),
        "--dt", str(dt),
        "--out", str(pos_file),
    ]
    start = time.monotonic()
    proc = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=1800,
    )
    wall = time.monotonic() - start

    metrics = _try_measure(pos_file)

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


def run_metal_native(
    *, config: str, timelimit: float, logstep: int, out_dir: Path,
    dt: float,
) -> dict:
    """Run the native Metal substrate via dump_metal_trajectory.py.

    Mirrors `run_cuda` but uses the Metal helper. Metal's helper exposes
    `--steps` (total XPBD steps, same role as CUDA's `--K`) and
    `--chunk` (snapshot stride, same role as `--logstep`).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    pos_file = out_dir / "position_buffer.txt"
    steps = max(1, int(round(timelimit / dt)))
    dump_script = REPO_ROOT / "src" / "metal_diff" / "dump_metal_trajectory.py"
    cmd = [
        sys.executable, str(dump_script),
        "--scenario", config,
        "--steps", str(steps),
        "--chunk", str(logstep),
        "--dt", str(dt),
        "--out", str(pos_file),
    ]
    start = time.monotonic()
    proc = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=1800,
    )
    wall = time.monotonic() - start

    metrics = _try_measure(pos_file)

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


def _try_measure(pos_file: Path) -> dict | None:
    """Run scripts/measure_cube_stability.measure() on `pos_file` if it exists."""
    if not pos_file.exists():
        return None
    sys.path.insert(0, str(SCRIPT_DIR))
    try:
        from measure_cube_stability import measure  # noqa: E402
        return measure(pos_file)
    except Exception as e:
        return {"error": repr(e)}


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
        notes.append(f"cube didn't fall (mean_y {init_y:.1f} -> {final_y:.1f})")
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
        "--dt", type=float, default=2e-5,
        help="Per-step timestep used to derive K = round(timelimit/dt) for "
             "the CUDA and metal-native substrates (default 2e-5, matches "
             "src/cuda_diff/dump_cuda_full_trajectory.py). Ignored for the "
             "opencl substrate, which reads dt from the config.",
    )
    p.add_argument(
        "--binary", action="append", default=[],
        help="A `NAME=PATH` (or bare PATH) of a Sibernetic-compatible "
             "binary to run. Pass multiple times to compare. The first "
             "successful run with PASS verdict becomes the baseline. For "
             "substrates other than `opencl`, the value is used as a "
             "display label only; the actual binary is located by the "
             "matching `src/*/dump_*_trajectory.py` helper.",
    )
    p.add_argument(
        "--substrate", action="append", default=[],
        choices=SUBSTRATE_CHOICES,
        help="Substrate kind for the paired `--binary` (index-aligned). "
             "One of: " + ", ".join(SUBSTRATE_CHOICES) + ". If omitted "
             "entirely, every --binary defaults to `opencl` (legacy "
             "behavior). If given at all, must be supplied once per "
             "--binary.",
    )
    p.add_argument(
        "--out-dir", type=Path,
        default=Path(os.path.join(
            os.environ.get("TMP", os.environ.get("TMPDIR", "/tmp")),
            f"cross_backend_{int(time.time())}",
        )),
        help="Where to put per-backend output buffers.",
    )
    args = p.parse_args(argv)

    if not args.binary:
        p.error("at least one --binary is required")

    # Pair --binary with --substrate. If --substrate is unset, default
    # everything to opencl. If set partially, error out.
    if not args.substrate:
        substrates = ["opencl"] * len(args.binary)
    elif len(args.substrate) != len(args.binary):
        p.error(
            f"--substrate given {len(args.substrate)} times but "
            f"--binary given {len(args.binary)} times; counts must match "
            f"(or omit --substrate entirely to default all to opencl)"
        )
    else:
        substrates = list(args.substrate)

    runs: list[dict] = []

    for spec, substrate in zip(args.binary, substrates):
        name, binary = parse_binary_arg(spec)
        if substrate == "opencl":
            if not binary.exists():
                print(f"binary not found: {binary}", file=sys.stderr)
                return 2
        else:
            # For non-opencl substrates the --binary value is a label /
            # convenience pointer; the helper script finds its own
            # binary HERE-relative. We still warn (not error) if the
            # path doesn't exist so the user knows to build it.
            if not binary.exists():
                print(
                    f"note: {name} ({substrate}) --binary path {binary} "
                    f"does not exist; the dump_*_trajectory.py helper "
                    f"will locate the real binary HERE-relative.",
                    file=sys.stderr,
                )
        print(f"==> {name} [{substrate}]: running ...", flush=True)
        t0 = time.monotonic()
        if substrate == "opencl":
            resp = run_opencl(
                binary=binary,
                config=args.config,
                timelimit=args.timelimit,
                logstep=args.logstep,
                out_dir=args.out_dir / name,
            )
        elif substrate == "cuda":
            resp = run_cuda(
                config=args.config,
                timelimit=args.timelimit,
                logstep=args.logstep,
                out_dir=args.out_dir / name,
                dt=args.dt,
            )
        elif substrate == "metal-native":
            resp = run_metal_native(
                config=args.config,
                timelimit=args.timelimit,
                logstep=args.logstep,
                out_dir=args.out_dir / name,
                dt=args.dt,
            )
        else:
            # Argparse `choices=` should make this unreachable.
            p.error(f"unknown substrate: {substrate}")
        wall = time.monotonic() - t0
        runs.append({
            "name": name, "substrate": substrate,
            "wall_clock_outer": wall, "resp": resp,
        })
        inner = resp["result"].get("wall_clock_seconds")
        print(f"    done in {wall:.1f}s outer (inner sim {inner:.1f}s)")

    # Summary
    print()
    print(f"=== Cross-backend regression: {args.config} / timelimit={args.timelimit}s ===")
    print(f"{'backend':28} {'sim wall':>10} {'extent':>10} {'mean_y init->final':>22} {'verdict':>8}")
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
        ret_str = f"{ret:.1%}" if ret is not None else "-"
        y_str = f"{i_y:.1f}->{f_y:.1f}" if m else "-"
        label = f"{name} [{run['substrate']}]"
        print(f"{label:28} {inner:>9.1f}s {ret_str:>10} {y_str:>22} {verdict:>8}  {why}")

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
            print(f"  {run['name']}: dretention {d_ret:+.1%}   dfinal_mean_y {d_y:+.2f}")

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
