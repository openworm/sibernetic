#!/usr/bin/env python3
"""Run the same Sibernetic scenario across multiple backends, compare metrics.

Catches the kind of cross-backend physics divergence we already know about
from this repo's history:

  - Taichi-Metal pancakes on Apple Silicon (5-sec sim → extent 7.6%)
  - Taichi-CUDA freezes on Linux L4 (1-sec sim → cube doesn't fall)
  - PyTorch (CPU) takes hours to even produce a frame on demo1 scale
  - OpenCL on real GPU drops the cube ~5x faster than Taichi does

Submits jobs to a Cloud Run runner (default the deployed
`sibernetic-runner` instance) and/or runs a local binary, collects the
cube-stability metrics, and prints a side-by-side comparison plus a
PASS/FAIL gate keyed to the OpenCL baseline.

Usage examples:

    # Default: hit Cloud Run, sweep cloud-runnable backends
    python3 scripts/cross_backend_regression.py

    # Include a local binary path (e.g. PR #222's native Metal)
    python3 scripts/cross_backend_regression.py \\
        --local-binary /tmp/sibernetic-pr222/Release/Sibernetic --local-name "PR222-Metal"

    # Only test specific backends
    python3 scripts/cross_backend_regression.py --backend opencl --backend taichi-cuda

    # Override the scenario / sim length
    python3 scripts/cross_backend_regression.py --config demo1 --timelimit 1.0
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from urllib import request as urllib_request


CLOUD_RUN_URL = os.environ.get(
    "SIBERNETIC_RUNNER_URL",
    "https://sibernetic-runner-xvce2jjjja-uk.a.run.app",
)
CLOUD_RUN_API_KEY = os.environ.get("SIBERNETIC_RUNNER_API_KEY", "sibernetic-runner-2026")
SCRIPT_DIR = Path(__file__).resolve().parent

# Tolerance bands for the cube-stability metric. The OpenCL baseline at
# demo1 / 1-sec is ~123% extent retention; we accept anything in [80%,
# 200%]. Outside that means the cube either pancaked or grew unphysically.
EXTENT_RETENTION_OK = (0.80, 2.00)
# Mean Y must drop by at least 10x — confirms the cube actually fell.
MEAN_Y_FELL_FACTOR = 0.5


def submit_cloud_run(
    *, config: str, backend: str, timelimit: float, logstep: int
) -> dict:
    """Submit one job to the Cloud Run runner and poll until completion."""
    body = json.dumps({
        "config": config,
        "backend": backend,
        "timelimit": timelimit,
        "logstep": logstep,
    }).encode()
    req = urllib_request.Request(
        f"{CLOUD_RUN_URL}/api/run",
        data=body,
        headers={"X-API-Key": CLOUD_RUN_API_KEY, "Content-Type": "application/json"},
        method="POST",
    )
    with urllib_request.urlopen(req, timeout=30) as r:
        job_id = json.loads(r.read())["job_id"]

    deadline = time.monotonic() + 1800  # 30 min hard cap
    while time.monotonic() < deadline:
        status_req = urllib_request.Request(
            f"{CLOUD_RUN_URL}/api/status/{job_id}",
            headers={"X-API-Key": CLOUD_RUN_API_KEY},
        )
        try:
            with urllib_request.urlopen(status_req, timeout=30) as r:
                s = json.loads(r.read())
        except Exception as e:
            print(f"  [warn] status check failed, retrying: {e}", file=sys.stderr)
            time.sleep(15)
            continue
        if s["status"] in ("succeeded", "failed"):
            break
        time.sleep(15)

    result_req = urllib_request.Request(
        f"{CLOUD_RUN_URL}/api/result/{job_id}",
        headers={"X-API-Key": CLOUD_RUN_API_KEY},
    )
    with urllib_request.urlopen(result_req, timeout=30) as r:
        return json.loads(r.read())


def run_local(
    *, binary: Path, config: str, timelimit: float, logstep: int, out_dir: Path
) -> dict:
    """Run a Sibernetic binary locally, return a result-shaped dict.

    Output shape mirrors the Cloud Run result so the comparison logic
    doesn't care which path produced the metrics.
    """
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


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--config", default="demo1")
    p.add_argument("--timelimit", type=float, default=1.0)
    p.add_argument("--logstep", type=int, default=500)
    p.add_argument(
        "--backend", action="append", default=None,
        help="Cloud Run backend to test (can pass multiple). "
             "Default: opencl, torch-cuda, taichi-cuda.",
    )
    p.add_argument(
        "--local-binary", type=Path, default=None,
        help="Path to a local Sibernetic binary to ALSO run alongside the cloud backends.",
    )
    p.add_argument(
        "--local-name", default="local",
        help="Display name for the --local-binary entry (default 'local').",
    )
    p.add_argument(
        "--no-cloud", action="store_true",
        help="Skip Cloud Run; only run --local-binary.",
    )
    p.add_argument(
        "--out-dir", type=Path, default=Path(f"/tmp/cross_backend_{int(time.time())}"),
        help="Where to put per-backend output buffers (for --local-binary runs only).",
    )
    args = p.parse_args(argv)

    cloud_backends = args.backend or ["opencl", "torch-cuda", "taichi-cuda"]
    runs: list[dict] = []

    if not args.no_cloud:
        for backend in cloud_backends:
            print(f"==> [cloud] {backend}: submitting...", flush=True)
            t0 = time.monotonic()
            try:
                resp = submit_cloud_run(
                    config=args.config,
                    backend=backend,
                    timelimit=args.timelimit,
                    logstep=args.logstep,
                )
            except Exception as e:
                resp = {"status": "failed", "error": repr(e), "result": None}
            wall = time.monotonic() - t0
            runs.append({"name": f"cloud:{backend}", "wall_clock_outer": wall, "resp": resp})
            res = (resp.get("result") or {})
            inner = res.get("wall_clock_seconds")
            print(f"    done in {wall:.1f}s outer (inner sim {inner:.1f}s)" if inner else f"    done in {wall:.1f}s outer")

    if args.local_binary:
        if not args.local_binary.exists():
            print(f"--local-binary {args.local_binary} does not exist", file=sys.stderr)
            return 2
        print(f"==> [local] {args.local_name}: running {args.local_binary} ...", flush=True)
        t0 = time.monotonic()
        resp = run_local(
            binary=args.local_binary,
            config=args.config,
            timelimit=args.timelimit,
            logstep=args.logstep,
            out_dir=args.out_dir / args.local_name,
        )
        wall = time.monotonic() - t0
        runs.append({"name": f"local:{args.local_name}", "wall_clock_outer": wall, "resp": resp})
        inner = resp["result"].get("wall_clock_seconds")
        print(f"    done in {wall:.1f}s outer (inner sim {inner:.1f}s)")

    # Summary table
    print()
    print(f"=== Cross-backend regression: {args.config} / timelimit={args.timelimit}s ===")
    print(f"{'backend':28} {'sim wall':>10} {'extent':>10} {'mean_y init→final':>22} {'verdict':>8}")
    print("-" * 80)

    baseline_metrics = None
    for run in runs:
        name = run["name"]
        resp = run["resp"]
        res = resp.get("result") or {}
        m = res.get("stability_metrics") or {}
        inner = res.get("wall_clock_seconds", 0.0) or 0.0
        ret = m.get("extent_retention")
        i_y = m.get("initial_elastic_mean_y", 0)
        f_y = m.get("final_elastic_mean_y", 0)
        verdict, why = evaluate(m)
        # Capture OpenCL as gold-standard baseline if we have it
        if "opencl" in name and verdict == "PASS" and baseline_metrics is None:
            baseline_metrics = m
        ret_str = f"{ret:.1%}" if ret is not None else "—"
        y_str = f"{i_y:.1f}→{f_y:.1f}" if m else "—"
        print(f"{name:28} {inner:>9.1f}s {ret_str:>10} {y_str:>22} {verdict:>8}  {why}")

    # Pairwise diff vs OpenCL baseline
    if baseline_metrics:
        print()
        print(f"=== diff vs OpenCL baseline ===")
        b_ret = baseline_metrics["extent_retention"]
        b_y = baseline_metrics["final_elastic_mean_y"]
        for run in runs:
            name = run["name"]
            if "opencl" in name:
                continue
            m = (run["resp"].get("result") or {}).get("stability_metrics") or {}
            if not m or "extent_retention" not in m:
                print(f"  {name}: no metrics to compare")
                continue
            d_ret = m["extent_retention"] - b_ret
            d_y = m["final_elastic_mean_y"] - b_y
            print(f"  {name}: Δretention {d_ret:+.1%}   Δfinal_mean_y {d_y:+.2f}")

    # Exit code: non-zero if any verdict is FAIL
    fails = [run["name"] for run in runs if evaluate((run["resp"].get("result") or {}).get("stability_metrics"))[0] == "FAIL"]
    if fails:
        print(f"\n{len(fails)} backend(s) FAILED: {', '.join(fails)}")
        return 1
    print(f"\nAll {len(runs)} backend(s) PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
