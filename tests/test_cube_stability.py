"""Cross-backend regression test: the elastic cube must keep its shape
after hitting the floor.

Runs the demo1 cube-drop scenario for ~3 seconds of simulated time on
each available backend and asserts that the elastic body retains at
least 50% of its initial vertical extent. The "pancake" failure mode
documented in the README (Taichi flattens to mean Y ~0.24 vs PyTorch's
1.45) registers as < 10% retention here.

This is the objective test backing what until now was a visual judgement
made from rendered MP4s.

Marked `slow` + `integration`: needs the compiled binary, takes minutes
to run on per-backend wall-clock. Skip in per-PR CI; opt in via
`pytest -m slow tests/test_cube_stability.py`.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SIBERNETIC_BIN = REPO_ROOT / "Release" / "Sibernetic"
SCRIPTS_DIR = REPO_ROOT / "scripts"

# Make the standalone metric script importable as a module so the test
# logic stays in lock-step with the human-facing CLI tool.
sys.path.insert(0, str(SCRIPTS_DIR))
from measure_cube_stability import measure  # noqa: E402

# Sim parameters: demo1 starts the cube at Y=[39.4, 49.4] and the floor
# is at Y=0. Falling under gravity, first impact happens around t~3-4s.
# We need post-impact dynamics to distinguish "cube lands intact and
# settles" from "cube collapses on contact" — running long enough that
# the elastic extent has had time to either recover or stay flat. 5s
# total sim matches the local validation buffer that surfaced the
# pancake at 8% retention.
DEMO1_TIMELIMIT_SECONDS = 5.0
DEMO1_LOGSTEP = 100  # ~50 logged frames per run; enough for first/last
DEMO1_RUN_TIMEOUT_SECONDS = 1800  # 30 min upper bound for slowest backend
EXTENT_RETENTION_THRESHOLD = 0.5  # 50% of initial elastic-Y extent


def _backend_runnable(backend: str) -> tuple[bool, str]:
    """Return (runnable, reason) for the given backend on this host."""
    arch = os.uname().machine
    sysname = os.uname().sysname

    if backend == "opencl":
        # Apple Silicon Sibernetic is built with -DOW_NO_OPENCL.
        if sysname == "Darwin" and arch == "arm64":
            return False, "OpenCL disabled at compile time on Apple Silicon"
        if not shutil.which("clinfo"):
            return True, ""  # clinfo isn't strictly required; let it run
        return True, ""

    if backend == "torch":
        try:
            import torch  # noqa: F401
            return True, ""
        except ImportError:
            return False, "torch not importable from system Python"

    if backend in ("taichi", "taichi-cpu", "taichi-cuda"):
        try:
            import taichi  # noqa: F401
        except ImportError:
            return False, "taichi not importable from system Python"
        if backend == "taichi" and not (sysname == "Darwin" and arch == "arm64"):
            return False, "backend=taichi defaults to Metal; Apple Silicon only"
        if backend == "taichi-cuda":
            # No reliable runtime check that doesn't require importing
            # taichi and probing arch=cuda; skip on non-Linux for now.
            if sysname != "Linux":
                return False, "Taichi-CUDA only attempted on Linux hosts"
        return True, ""

    return False, f"unknown backend {backend!r}"


def _run_demo1(backend: str, out_dir: Path) -> Path:
    """Run demo1 with the given backend, return path to position_buffer.txt."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(SIBERNETIC_BIN),
        "-no_g",
        "-f", "demo1",
        f"backend={backend}",
        f"timelimit={DEMO1_TIMELIMIT_SECONDS}",
        f"logstep={DEMO1_LOGSTEP}",
        "-l_to",
        f"lpath={out_dir}",
    ]
    subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        check=True,
        timeout=DEMO1_RUN_TIMEOUT_SECONDS,
    )
    pos_file = out_dir / "position_buffer.txt"
    assert pos_file.exists(), f"binary completed but no buffer at {pos_file}"
    return pos_file


@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.parametrize("backend", ["opencl", "torch", "taichi-cpu", "taichi"])
def test_demo1_cube_retains_structure_after_floor_collision(backend, tmp_path):
    """Cube must keep ≥50% of its initial vertical extent after impact.

    Catches the 'pancake' failure mode where the elastic body collapses
    flat instead of maintaining its 3D shape via internal pressure.
    Reference numbers from the README (3s sim, demo1):
      - PyTorch (working) : elastic mean Y ≈ 1.45 (kept structure)
      - Taichi  (broken)  : elastic mean Y ≈ 0.24 (pancake)
    """
    if not SIBERNETIC_BIN.exists():
        pytest.skip(f"Sibernetic binary not built at {SIBERNETIC_BIN}")
    runnable, reason = _backend_runnable(backend)
    if not runnable:
        pytest.skip(f"backend={backend}: {reason}")

    pos_file = _run_demo1(backend, tmp_path / "buffers")
    m = measure(pos_file)

    # Print metrics so a failure has full context in the test log.
    print(
        f"\nbackend={backend}: "
        f"extent {m['initial_elastic_extent_y']:.2f} → {m['final_elastic_extent_y']:.2f} "
        f"(retention {m['extent_retention']:.1%}), "
        f"mean_Y {m['initial_elastic_mean_y']:.2f} → {m['final_elastic_mean_y']:.2f}, "
        f"min_Y {m['initial_elastic_min_y']:.2f} → {m['final_elastic_min_y']:.2f}, "
        f"frames={m['n_frames']}"
    )

    assert m["extent_retention"] >= EXTENT_RETENTION_THRESHOLD, (
        f"backend={backend}: cube collapsed after floor collision. "
        f"Initial elastic-Y extent {m['initial_elastic_extent_y']:.2f}, "
        f"final {m['final_elastic_extent_y']:.2f}, "
        f"retention {m['extent_retention']:.1%} "
        f"(threshold {EXTENT_RETENTION_THRESHOLD:.0%}). "
        f"Final mean_Y={m['final_elastic_mean_y']:.2f}, "
        f"min_Y={m['final_elastic_min_y']:.2f}. "
        f"This is the 'pancake' failure mode described in the README — "
        f"elastic forces aren't strong enough to resist gravity + SPH "
        f"pressure on impact, so the cube flattens against the floor."
    )
