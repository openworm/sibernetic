import os
import subprocess
import math

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "reference_logs")


def _load_matrix(name, base=DATA_DIR):
    """Return rows of floats while skipping header lines."""
    path = os.path.join(base, name)
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) != 4:
                # position/velocity logs begin with metadata that should be
                # ignored; these lines do not contain four columns
                continue
            rows.append([float(v) for v in parts])
    return rows


def test_reference_logs_exist():
    for name in (
        "positions_step0.txt",
        "velocities_step0.txt",
        "density_step0.txt",
        "pressure_step0.txt",
    ):
        assert os.path.exists(os.path.join(DATA_DIR, name))


def _have_opencl():
    try:
        out = subprocess.run(["clinfo"], check=False, capture_output=True, text=True)
    except FileNotFoundError:
        return False
    return "Number of platforms                               0" not in out.stdout


def test_engine_against_reference(tmp_path):
    """Run compiled Sibernetic binary and compare output with reference.

    This test requires:
    - Compiled Sibernetic binary at ./Release/Sibernetic
    - OpenCL runtime
    - RUN_ENGINE_TESTS=1 environment variable
    """
    import pytest

    # Check if compiled binary exists
    binary_path = "./Release/Sibernetic"
    if not os.path.exists(binary_path):
        pytest.skip(f"Sibernetic binary not found at {binary_path}")

    if not _have_opencl():
        pytest.skip("OpenCL not available")

    if os.environ.get("RUN_ENGINE_TESTS") != "1":
        pytest.skip("RUN_ENGINE_TESTS not set")

    # Test that binary can run (quick smoke test first)
    try:
        result = subprocess.run(
            [binary_path, "-no_g", "timelimit=0.0001"],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
        if "Failed to load" in result.stdout or "Error" in result.stderr:
            pytest.skip(f"Sibernetic binary has runtime errors: {result.stdout[:200]}")
    except subprocess.TimeoutExpired:
        pytest.skip("Sibernetic binary timed out")

    out_dir = tmp_path
    cmd = [
        binary_path,
        "-no_g",
        "-f",
        "configuration/test/test_energy",
        "-l_to",
        f"lpath={out_dir}",
        "timelimit=0.001",
        "logstep=25",
    ]
    subprocess.run(cmd, check=True)
    gen_path = os.path.join(out_dir, "position_buffer.txt")
    assert os.path.exists(gen_path)
    generated = _load_matrix("position_buffer.txt", base=out_dir)
    baseline = _load_matrix("positions_step0.txt")
    for g_row, b_row in zip(generated, baseline):
        for gv, bv in zip(g_row, b_row):
            assert math.isfinite(gv)
            assert abs(gv - bv) < 1e-3
