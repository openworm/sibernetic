import os
import subprocess
import math

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "reference_logs")


def _load_matrix(name, base=DATA_DIR):
    path = os.path.join(base, name)
    with open(path) as f:
        return [list(map(float, line.split())) for line in f if line.strip()]


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
    if not _have_opencl() or os.environ.get("RUN_ENGINE_TESTS") != "1":
        print("Skipping engine comparison test due to missing OpenCL")
        return
    out_dir = tmp_path
    cmd = [
        "./Release/Sibernetic",
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
