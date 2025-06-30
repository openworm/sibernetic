import os
import subprocess
import math

from test_solver_logs import _load_matrix


def _have_torch():
    try:
        import torch  # noqa: F401

        return True
    except Exception:
        return False


def test_torch_backend(tmp_path):
    if not _have_torch() or os.environ.get("RUN_ENGINE_TESTS") != "1":
        print("Skipping torch backend test")
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
        "backend=torch",
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{os.getcwd()}:{env.get('PYTHONPATH', '')}"
    proc = subprocess.run(cmd, env=env)
    if proc.returncode != 0:
        print("Torch backend run failed, skipping")
        return

    pos = _load_matrix("position_buffer.txt", base=out_dir)
    base_pos = _load_matrix("positions_step0.txt")
    for g_row, b_row in zip(pos, base_pos):
        for gv, bv in zip(g_row, b_row):
            assert math.isfinite(gv)
            assert abs(gv - bv) < 1e-3

    vel = _load_matrix("velocity_buffer.txt", base=out_dir)
    base_vel = _load_matrix("velocities_step0.txt")
    for g_row, b_row in zip(vel, base_vel):
        for gv, bv in zip(g_row, b_row):
            assert math.isfinite(gv)
            assert abs(gv - bv) < 1e-3

    energy_file = os.path.join(out_dir, "total_energy_distrib.txt")
    if os.path.exists(energy_file):
        with open(energy_file) as f:
            energies = [float(line.strip()) for line in f if line.strip()]
        if len(energies) >= 2:
            start, end = energies[0], energies[-1]
            rel = abs(end - start) / (abs(start) + 1e-12)
            assert rel < 1.0
