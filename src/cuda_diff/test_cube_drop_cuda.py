"""Phase 3 acceptance test — cube-drop with the CUDA xpbd_step driver.

Generates a 3x3x3 elastic cube, drops it onto a floor plane under gravity,
runs N XPBD steps, dumps a Sibernetic-format position_buffer.txt, and
invokes scripts/measure_cube_stability.py to check that the cube didn't
pancake (extent retention >= 0.5).

Local-only. No Metal needed — the cube doesn't pancake or it does, and
that single number is the Phase 3 acceptance gate per src/cuda/README.md.
"""
import os
import platform
import struct
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
BINARY_NAME = "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda"
BINARY = os.path.join(HERE, BINARY_NAME)
TMP = tempfile.gettempdir()
MEASURE = os.path.join(REPO, "scripts", "measure_cube_stability.py")


def make_cube(side=3, spacing=1.0, y_offset=4.0):
    """Return (positions [N,3], bonds list of (i,j,rest_len)).

    Bonds include nearest neighbours + face diagonals + body diagonals
    (every pair within sqrt(3)*spacing). Axis-aligned bonds alone give a
    hinged frame that collapses under gravity; adding the diagonal struts
    gives full 3D rigidity.
    """
    coords = []
    for k in range(side):
        for j in range(side):
            for i in range(side):
                x = (i - (side - 1) / 2.0) * spacing
                z = (k - (side - 1) / 2.0) * spacing
                y = j * spacing + y_offset
                coords.append((x, y, z))
    positions = np.array(coords, dtype=np.float32)

    bonds = []
    cutoff = spacing * np.sqrt(3) + 1e-3
    n = len(positions)
    for a in range(n):
        for b in range(a + 1, n):
            d = float(np.linalg.norm(positions[a] - positions[b]))
            if d <= cutoff:
                bonds.append((a, b, d))
    return positions, bonds


def write_bonds_bin(path, bonds):
    with open(path, "wb") as f:
        for i, j, L in bonds:
            f.write(struct.pack("ii f f", int(i), int(j), float(L), 0.0))


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1

    positions, bonds = make_cube(side=3, spacing=1.0, y_offset=4.0)
    n_active = len(positions)
    n_bonds = len(bonds)
    print(f"Cube: {n_active} particles, {n_bonds} bonds, "
          f"y in [{positions[:,1].min()}, {positions[:,1].max()}]")

    velocities = np.zeros_like(positions)
    static_p = np.zeros((0, 3), dtype=np.float32)  # no boundary particles

    p_pos = os.path.join(TMP, "sib_cuda_xpbd_pos.bin")
    p_vel = os.path.join(TMP, "sib_cuda_xpbd_vel.bin")
    p_static = os.path.join(TMP, "sib_cuda_xpbd_static.bin")
    p_bonds = os.path.join(TMP, "sib_cuda_xpbd_bonds.bin")
    p_pos_out = os.path.join(TMP, "sib_cuda_xpbd_pos_out.bin")
    p_vel_out = os.path.join(TMP, "sib_cuda_xpbd_vel_out.bin")
    p_traj = os.path.join(TMP, "sib_cuda_xpbd_traj.txt")

    positions.astype(np.float32).tofile(p_pos)
    velocities.astype(np.float32).tofile(p_vel)
    static_p.astype(np.float32).tofile(p_static)
    write_bonds_bin(p_bonds, bonds)

    # Physics knobs. Stiff bonds (alpha small) hold the cube; gravity is
    # gentle so we don't need very small dt for stability.
    args = [BINARY, "xpbd_step",
            str(n_active),             # n_active
            "0",                       # n_static
            "1.0",                     # h (smoothing radius — unused, density skipped)
            "1.0",                     # mass
            "10.0",                    # rho_rest (unused; skip_density=1)
            "0.05",                    # dt
            "-1.0",                    # gravity_y
            "0.0",                     # floor_y
            "1e-6",                    # alpha_density (unused)
            "10",                      # n_iters (inner XPBD)
            "200",                     # n_steps
            p_pos, p_vel, p_static,
            str(n_bonds), p_bonds,
            "1e-6",                    # alpha_dist (compliance, stiff)
            p_pos_out, p_vel_out,
            p_traj,                    # traj_out
            "5",                       # traj_every steps
            "1.0",                     # sim_scale
            "0.0",                     # vel_clamp (off)
            "0.0",                     # restitution
            "1",                       # skip_density = 1
            ]
    print("Running xpbd_step ...")
    subprocess.run(args, check=True)

    print(f"\nTrajectory written: {p_traj}")
    print(f"Invoking measure_cube_stability.py ...\n")
    rc = subprocess.run(
        [sys.executable, MEASURE, p_traj, "--threshold", "0.5"]
    ).returncode
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
