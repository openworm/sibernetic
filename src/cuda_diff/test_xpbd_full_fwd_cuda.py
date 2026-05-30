"""Structural smoke test for the CUDA xpbd_full_fwd driver.

Not a numerical correctness check (that comes via sub-task F's FD test).
Goal: prove the orchestrator runs end-to-end without crashes/NaNs and the
tape has the correct shape across the four optional-feature combinations.

Scenarios (4 active particles in a 2x2 grid above 8 boundary particles):
  (a) Required-args only (no pair, no springs, no floor).
  (b) visc_pair_coef > 0 (use_pair → +n per step on the tape).
  (c) spring_K > 0 with a single bond (no tape extra).
  (d) floor_y set with restitution > 0 (use_floor → +n per step on tape).
"""
import os
import platform
import struct
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY_NAME = "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda"
BINARY = os.path.join(HERE, BINARY_NAME)
TMP = tempfile.gettempdir()


# ── Scenario fixture ──────────────────────────────────────────────────
N_ACTIVE = 4
N_STATIC = 8
K = 3
H = 0.5
MASS = 1.0
RHO_REST = 4.0
DT = 0.01
G_Y = -9.81
ALPHA_DENS = 1e3


def make_scenario():
    # 4 active particles in a 2x2 grid at y=1.5.
    pos0 = np.array([
        [-0.25,  1.5, -0.25],
        [ 0.25,  1.5, -0.25],
        [-0.25,  1.5,  0.25],
        [ 0.25,  1.5,  0.25],
    ], dtype=np.float32)
    vel0 = np.zeros((N_ACTIVE, 3), dtype=np.float32)

    # 8 boundary particles forming the corners of a small box at y=0.
    pos_static = np.array([
        [-0.5, 0.0, -0.5], [0.5, 0.0, -0.5],
        [-0.5, 0.0,  0.5], [0.5, 0.0,  0.5],
        [-0.5, 0.5, -0.5], [0.5, 0.5, -0.5],
        [-0.5, 0.5,  0.5], [0.5, 0.5,  0.5],
    ], dtype=np.float32)
    assert pos_static.shape[0] == N_STATIC
    return pos0, vel0, pos_static


def write_bonds_bin(path, bonds):
    with open(path, "wb") as f:
        for i, j, L in bonds:
            f.write(struct.pack("ii f f", int(i), int(j), float(L), 0.0))


def run_full_fwd(extra_args, label):
    """Run xpbd_full_fwd with the given optional-arg list appended."""
    pos0, vel0, pos_static = make_scenario()
    p_pos = os.path.join(TMP, f"cuda_xfull_{label}_pos.bin")
    p_vel = os.path.join(TMP, f"cuda_xfull_{label}_vel.bin")
    p_static = os.path.join(TMP, f"cuda_xfull_{label}_static.bin")
    p_state = os.path.join(TMP, f"cuda_xfull_{label}_state.bin")
    pos0.tofile(p_pos)
    vel0.tofile(p_vel)
    pos_static.tofile(p_static)

    base = [BINARY, "xpbd_full_fwd",
            str(N_ACTIVE), str(N_STATIC), str(K),
            str(H), str(MASS), str(RHO_REST), str(DT),
            str(G_Y), str(ALPHA_DENS),
            p_pos, p_vel, p_static, p_state]
    cmd = base + [str(a) for a in extra_args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[{label}] STDERR:\n{r.stderr}", file=sys.stderr)
        print(f"[{label}] STDOUT:\n{r.stdout}", file=sys.stderr)
        sys.exit(1)
    return p_state


def parse_state(path, n_active, K_, use_pair, use_floor):
    """Parse the state file into (state_per_step, traj, final_v)."""
    raw = np.fromfile(path, dtype=np.float32)
    extra = (1 if use_pair else 0) + (1 if use_floor else 0)
    per_step = n_active * (3 + 3 + 1 + 3 + 1 + extra)
    state_n = K_ * per_step
    traj_n = (K_ + 1) * n_active * 3
    final_v_n = n_active * 3
    expected_n = state_n + traj_n + final_v_n
    if raw.size != expected_n:
        return None, None, None, raw.size, expected_n
    state = raw[:state_n].reshape(K_, per_step)
    traj = raw[state_n:state_n + traj_n].reshape(K_ + 1, n_active, 3)
    final_v = raw[state_n + traj_n:].reshape(n_active, 3)
    return state, traj, final_v, raw.size, expected_n


def expected_size_bytes(n_active, K_, use_pair, use_floor):
    extra = (1 if use_pair else 0) + (1 if use_floor else 0)
    per_step = n_active * (3 + 3 + 1 + 3 + 1 + extra)
    floats = K_ * per_step + (K_ + 1) * n_active * 3 + n_active * 3
    return 4 * floats


def check(label, use_pair, use_floor, state_path):
    file_bytes = os.path.getsize(state_path)
    expected = expected_size_bytes(N_ACTIVE, K, use_pair, use_floor)
    print(f"[{label}] file size: {file_bytes} bytes (expected {expected})")
    if file_bytes != expected:
        print(f"[{label}] [FAIL] size mismatch")
        return False
    state, traj, final_v, got, exp = parse_state(
        state_path, N_ACTIVE, K, use_pair, use_floor)
    if state is None:
        print(f"[{label}] [FAIL] parse: got {got} floats, expected {exp}")
        return False
    # No NaN/Inf anywhere.
    raw = np.fromfile(state_path, dtype=np.float32)
    if np.isnan(raw).any() or np.isinf(raw).any():
        print(f"[{label}] [FAIL] NaN/Inf in state file")
        return False
    # Trajectory must evolve — positions shouldn't all be identical across
    # frames. Compare any two consecutive frames; total motion > eps.
    motion = float(np.linalg.norm(traj[-1] - traj[0]))
    print(f"[{label}] total trajectory motion: {motion:.6f}")
    if motion < 1e-8:
        print(f"[{label}] [FAIL] trajectory frames are identical (no evolution)")
        return False
    print(f"[{label}] [PASS]")
    return True


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1

    print("Test: xpbd_full_fwd structural smoke test\n")

    overall = True

    # (a) Required args only.
    label = "a_required"
    print(f"--- Scenario (a): required args only ---")
    p = run_full_fwd([], label)
    overall &= check(label, use_pair=False, use_floor=False, state_path=p)

    # (b) visc_pair_coef > 0. Need sim_scale (arg 15) before visc_pair_coef.
    label = "b_pair"
    print(f"\n--- Scenario (b): visc_pair_coef = 1e-4 ---")
    p = run_full_fwd([1.0, 1e-4], label)  # sim_scale=1.0, visc=1e-4
    overall &= check(label, use_pair=True, use_floor=False, state_path=p)

    # (c) spring_K > 0 with a single bond. Arg layout: sim_scale=1.0,
    # visc=0, spring_K=100, bonds.bin.
    label = "c_springs"
    print(f"\n--- Scenario (c): spring_K = 100 with 1 bond ---")
    bonds = [(0, 1, 0.5)]  # particles 0 and 1 are at distance 0.5 in scenario
    p_bonds = os.path.join(TMP, f"cuda_xfull_{label}_bonds.bin")
    write_bonds_bin(p_bonds, bonds)
    p = run_full_fwd([1.0, 0.0, 100.0, p_bonds], label)
    overall &= check(label, use_pair=False, use_floor=False, state_path=p)

    # (d) floor_y = 0.0 + restitution = 0.5. Arg layout requires sim_scale,
    # visc, spring_K, bonds.bin slots to be filled (even if unused) — Metal
    # uses positional args. Pass spring_K=0 and a placeholder bonds path.
    label = "d_floor"
    print(f"\n--- Scenario (d): floor_y = 0.0, restitution = 0.5 ---")
    p_dummy_bonds = os.path.join(TMP, f"cuda_xfull_{label}_dummy_bonds.bin")
    # Write empty bonds file so the slot is positionally filled. Driver
    # only loads it when spring_K > 0, so an empty file is fine.
    open(p_dummy_bonds, "wb").close()
    p = run_full_fwd([1.0, 0.0, 0.0, p_dummy_bonds, 0.0, 0.5], label)
    overall &= check(label, use_pair=False, use_floor=True, state_path=p)

    print()
    if overall:
        print("[OVERALL PASS] xpbd_full_fwd structural smoke test")
        return 0
    print("[OVERALL FAIL]")
    return 1


if __name__ == "__main__":
    sys.exit(main())
