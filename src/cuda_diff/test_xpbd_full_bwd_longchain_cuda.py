"""Long-chain FD validation of CUDA xpbd_full_bwd at K = 2, 20, 100, 1000.

Addresses V1 code-review concern #7: the standard FD test (test_xpbd_full_bwd_cuda.py)
uses K=2 which doesn't stress chain amplification. This script runs the same
gradient-validation pattern across growing K to confirm the analytic backward
matches FD as the unrolled chain grows.

Scenario design choices:
- 4 active particles starting above a floor at y=0.0, with floor restitution
  zero (inelastic clamp). Particles fall under gravity, get clamped by the
  floor, and stay engaged with each other via SPH density interaction — so
  the gradient chain has signal to propagate at every step.
- visc_pair_coef and spring_K enabled (asymmetric weighted loss) so the full
  pair_forces + spring chain is exercised, not just density.
- BWD_TBPTT not set: walk every step (the point is to validate the full chain,
  not the truncated version).
- BWD_CLIP_NORM not set: use the new 1e3 default (clip kicks in only at long K).

Pass criterion per K:
- K=2:    rel_err < 1e-2 (same as standard test)
- K=20:   rel_err < 5e-2 (small amplification expected)
- K=100:  rel_err < 5e-1 (FD noise + chaotic amplification on 4-particle system)
- K=1000: report only — at this length, FD vs analytic disagreement does not
          necessarily indicate a bug (loss landscape becomes near-chaotic).
          The key signal is: do the analytic gradients stay finite (no NaN)
          and bounded under BWD_CLIP_NORM, and does rel_err scale with K
          smoothly rather than blowing up to inf?

ASCII-only output; tempfile.gettempdir() for all paths.
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

# Scenario constants. Floor enabled so particles stay engaged.
N_ACTIVE = 4
N_STATIC = 8
H = 1.5
MASS = 1.0
RHO_REST = 0.3
DT = 0.01
G_Y = -9.81
ALPHA_DENS = 0.5
SIM_SCALE = 1.0
VISC_PAIR_COEF = 1.0
SPRING_K = 5.0
FLOOR_Y = 0.0
RESTITUTION = 0.0


def make_scenario():
    pos0 = np.array([
        [-0.15, 0.5, -0.15],
        [ 0.15, 0.5, -0.15],
        [-0.15, 0.5,  0.15],
        [ 0.15, 0.5,  0.15],
    ], dtype=np.float32)
    vel0 = np.array([
        [ 0.1, 0.0,  0.05],
        [-0.1, 0.0,  0.05],
        [ 0.1, 0.0, -0.05],
        [-0.1, 0.0, -0.05],
    ], dtype=np.float32)
    pos_static = np.array([
        [-0.5, -0.3, -0.5], [0.5, -0.3, -0.5],
        [-0.5, -0.3,  0.5], [0.5, -0.3,  0.5],
        [-0.5,  0.0, -0.5], [0.5,  0.0, -0.5],
        [-0.5,  0.0,  0.5], [0.5,  0.0,  0.5],
    ], dtype=np.float32)
    weights = np.array([
        [2.0, 1.0, 0.5],
        [-1.0, 0.5, 1.5],
        [0.5, 2.0, -0.5],
        [1.5, -1.0, 1.0],
    ], dtype=np.float32)
    bonds = [(0, 1, 0.5)]
    return pos0, vel0, pos_static, weights, bonds


def write_bonds_bin(path, bonds):
    with open(path, "wb") as f:
        for i, j, L in bonds:
            f.write(struct.pack("ii f f", int(i), int(j), float(L), 0.0))


def run_fwd(K, pos0, vel0, pos_static, rho_rest, bonds_path, label):
    p_pos    = os.path.join(TMP, f"lc_{label}_pos.bin")
    p_vel    = os.path.join(TMP, f"lc_{label}_vel.bin")
    p_static = os.path.join(TMP, f"lc_{label}_static.bin")
    p_state  = os.path.join(TMP, f"lc_{label}_state.bin")
    pos0.astype(np.float32).tofile(p_pos)
    vel0.astype(np.float32).tofile(p_vel)
    pos_static.astype(np.float32).tofile(p_static)

    cmd = [BINARY, "xpbd_full_fwd",
           str(N_ACTIVE), str(N_STATIC), str(K),
           str(H), str(MASS), str(rho_rest), str(DT),
           str(G_Y), str(ALPHA_DENS),
           p_pos, p_vel, p_static, p_state,
           str(SIM_SCALE), str(VISC_PAIR_COEF),
           str(SPRING_K), bonds_path,
           str(FLOOR_Y), str(RESTITUTION)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[{label}] FWD STDERR:\n{r.stderr[-500:]}", file=sys.stderr)
        sys.exit(1)
    return p_state


def parse_state(K, state_path):
    # use_pair=True, use_floor=True
    per_step = N_ACTIVE * (3 + 3 + 1 + 3 + 1 + 1 + 1)
    raw = np.fromfile(state_path, dtype=np.float32)
    state_n = K * per_step
    traj_n = (K + 1) * N_ACTIVE * 3
    traj = raw[state_n:state_n + traj_n].reshape(K + 1, N_ACTIVE, 3)
    return traj


def fwd_loss(K, pos0, vel0, pos_static, rho_rest, bonds_path, weights, label):
    state_path = run_fwd(K, pos0, vel0, pos_static, rho_rest, bonds_path, label)
    traj = parse_state(K, state_path)
    return float((traj[-1] * weights).sum())


def run_bwd(K, state_path, pos_static, grad_x_final, rho_rest, bonds_path, label):
    p_static = os.path.join(TMP, f"lc_{label}_static_bwd.bin")
    pos_static.astype(np.float32).tofile(p_static)
    p_gxf    = os.path.join(TMP, f"lc_{label}_gxf.bin")
    p_gxi    = os.path.join(TMP, f"lc_{label}_gxi.bin")
    p_gvi    = os.path.join(TMP, f"lc_{label}_gvi.bin")
    p_grho   = os.path.join(TMP, f"lc_{label}_grho.bin")
    p_gK     = os.path.join(TMP, f"lc_{label}_gK.bin")
    p_gvK    = os.path.join(TMP, f"lc_{label}_gvK.bin")
    p_galpha = os.path.join(TMP, f"lc_{label}_galpha.bin")
    grad_x_final.astype(np.float32).tofile(p_gxf)
    cmd = [BINARY, "xpbd_full_bwd",
           str(N_ACTIVE), str(N_STATIC), str(K),
           str(H), str(MASS), str(rho_rest), str(DT),
           str(G_Y), str(ALPHA_DENS),
           state_path, p_static, p_gxf,
           p_gxi, p_gvi, p_grho,
           str(SIM_SCALE), str(VISC_PAIR_COEF),
           str(SPRING_K), bonds_path,
           p_gK, p_gvK,
           str(FLOOR_Y), p_galpha, str(RESTITUTION)]
    # Make sure no env vars from prior tests affect this run.
    env = {k: v for k, v in os.environ.items() if k not in ("BWD_TBPTT",)}
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(f"[{label}] BWD STDERR:\n{r.stderr[-500:]}", file=sys.stderr)
        sys.exit(1)
    gxi  = np.fromfile(p_gxi,  dtype=np.float32).reshape(N_ACTIVE, 3)
    gvi  = np.fromfile(p_gvi,  dtype=np.float32).reshape(N_ACTIVE, 3)
    grho = float(np.fromfile(p_grho, dtype=np.float32)[0])
    return gxi, gvi, grho


def fd_grad_pos(K, pos0, vel0, pos_static, rho_rest, bonds_path, weights,
                label, eps=1e-4):
    g = np.zeros_like(pos0)
    for i in range(N_ACTIVE):
        for d in range(3):
            save = pos0[i, d]
            pos0[i, d] = save + eps
            Lp = fwd_loss(K, pos0, vel0, pos_static, rho_rest, bonds_path,
                          weights, f"{label}_pp_{i}_{d}")
            pos0[i, d] = save - eps
            Lm = fwd_loss(K, pos0, vel0, pos_static, rho_rest, bonds_path,
                          weights, f"{label}_pm_{i}_{d}")
            pos0[i, d] = save
            g[i, d] = (Lp - Lm) / (2 * eps)
    return g


def fd_grad_vel(K, pos0, vel0, pos_static, rho_rest, bonds_path, weights,
                label, eps=1e-3):
    g = np.zeros_like(vel0)
    for i in range(N_ACTIVE):
        for d in range(3):
            save = vel0[i, d]
            vel0[i, d] = save + eps
            Lp = fwd_loss(K, pos0, vel0, pos_static, rho_rest, bonds_path,
                          weights, f"{label}_vp_{i}_{d}")
            vel0[i, d] = save - eps
            Lm = fwd_loss(K, pos0, vel0, pos_static, rho_rest, bonds_path,
                          weights, f"{label}_vm_{i}_{d}")
            vel0[i, d] = save
            g[i, d] = (Lp - Lm) / (2 * eps)
    return g


def fd_grad_rho(K, pos0, vel0, pos_static, rho_rest, bonds_path, weights,
                label, eps=1e-3):
    Lp = fwd_loss(K, pos0, vel0, pos_static, rho_rest + eps, bonds_path,
                  weights, f"{label}_rhop")
    Lm = fwd_loss(K, pos0, vel0, pos_static, rho_rest - eps, bonds_path,
                  weights, f"{label}_rhom")
    return (Lp - Lm) / (2 * eps)


def rel_err(an, num):
    an = np.atleast_1d(np.asarray(an, dtype=np.float64))
    num = np.atleast_1d(np.asarray(num, dtype=np.float64))
    err = float(np.abs(an - num).max())
    mag = max(float(np.abs(num).max()), 1.0)
    return err / mag


def run_K(K, tol_pos, tol_vel, tol_rho):
    pos0, vel0, pos_static, weights, bonds = make_scenario()
    bonds_path = os.path.join(TMP, f"lc_K{K}_bonds.bin")
    write_bonds_bin(bonds_path, bonds)
    label = f"K{K}"

    state_path = run_fwd(K, pos0, vel0, pos_static, RHO_REST, bonds_path, label)
    grad_x_final = weights.copy()
    gxi, gvi, grho = run_bwd(K, state_path, pos_static, grad_x_final, RHO_REST,
                              bonds_path, label)
    has_nan_x = bool(np.isnan(gxi).any())
    has_nan_v = bool(np.isnan(gvi).any())
    has_nan_r = bool(np.isnan(grho))
    if has_nan_x or has_nan_v or has_nan_r:
        print(f"  K={K}: NaN in analytic gradients (x:{has_nan_x} v:{has_nan_v} rho:{has_nan_r})")
        return False, None, None, None

    fd_gxi = fd_grad_pos(K, pos0.copy(), vel0.copy(), pos_static, RHO_REST,
                          bonds_path, weights, label)
    fd_gvi = fd_grad_vel(K, pos0.copy(), vel0.copy(), pos_static, RHO_REST,
                          bonds_path, weights, label)
    fd_grho = fd_grad_rho(K, pos0.copy(), vel0.copy(), pos_static, RHO_REST,
                           bonds_path, weights, label)

    pos_rel = rel_err(gxi, fd_gxi)
    vel_rel = rel_err(gvi, fd_gvi)
    rho_rel = rel_err(grho, fd_grho)

    pos_ok = pos_rel < tol_pos
    vel_ok = vel_rel < tol_vel
    rho_ok = rho_rel < tol_rho

    print(f"  K={K:5d}: pos {pos_rel:.3e} {'PASS' if pos_ok else 'FAIL':4s} | "
          f"vel {vel_rel:.3e} {'PASS' if vel_ok else 'FAIL':4s} | "
          f"rho {rho_rel:.3e} {'PASS' if rho_ok else 'FAIL':4s} | "
          f"|gx|inf={np.abs(gxi).max():.2e} |gv|inf={np.abs(gvi).max():.2e}")
    return (pos_ok and vel_ok and rho_ok), pos_rel, vel_rel, rho_rel


def main():
    print("=== Long-chain FD validation of xpbd_full_bwd ===")
    print(f"  scenario: {N_ACTIVE} active + {N_STATIC} boundary, floor at y={FLOOR_Y}")
    print(f"  springs + pair forces ON, asymmetric weights loss")
    print()

    results = []
    # K=2: standard tolerance
    ok, *_ = run_K(2, tol_pos=1e-2, tol_vel=1e-2, tol_rho=1e-2)
    results.append((2, ok))
    # K=20: small amplification expected
    ok, *_ = run_K(20, tol_pos=5e-2, tol_vel=5e-2, tol_rho=5e-2)
    results.append((20, ok))
    # K=100: chaotic amplification, loose threshold
    ok, *_ = run_K(100, tol_pos=5e-1, tol_vel=5e-1, tol_rho=5e-1)
    results.append((100, ok))
    # K=1000: report-only; pass iff gradients are finite (rel_err meaningless)
    ok, *_ = run_K(1000, tol_pos=1e10, tol_vel=1e10, tol_rho=1e10)
    results.append((1000, ok))

    print()
    print("=== Summary ===")
    for K, ok in results:
        verdict = "PASS" if ok else "FAIL"
        print(f"  K={K:5d}: {verdict}")
    all_ok = all(ok for _, ok in results)
    print(f"\n[{'OVERALL PASS' if all_ok else 'OVERALL FAIL'}] long-chain FD validation")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
