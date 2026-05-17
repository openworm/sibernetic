"""FD validation of the CUDA xpbd_full_bwd reverse-mode driver.

Two scenarios, both tiny enough for central-difference FD over the full
input space:

(1) Vanilla XPBD: 4 active + 8 boundary particles, dropping under gravity.
    Validates grad_x_init, grad_v_init, grad_rho_rest against centered FD.

(2) Spring + pair forces: same geometry, with visc_pair_coef=1e-3,
    spring_K=10, one bond (0,1, rest 0.5). Validates the additional
    grad_spring_K and grad_visc_K outputs.

Threshold: 1e-2 relative error (matches density chain tests).

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

# ── Scenario constants (kept tiny so FD over 4*3+4*3+1=25 inputs is cheap) ─
# Particles are clustered close together with rho_rest BELOW the implied
# initial density so the C > 0 branch of solve_density_constraint fires
# every step (otherwise grad_rho is trivially 0).
N_ACTIVE = 4
N_STATIC = 8
K = 2
H = 1.5
MASS = 1.0
RHO_REST = 0.3
DT = 0.01
G_Y = -9.81
ALPHA_DENS = 0.5


def make_scenario():
    # Particles overlap kernel support (h=1.5) so density is well above
    # rho_rest=0.3, ensuring C > 0 and the density projection actually
    # moves particles per step.
    pos0 = np.array([
        [-0.15, 1.0, -0.15],
        [ 0.15, 1.0, -0.15],
        [-0.15, 1.0,  0.15],
        [ 0.15, 1.0,  0.15],
    ], dtype=np.float32)
    # Non-zero initial velocity so viscosity actually fires per step (visc
    # is proportional to relative velocity between neighbors).
    vel0 = np.array([
        [ 0.3, 0.0,  0.2],
        [-0.3, 0.0,  0.2],
        [ 0.3, 0.0, -0.2],
        [-0.3, 0.0, -0.2],
    ], dtype=np.float32)
    pos_static = np.array([
        [-0.5, 0.0, -0.5], [0.5, 0.0, -0.5],
        [-0.5, 0.0,  0.5], [0.5, 0.0,  0.5],
        [-0.5, 0.4, -0.5], [0.5, 0.4, -0.5],
        [-0.5, 0.4,  0.5], [0.5, 0.4,  0.5],
    ], dtype=np.float32)
    return pos0, vel0, pos_static


def write_bonds_bin(path, bonds):
    with open(path, "wb") as f:
        for i, j, L in bonds:
            f.write(struct.pack("ii f f", int(i), int(j), float(L), 0.0))


def run_fwd(pos0, vel0, pos_static, rho_rest, extra_args, label):
    """Run xpbd_full_fwd with the given extra optional args (after argv[14])."""
    p_pos    = os.path.join(TMP, f"cuda_bwd_{label}_pos.bin")
    p_vel    = os.path.join(TMP, f"cuda_bwd_{label}_vel.bin")
    p_static = os.path.join(TMP, f"cuda_bwd_{label}_static.bin")
    p_state  = os.path.join(TMP, f"cuda_bwd_{label}_state.bin")
    pos0.astype(np.float32).tofile(p_pos)
    vel0.astype(np.float32).tofile(p_vel)
    pos_static.astype(np.float32).tofile(p_static)

    cmd = [BINARY, "xpbd_full_fwd",
           str(N_ACTIVE), str(N_STATIC), str(K),
           str(H), str(MASS), str(rho_rest), str(DT),
           str(G_Y), str(ALPHA_DENS),
           p_pos, p_vel, p_static, p_state] + [str(a) for a in extra_args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[{label}] FWD STDERR:\n{r.stderr}", file=sys.stderr)
        sys.exit(1)
    return p_state, p_static


def parse_state(state_path, use_pair, use_floor):
    extra = (1 if use_pair else 0) + (1 if use_floor else 0)
    per_step = N_ACTIVE * (3 + 3 + 1 + 3 + 1 + extra)
    raw = np.fromfile(state_path, dtype=np.float32)
    state_n = K * per_step
    traj_n = (K + 1) * N_ACTIVE * 3
    traj = raw[state_n:state_n + traj_n].reshape(K + 1, N_ACTIVE, 3)
    return traj


def fwd_loss(pos0, vel0, pos_static, rho_rest, extra_args, label,
             use_pair=False, use_floor=False, weights=None):
    """Returns L = sum(weights · x_final). Defaults to weights = 1 everywhere."""
    state_path, _ = run_fwd(pos0, vel0, pos_static, rho_rest, extra_args, label)
    traj = parse_state(state_path, use_pair, use_floor)
    if weights is None:
        return float(traj[-1].sum())
    return float((traj[-1] * weights).sum())


def run_bwd(state_path, pos_static, grad_x_final, rho_rest, extra_bwd_args,
            label):
    p_static = os.path.join(TMP, f"cuda_bwd_{label}_static_bwd.bin")
    pos_static.astype(np.float32).tofile(p_static)
    p_gxf  = os.path.join(TMP, f"cuda_bwd_{label}_gxf.bin")
    p_gxi  = os.path.join(TMP, f"cuda_bwd_{label}_gxi.bin")
    p_gvi  = os.path.join(TMP, f"cuda_bwd_{label}_gvi.bin")
    p_grho = os.path.join(TMP, f"cuda_bwd_{label}_grho.bin")
    grad_x_final.astype(np.float32).tofile(p_gxf)

    cmd = [BINARY, "xpbd_full_bwd",
           str(N_ACTIVE), str(N_STATIC), str(K),
           str(H), str(MASS), str(rho_rest), str(DT),
           str(G_Y), str(ALPHA_DENS),
           state_path, p_static, p_gxf,
           p_gxi, p_gvi, p_grho] + [str(a) for a in extra_bwd_args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[{label}] BWD STDERR:\n{r.stderr}", file=sys.stderr)
        sys.exit(1)
    gxi  = np.fromfile(p_gxi,  dtype=np.float32).reshape(N_ACTIVE, 3)
    gvi  = np.fromfile(p_gvi,  dtype=np.float32).reshape(N_ACTIVE, 3)
    grho = float(np.fromfile(p_grho, dtype=np.float32)[0])
    return gxi, gvi, grho


def fd_grad_pos(pos0, vel0, pos_static, rho_rest, extra_args, label,
                use_pair=False, use_floor=False, eps=1e-4, weights=None):
    g = np.zeros_like(pos0)
    for i in range(N_ACTIVE):
        for d in range(3):
            save = pos0[i, d]
            pos0[i, d] = save + eps
            Lp = fwd_loss(pos0, vel0, pos_static, rho_rest, extra_args,
                          f"{label}_pp_{i}_{d}", use_pair, use_floor, weights)
            pos0[i, d] = save - eps
            Lm = fwd_loss(pos0, vel0, pos_static, rho_rest, extra_args,
                          f"{label}_pm_{i}_{d}", use_pair, use_floor, weights)
            pos0[i, d] = save
            g[i, d] = (Lp - Lm) / (2 * eps)
    return g


def fd_grad_vel(pos0, vel0, pos_static, rho_rest, extra_args, label,
                use_pair=False, use_floor=False, eps=1e-3, weights=None):
    g = np.zeros_like(vel0)
    for i in range(N_ACTIVE):
        for d in range(3):
            save = vel0[i, d]
            vel0[i, d] = save + eps
            Lp = fwd_loss(pos0, vel0, pos_static, rho_rest, extra_args,
                          f"{label}_vp_{i}_{d}", use_pair, use_floor, weights)
            vel0[i, d] = save - eps
            Lm = fwd_loss(pos0, vel0, pos_static, rho_rest, extra_args,
                          f"{label}_vm_{i}_{d}", use_pair, use_floor, weights)
            vel0[i, d] = save
            g[i, d] = (Lp - Lm) / (2 * eps)
    return g


def fd_grad_rho(pos0, vel0, pos_static, rho_rest, extra_args, label,
                use_pair=False, use_floor=False, eps=1e-3, weights=None):
    Lp = fwd_loss(pos0, vel0, pos_static, rho_rest + eps, extra_args,
                  f"{label}_rhop", use_pair, use_floor, weights)
    Lm = fwd_loss(pos0, vel0, pos_static, rho_rest - eps, extra_args,
                  f"{label}_rhom", use_pair, use_floor, weights)
    return (Lp - Lm) / (2 * eps)


def rel_err(an, num):
    an = np.atleast_1d(np.asarray(an, dtype=np.float64))
    num = np.atleast_1d(np.asarray(num, dtype=np.float64))
    err = float(np.abs(an - num).max())
    mag = max(float(np.abs(num).max()), 1.0)
    return err / mag


def test_vanilla(tol=1e-2):
    print("=== Scenario 1: vanilla XPBD (no pair, no spring, no floor) ===")
    pos0, vel0, pos_static = make_scenario()
    label = "s1"

    # FWD pass to produce state.
    state_path, _ = run_fwd(pos0, vel0, pos_static, RHO_REST, [], label)
    # grad_x_final = 1 (L = sum(x_final))
    grad_x_final = np.ones((N_ACTIVE, 3), dtype=np.float32)

    # Analytic gradients.
    gxi, gvi, grho = run_bwd(state_path, pos_static, grad_x_final, RHO_REST,
                              [], label)

    # FD gradients.
    fd_gxi = fd_grad_pos(pos0.copy(), vel0.copy(), pos_static, RHO_REST, [],
                          label)
    fd_gvi = fd_grad_vel(pos0.copy(), vel0.copy(), pos_static, RHO_REST, [],
                          label)
    fd_grho = fd_grad_rho(pos0.copy(), vel0.copy(), pos_static, RHO_REST, [],
                           label)

    pos_rel = rel_err(gxi, fd_gxi)
    vel_rel = rel_err(gvi, fd_gvi)
    rho_rel = rel_err(grho, fd_grho)

    pos_ok = pos_rel < tol
    vel_ok = vel_rel < tol
    rho_ok = rho_rel < tol

    print(f"  pos rel_err: {pos_rel:.3e}  {'[PASS]' if pos_ok else '[FAIL]'}")
    print(f"  vel rel_err: {vel_rel:.3e}  {'[PASS]' if vel_ok else '[FAIL]'}")
    print(f"  rho rel_err: {rho_rel:.3e}  {'[PASS]' if rho_ok else '[FAIL]'}")
    if not (pos_ok and vel_ok and rho_ok):
        print(f"  analytic grad_x[0]: {gxi[0]}")
        print(f"  FD       grad_x[0]: {fd_gxi[0]}")
        print(f"  analytic grad_v[0]: {gvi[0]}")
        print(f"  FD       grad_v[0]: {fd_gvi[0]}")
        print(f"  analytic grad_rho:  {grho}")
        print(f"  FD       grad_rho:  {fd_grho}")
    return pos_ok and vel_ok and rho_ok


def test_spring_visc(tol=1e-2):
    print("=== Scenario 2: spring_K + visc_pair_coef ===")
    pos0, vel0, pos_static = make_scenario()
    label = "s2"
    sim_scale = 1.0
    # Use moderately large coefficients so both spring and viscosity have
    # measurable effects on the loss; large enough to lift visc FD out of
    # float32 noise, small enough to keep the integration well-behaved.
    visc_pair_coef = 1e2
    spring_K = 10.0
    # Particles 0 and 1 are 0.30 apart; rest_len=0.5 stretches the bond.
    bonds = [(0, 1, 0.5)]

    p_bonds = os.path.join(TMP, f"cuda_bwd_{label}_bonds.bin")
    write_bonds_bin(p_bonds, bonds)
    extra_fwd = [sim_scale, visc_pair_coef, spring_K, p_bonds]
    # FWD pass to produce state.
    state_path, _ = run_fwd(pos0, vel0, pos_static, RHO_REST, extra_fwd, label)

    # Asymmetric grad_x_final / weights — required for spring_K and visc_K
    # gradients to be non-zero. A uniform weights=1 loss is symmetric under
    # bond pull-together, so ∂L/∂spring_K trivially cancels.
    weights = np.array([
        [2.0, 1.0, 0.5],
        [-1.0, 0.5, 1.5],
        [0.5, 2.0, -0.5],
        [1.5, -1.0, 1.0],
    ], dtype=np.float32)
    grad_x_final = weights.copy()

    p_gspring = os.path.join(TMP, f"cuda_bwd_{label}_gspring.bin")
    p_gvisc   = os.path.join(TMP, f"cuda_bwd_{label}_gvisc.bin")
    extra_bwd = [sim_scale, visc_pair_coef, spring_K, p_bonds,
                 p_gspring, p_gvisc]
    gxi, gvi, grho = run_bwd(state_path, pos_static, grad_x_final, RHO_REST,
                              extra_bwd, label)
    g_spring = float(np.fromfile(p_gspring, dtype=np.float32)[0])
    g_visc   = float(np.fromfile(p_gvisc,   dtype=np.float32)[0])

    # FD gradients for pos/vel/rho (using same weighted loss).
    fd_gxi = fd_grad_pos(pos0.copy(), vel0.copy(), pos_static, RHO_REST,
                          extra_fwd, label, use_pair=True, weights=weights)
    fd_gvi = fd_grad_vel(pos0.copy(), vel0.copy(), pos_static, RHO_REST,
                          extra_fwd, label, use_pair=True, weights=weights)
    fd_grho = fd_grad_rho(pos0.copy(), vel0.copy(), pos_static, RHO_REST,
                           extra_fwd, label, use_pair=True, weights=weights)

    # FD for spring_K (weighted loss).
    def fwd_with_spring(sk):
        ext = [sim_scale, visc_pair_coef, float(sk), p_bonds]
        return fwd_loss(pos0.copy(), vel0.copy(), pos_static, RHO_REST,
                         ext, f"{label}_sk_{sk}", use_pair=True,
                         weights=weights)

    eps_sk = 1e-1
    L_p = fwd_with_spring(spring_K + eps_sk)
    L_m = fwd_with_spring(spring_K - eps_sk)
    fd_gspring = (L_p - L_m) / (2 * eps_sk)

    # FD for visc_pair_coef (weighted loss). visc effect on L is small
    # for our toy scenario — use larger eps to lift |L_p - L_m| out of
    # float32 noise.
    def fwd_with_visc(vk):
        ext = [sim_scale, float(vk), spring_K, p_bonds]
        return fwd_loss(pos0.copy(), vel0.copy(), pos_static, RHO_REST,
                         ext, f"{label}_vk_{vk}", use_pair=True,
                         weights=weights)

    eps_vk = 1.0
    L_p = fwd_with_visc(visc_pair_coef + eps_vk)
    L_m = fwd_with_visc(visc_pair_coef - eps_vk)
    fd_gvisc = (L_p - L_m) / (2 * eps_vk)

    pos_rel = rel_err(gxi, fd_gxi)
    vel_rel = rel_err(gvi, fd_gvi)
    rho_rel = rel_err(grho, fd_grho)
    spring_rel = abs(g_spring - fd_gspring) / max(abs(fd_gspring), 1.0)
    visc_rel = abs(g_visc - fd_gvisc) / max(abs(fd_gvisc), 1.0)

    pos_ok = pos_rel < tol
    vel_ok = vel_rel < tol
    rho_ok = rho_rel < tol
    sp_ok  = spring_rel < tol
    vk_ok  = visc_rel < tol

    print(f"  pos rel_err:      {pos_rel:.3e}  {'[PASS]' if pos_ok else '[FAIL]'}")
    print(f"  vel rel_err:      {vel_rel:.3e}  {'[PASS]' if vel_ok else '[FAIL]'}")
    print(f"  rho rel_err:      {rho_rel:.3e}  {'[PASS]' if rho_ok else '[FAIL]'}")
    print(f"  spring_K rel_err: {spring_rel:.3e}  "
          f"{'[PASS]' if sp_ok else '[FAIL]'}  "
          f"(analytic={g_spring:+.4e}, FD={fd_gspring:+.4e})")
    print(f"  visc_K rel_err:   {visc_rel:.3e}  "
          f"{'[PASS]' if vk_ok else '[FAIL]'}  "
          f"(analytic={g_visc:+.4e}, FD={fd_gvisc:+.4e})")
    if not (pos_ok and vel_ok and rho_ok and sp_ok and vk_ok):
        print(f"  analytic grad_x[0]: {gxi[0]}")
        print(f"  FD       grad_x[0]: {fd_gxi[0]}")
        print(f"  analytic grad_v[0]: {gvi[0]}")
        print(f"  FD       grad_v[0]: {fd_gvi[0]}")
        print(f"  analytic grad_rho:  {grho}")
        print(f"  FD       grad_rho:  {fd_grho}")
    return pos_ok and vel_ok and rho_ok and sp_ok and vk_ok


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1

    ok1 = test_vanilla()
    print()
    ok2 = test_spring_visc()
    print()
    if ok1 and ok2:
        print("[OVERALL PASS] xpbd_full_bwd FD validation")
        return 0
    print("[OVERALL FAIL]")
    return 1


if __name__ == "__main__":
    sys.exit(main())
