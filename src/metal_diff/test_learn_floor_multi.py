"""Multi-particle gradient learning: 64-particle cube falls, learn floor_y.

No bond physics yet (distance constraint backward is the next milestone).
This validates that the M7.D/E framework scales from 1 to 64 particles
and that the gradient on a global parameter (floor_y) accumulates
correctly across all particles' clamp events over many timesteps.

Compare to test_learn_floor.py (1 particle) to confirm the framework
generalizes without changes to the kernels.
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def fwd(x0, v0, K, dt, g_y, fl_y):
    n = len(x0)
    p_x = "/tmp/lfm_x0.bin"; p_v = "/tmp/lfm_v0.bin"
    p_tp = "/tmp/lfm_tpos.bin"; p_tc = "/tmp/lfm_tclmp.bin"
    p_vf = "/tmp/lfm_vfinal.bin"
    x0.astype(np.float32).tofile(p_x)
    v0.astype(np.float32).tofile(p_v)
    subprocess.run(
        [BINARY, "step_floor_fwd", str(n), str(K),
         str(dt), str(g_y), str(fl_y),
         p_x, p_v, p_tp, p_tc, p_vf],
        check=True,
    )
    traj_pos = np.fromfile(p_tp, dtype=np.float32).reshape(K + 1, n, 3)
    traj_clamp = np.fromfile(p_tc, dtype=np.int32).reshape(K, n)
    v_final = np.fromfile(p_vf, dtype=np.float32).reshape(n, 3)
    x_final = traj_pos[-1]
    return x_final, v_final, traj_pos, traj_clamp


def bwd(grad_x_final, traj_pos, traj_clamp, dt):
    K = traj_clamp.shape[0]
    n = traj_clamp.shape[1]
    p_tp = "/tmp/lfm_tpos.bin"; p_tc = "/tmp/lfm_tclmp.bin"
    p_gxf = "/tmp/lfm_gxfin.bin"
    p_gxi = "/tmp/lfm_gxinit.bin"; p_gvi = "/tmp/lfm_gvinit.bin"
    p_ggy = "/tmp/lfm_ggy.bin"; p_gfy = "/tmp/lfm_gfy.bin"
    traj_pos.astype(np.float32).tofile(p_tp)
    traj_clamp.astype(np.int32).tofile(p_tc)
    grad_x_final.astype(np.float32).tofile(p_gxf)
    subprocess.run(
        [BINARY, "step_floor_bwd", str(n), str(K), str(dt),
         p_tp, p_tc, p_gxf, p_gxi, p_gvi, p_ggy, p_gfy],
        check=True,
    )
    g_x = np.fromfile(p_gxi, dtype=np.float32).reshape(n, 3)
    g_v = np.fromfile(p_gvi, dtype=np.float32).reshape(n, 3)
    g_g = float(np.fromfile(p_ggy, dtype=np.float32)[0])
    g_f = float(np.fromfile(p_gfy, dtype=np.float32)[0])
    return g_x, g_v, g_g, g_f


def make_cube(n_side, spacing, center_y):
    n = n_side ** 3
    pos = np.zeros((n, 3), dtype=np.float32)
    half = (n_side - 1) * spacing / 2
    for ix in range(n_side):
        for iy in range(n_side):
            for iz in range(n_side):
                p = ix * n_side * n_side + iy * n_side + iz
                pos[p, 0] = ix * spacing - half
                pos[p, 1] = center_y + iy * spacing
                pos[p, 2] = iz * spacing - half
    return pos


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    # 64-particle cube, drop from height, clamp on floor.
    n_side = 4
    n = n_side ** 3
    K = 600           # 3.0s sim — long enough to fully settle
    dt = 0.005
    g_y = -9.8
    true_floor = 0.5
    spacing = 0.6

    x0 = make_cube(n_side, spacing, center_y=5.0)
    v0 = np.zeros_like(x0)

    # Generate "observed" final positions with the true floor.
    x_obs, _, _, _ = fwd(x0, v0, K, dt, g_y, true_floor)
    print(f"=== Multi-particle (n={n}) gradient learning ===")
    print(f"  Setup: 4×4×4 cube dropped from y=5, K={K} steps × dt={dt}={K*dt:.1f}s sim")
    print(f"  True floor_y = {true_floor:.3f}")
    print(f"  Observed final state: mean_y={x_obs[:,1].mean():.4f}  "
          f"min_y={x_obs[:,1].min():.4f}  max_y={x_obs[:,1].max():.4f}")
    print()

    # ── M7.D scale-up: gradient check on multi-particle ──
    bad_floor = 1.0
    x_bad, _, traj_pos, traj_clamp = fwd(x0, v0, K, dt, g_y, bad_floor)
    grad_x_final = (x_bad - x_obs).astype(np.float32)
    _, _, _, g_floor_kernel = bwd(grad_x_final, traj_pos, traj_clamp, dt)

    eps = 1e-3
    x_plus, _, _, _ = fwd(x0, v0, K, dt, g_y, bad_floor + eps)
    x_minus, _, _, _ = fwd(x0, v0, K, dt, g_y, bad_floor - eps)
    L_plus = 0.5 * ((x_plus - x_obs) ** 2).sum()
    L_minus = 0.5 * ((x_minus - x_obs) ** 2).sum()
    g_floor_numerical = (L_plus - L_minus) / (2 * eps)

    print(f"  Gradient check at floor_y={bad_floor:.3f} ({n} particles, K={K} steps):")
    print(f"    ∂L/∂floor_y kernel    = {g_floor_kernel:+.6e}")
    print(f"    ∂L/∂floor_y numerical = {g_floor_numerical:+.6e}")
    rel = abs(g_floor_kernel - g_floor_numerical) / max(abs(g_floor_numerical), 1e-9)
    print(f"    relative error        = {rel:.3e}")
    grad_pass = rel < 5e-2
    print(f"    {'[PASS]' if grad_pass else '[FAIL]'} multi-particle gradient agrees with finite-diff")
    print()

    # ── M7.E scale-up: SGD on floor_y ──
    print(f"  SGD: learning floor_y from observed multi-particle state")
    print(f"       (start floor_y={bad_floor:.3f}, target {true_floor:.3f})")
    learning_rate = 0.5 / n  # scale by n since gradient sums across particles
    cur_floor = float(bad_floor)
    losses = []
    for it in range(40):
        x_cur, _, tp, tc = fwd(x0, v0, K, dt, g_y, cur_floor)
        L = 0.5 * ((x_cur - x_obs) ** 2).sum()
        grad_x_final = (x_cur - x_obs).astype(np.float32)
        _, _, _, g_f = bwd(grad_x_final, tp, tc, dt)
        cur_floor -= learning_rate * g_f
        losses.append(L)
        if it % 5 == 0 or it == 39:
            print(f"    iter {it:3d}: floor_y={cur_floor:+.4f}  loss={L:.4e}  grad={g_f:+.3e}")

    err = abs(cur_floor - true_floor)
    print(f"    learned floor_y = {cur_floor:.4f}, true = {true_floor:.4f}, error = {err:.4e}")
    sgd_pass = err < 1e-2
    print(f"    {'[PASS]' if sgd_pass else '[FAIL]'} multi-particle SGD converged")
    print()

    if grad_pass and sgd_pass:
        print("  [OVERALL PASS] Multi-particle gradient learning works on hand-written Metal")
    else:
        print("  [OVERALL FAIL] see failures above")
    return 0 if (grad_pass and sgd_pass) else 1


if __name__ == "__main__":
    raise SystemExit(main())
