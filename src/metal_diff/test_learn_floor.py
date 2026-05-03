"""M7.D + M7.E demo: numerical gradient check + SGD on floor_y.

Drops a particle from height under gravity. After K steps it lands on
the floor. The loss is squared distance from a target final position
(below the actual floor → encourages floor_y to drop).

Tests:
  1. Numerical gradient validation: compare ∂L/∂floor_y from kernel
     backward to finite-difference. PASS if they agree.
  2. SGD parameter learning: start with wrong floor_y, run SGD using
     kernel gradients, verify floor_y converges to the value that
     minimizes the loss.

Run:
    .venv/bin/python src/metal_diff/test_learn_floor.py
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def fwd(x0, v0, K, dt, g_y, fl_y):
    """Run K-step forward. Returns (x_final, v_final, traj_pos, traj_clamp)."""
    n = len(x0)
    p_x = "/tmp/lf_x0.bin"; p_v = "/tmp/lf_v0.bin"
    p_tp = "/tmp/lf_tpos.bin"; p_tc = "/tmp/lf_tclmp.bin"
    p_vf = "/tmp/lf_vfinal.bin"
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
    """Run backward. Returns (grad_x_init, grad_v_init, grad_g_y, grad_floor_y)."""
    K = traj_clamp.shape[0]
    n = traj_clamp.shape[1]
    p_tp = "/tmp/lf_tpos.bin"; p_tc = "/tmp/lf_tclmp.bin"
    p_gxf = "/tmp/lf_gxfin.bin"
    p_gxi = "/tmp/lf_gxinit.bin"; p_gvi = "/tmp/lf_gvinit.bin"
    p_ggy = "/tmp/lf_ggy.bin"; p_gfy = "/tmp/lf_gfy.bin"
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


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    # ── Setup: 1 particle dropped from y=5, K steps long enough to ──
    # ── settle on the floor. With floor_y=true_floor, particle ends ──
    # ── at exactly true_floor (clamped). SGD should recover that.   ──
    n = 1
    K = 500          # 500 steps of dt=0.005 = 2.5s — long after impact
    dt = 0.005
    g_y = -9.8
    true_floor = 0.3  # clearly above 0 so particle is clamped from the start

    x0 = np.array([[0.0, 5.0, 0.0]], dtype=np.float32)
    v0 = np.zeros((n, 3), dtype=np.float32)

    # First, run forward with the TRUE floor and capture the "observed" final position.
    x_obs, _, _, _ = fwd(x0, v0, K, dt, g_y, true_floor)
    print(f"=== M7.D/E: learning floor_y by SGD ===")
    print(f"  Setup: 1 particle dropped from y=5, K={K} steps × dt={dt} = 2.5s sim")
    print(f"  True floor_y = {true_floor:.3f}, observed final y = {x_obs[0,1]:.4f}")
    print()

    # ── M7.D: numerical gradient validation ──
    # Use a wrong floor_y. After 2.5s settling, observed_final ≈ true_floor.
    bad_floor = 0.7
    x_bad, _, traj_pos, traj_clamp = fwd(x0, v0, K, dt, g_y, bad_floor)
    # Loss = 0.5 * sum((x_final - x_obs)²)  → ∂L/∂x_final = (x_final - x_obs)
    grad_x_final = (x_bad - x_obs).astype(np.float32)
    g_x, g_v, g_g_y, g_floor = bwd(grad_x_final, traj_pos, traj_clamp, dt)

    # Numerical gradient on floor_y via finite-diff
    eps = 1e-3
    x_plus, _, _, _ = fwd(x0, v0, K, dt, g_y, bad_floor + eps)
    x_minus, _, _, _ = fwd(x0, v0, K, dt, g_y, bad_floor - eps)
    L_plus = 0.5 * ((x_plus - x_obs) ** 2).sum()
    L_minus = 0.5 * ((x_minus - x_obs) ** 2).sum()
    g_floor_numerical = (L_plus - L_minus) / (2 * eps)

    print(f"  M7.D — gradient check at floor_y={bad_floor:.3f}")
    print(f"    L = 0.5·||x_final - x_obs||²  =  {0.5*((x_bad-x_obs)**2).sum():.6e}")
    print(f"    ∂L/∂floor_y kernel    = {g_floor:+.6e}")
    print(f"    ∂L/∂floor_y numerical = {g_floor_numerical:+.6e}")
    rel_err = abs(g_floor - g_floor_numerical) / max(abs(g_floor_numerical), 1e-9)
    print(f"    relative error        = {rel_err:.3e}")
    grad_pass = rel_err < 5e-2
    print(f"    {'[PASS]' if grad_pass else '[FAIL]'} M7.D: gradient agrees with finite-diff")
    print()

    # ── M7.E: SGD on floor_y ──
    print(f"  M7.E — SGD on floor_y (start at {bad_floor:.3f}, target {true_floor:.3f})")
    learning_rate = 1.0
    cur_floor = float(bad_floor)
    for it in range(20):
        x_cur, _, tp, tc = fwd(x0, v0, K, dt, g_y, cur_floor)
        L = 0.5 * ((x_cur - x_obs) ** 2).sum()
        grad_x_final = (x_cur - x_obs).astype(np.float32)
        _, _, _, g_f = bwd(grad_x_final, tp, tc, dt)
        cur_floor -= learning_rate * g_f
        if it % 4 == 0 or it == 19:
            print(f"    iter {it:2d}: floor_y={cur_floor:+.4f}  loss={L:.4e}  grad={g_f:+.3e}")

    learned_err = abs(cur_floor - true_floor)
    print(f"    learned floor_y = {cur_floor:.4f}, true = {true_floor:.4f}, error = {learned_err:.4e}")
    sgd_pass = learned_err < 1e-2
    print(f"    {'[PASS]' if sgd_pass else '[FAIL]'} M7.E: SGD converged to true floor_y")
    print()

    all_pass = grad_pass and sgd_pass
    if all_pass:
        print("  [OVERALL PASS] M7.D+E: end-to-end gradient learning works on hand-written Metal")
    else:
        print("  [OVERALL FAIL] see failures above")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
