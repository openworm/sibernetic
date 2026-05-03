"""M9.D — validate solve_density_constraint_backward.

Tests every gradient output (∂L/∂{density, grad_C, denom_helper,
lambda_pre, rho_rest, pos_pre}) against finite-difference on the
forward kernel.

Forward (per active i, when C > 0):
  C = density/ρ_rest - 1
  D = |grad_C|²/m + (m/ρ²)·denom_helper + α/dt²
  Δλ = -(C + α/dt²·λ_pre) / D
  pos_post = pos_pre + grad_C·Δλ/m
  λ_post = λ_pre + Δλ

Loss: L = Σ (ω_i · pos_post[i]) + Σ (ψ_i · λ_post[i])
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def _w(p, a):
    a.astype(np.float32).tofile(p)


def fwd(pos_pre, lam_pre, density, grad_C, denom_h,
        rho_rest, mass, alpha_inv_dt2):
    n = len(pos_pre)
    p_pos = "/tmp/sdc_pos.bin"; p_lam = "/tmp/sdc_lam.bin"
    p_d = "/tmp/sdc_d.bin"; p_gc = "/tmp/sdc_gc.bin"; p_dh = "/tmp/sdc_dh.bin"
    p_pos_o = "/tmp/sdc_pos_o.bin"; p_lam_o = "/tmp/sdc_lam_o.bin"
    _w(p_pos, pos_pre); _w(p_lam, lam_pre)
    _w(p_d, density); _w(p_gc, grad_C); _w(p_dh, denom_h)
    subprocess.run(
        [BINARY, "solve_density_constraint_fwd",
         str(n), str(rho_rest), str(mass), str(alpha_inv_dt2),
         p_pos, p_lam, p_d, p_gc, p_dh, p_pos_o, p_lam_o],
        check=True,
    )
    pos_post = np.fromfile(p_pos_o, dtype=np.float32).reshape(n, 3)
    lam_post = np.fromfile(p_lam_o, dtype=np.float32)
    return pos_post, lam_post


def bwd(density, grad_C, denom_h, lam_pre, grad_pos_post, grad_lam_post,
        rho_rest, mass, alpha_inv_dt2):
    n = len(density)
    p_d = "/tmp/sdc_d.bin"; p_gc = "/tmp/sdc_gc.bin"; p_dh = "/tmp/sdc_dh.bin"
    p_lam = "/tmp/sdc_lam.bin"
    p_gP = "/tmp/sdc_gP.bin"; p_gL = "/tmp/sdc_gL.bin"
    out_pos = "/tmp/sdc_o_pos.bin"; out_lam = "/tmp/sdc_o_lam.bin"
    out_d = "/tmp/sdc_o_d.bin"; out_gc = "/tmp/sdc_o_gc.bin"
    out_dh = "/tmp/sdc_o_dh.bin"; out_r = "/tmp/sdc_o_r.bin"
    _w(p_d, density); _w(p_gc, grad_C); _w(p_dh, denom_h); _w(p_lam, lam_pre)
    _w(p_gP, grad_pos_post); _w(p_gL, grad_lam_post)
    subprocess.run(
        [BINARY, "solve_density_constraint_bwd",
         str(n), str(rho_rest), str(mass), str(alpha_inv_dt2),
         p_d, p_gc, p_dh, p_lam, p_gP, p_gL,
         out_pos, out_lam, out_d, out_gc, out_dh, out_r],
        check=True,
    )
    g_pos = np.fromfile(out_pos, dtype=np.float32).reshape(n, 3)
    g_lam = np.fromfile(out_lam, dtype=np.float32)
    g_d = np.fromfile(out_d, dtype=np.float32)
    g_gc = np.fromfile(out_gc, dtype=np.float32).reshape(n, 3)
    g_dh = np.fromfile(out_dh, dtype=np.float32)
    g_r_per = np.fromfile(out_r, dtype=np.float32)
    return g_pos, g_lam, g_d, g_gc, g_dh, g_r_per.sum()


def loss(omega, psi, pos, lam):
    return float((omega * pos).sum() + (psi * lam).sum())


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    np.random.seed(42)
    n = 6
    rho_rest = 1.5
    mass = 1.0
    alpha_inv_dt2 = 0.5
    eps = 1e-3

    # Set up state where C > 0 for some particles, C ≤ 0 for others
    # so we exercise both branches.
    pos_pre = np.random.randn(n, 3).astype(np.float32) * 0.5
    lam_pre = np.random.randn(n).astype(np.float32) * 0.1
    density = np.random.uniform(0.5, 3.0, n).astype(np.float32)  # mix of C+ and C-
    grad_C = np.random.randn(n, 3).astype(np.float32) * 0.5
    denom_h = np.random.uniform(0.1, 1.0, n).astype(np.float32)

    print(f"=== M9.D: solve_density_constraint backward ===")
    print(f"  n={n}, ρ_rest={rho_rest}, mass={mass}, α/dt²={alpha_inv_dt2}")
    print(f"  C values (density/ρ - 1): {density / rho_rest - 1.0}")
    print(f"  → {(density / rho_rest - 1.0 > 0).sum()} particles have constraint firing")

    omega = np.random.randn(n, 3).astype(np.float32)
    psi = np.random.randn(n).astype(np.float32)

    # Kernel backward
    pos_post0, lam_post0 = fwd(pos_pre, lam_pre, density, grad_C, denom_h,
                                rho_rest, mass, alpha_inv_dt2)
    L0 = loss(omega, psi, pos_post0, lam_post0)
    print(f"  Forward L = {L0:.4e}")

    g_pos, g_lam, g_d, g_gc, g_dh, g_r = bwd(
        density, grad_C, denom_h, lam_pre, omega, psi,
        rho_rest, mass, alpha_inv_dt2,
    )

    def fd(perturb_fn):
        L_p = loss(omega, psi, *fwd(*perturb_fn(+eps),
                                     rho_rest, mass, alpha_inv_dt2))
        L_m = loss(omega, psi, *fwd(*perturb_fn(-eps),
                                     rho_rest, mass, alpha_inv_dt2))
        return (L_p - L_m) / (2 * eps)

    # Helper to perturb one input
    def make_perturb(target, i, ax=None):
        def _p(e):
            args = (pos_pre.copy(), lam_pre.copy(), density.copy(),
                    grad_C.copy(), denom_h.copy())
            idx = {'pos': 0, 'lam': 1, 'd': 2, 'gc': 3, 'dh': 4}[target]
            arr = args[idx]
            if ax is not None:
                arr[i, ax] += e
            else:
                arr[i] += e
            return args
        return _p

    # Sample one element of each input gradient for the comparison.
    print()
    print(f"  {'input':>10}  {'i':>3} {'ax':>4}  {'kernel':>14}  {'numerical':>14}  {'rel err':>10}")
    fail = False
    # Pick one firing (C > 0) and one non-firing (C ≤ 0) particle
    # so both kernel branches get exercised.
    C_vals = density / rho_rest - 1.0
    firing_idx = int(np.where(C_vals > 0)[0][0])
    nonfiring_idx = int(np.where(C_vals <= 0)[0][0])
    test_particles = [firing_idx, nonfiring_idx]
    print(f"  Testing firing particle {firing_idx} (C={C_vals[firing_idx]:+.3f}) "
          f"and non-firing {nonfiring_idx} (C={C_vals[nonfiring_idx]:+.3f})")
    print()

    for target, name in [('pos', '∂L/∂pos'), ('lam', '∂L/∂λ'),
                          ('d', '∂L/∂dens'), ('gc', '∂L/∂grad_C'),
                          ('dh', '∂L/∂denom_h')]:
        for i in test_particles:
            if target in ('pos', 'gc'):
                axes = [0, 1, 2]
            else:
                axes = [None]
            for ax in axes:
                k_grad = (g_pos if target == 'pos'
                          else g_lam if target == 'lam'
                          else g_d if target == 'd'
                          else g_gc if target == 'gc'
                          else g_dh)
                if ax is not None:
                    k_val = k_grad[i, ax]
                else:
                    k_val = k_grad[i]
                num_val = fd(make_perturb(target, i, ax))
                rel = abs(k_val - num_val) / max(abs(num_val), 1e-9)
                axstr = str(ax) if ax is not None else '-'
                ok = rel < 5e-2
                if not ok and abs(num_val) > 1e-6:
                    fail = True
                marker = ' ' if ok or abs(num_val) < 1e-6 else '*'
                print(f"  {name:>10}  {i:>3} {axstr:>4}  "
                      f"{k_val:>+14.6e}  {num_val:>+14.6e}  {rel:>10.3e}{marker}")

    # ρ_rest gradient (scalar, sum across particles)
    print()
    L_p = loss(omega, psi, *fwd(pos_pre, lam_pre, density, grad_C, denom_h,
                                 rho_rest + eps, mass, alpha_inv_dt2))
    L_m = loss(omega, psi, *fwd(pos_pre, lam_pre, density, grad_C, denom_h,
                                 rho_rest - eps, mass, alpha_inv_dt2))
    g_r_num = (L_p - L_m) / (2 * eps)
    rel_r = abs(g_r - g_r_num) / max(abs(g_r_num), 1e-9)
    print(f"  ∂L/∂ρ_rest:  kernel={g_r:+.6e}  num={g_r_num:+.6e}  rel err={rel_r:.3e}")

    overall = not fail and rel_r < 5e-2
    print()
    print(f"  {'[PASS]' if overall else '[FAIL]'} M9.D: solve_density_constraint backward")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
