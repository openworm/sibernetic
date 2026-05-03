"""M9.E — end-to-end SGD on rho_rest using the M9 backward chain.

The capstone demo: given an observed final particle configuration
(generated with a TRUE rho_rest), use SGD to recover that rho_rest by
backpropagating through the XPBD density constraint solve.

Pipeline (one constraint-solve iteration, no time stepping):
  Forward:
    pos → r²(pos) → W(r²) → density(W)
        → grad_C(pos, ρ) [via density_constraint_grad]
        → denom_helper(pos)
        → solve_density_constraint(pos, λ=0, density, grad_C, denom_h, ρ)
        → pos_post

  Loss: L = ½ · ||pos_post − obs||²

  Backward (∂L/∂ρ chain):
    ∂L/∂pos_post = pos_post − obs
    solve_density_constraint_bwd → {∂L/∂grad_C, ∂L/∂ρ_partial}
    Add the IMPLICIT dependency: ∂L/∂ρ_via_grad_C = -(∂L/∂grad_C · grad_C)/ρ
      (since grad_C = (1/ρ) · Σ ..., so ∂grad_C/∂ρ = -grad_C/ρ)
    Total ∂L/∂ρ = ∂L/∂ρ_partial + ∂L/∂ρ_via_grad_C

We don't propagate further into pos here (that would need M9.A/B/C
chain wired through). For LEARNING ρ alone, the local gradient at
this step suffices.
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def _w(p, a):
    a.astype(np.float32).tofile(p)


def fwd_density_chain(active, static_p, h, mass, rho_rest):
    """Compute density, grad_C, denom_helper from positions."""
    n_a = len(active); n_s = len(static_p)
    p_a = "/tmp/lr_a.bin"; p_s = "/tmp/lr_s.bin"
    p_r2aa = "/tmp/lr_r2aa.bin"; p_r2as = "/tmp/lr_r2as.bin"
    p_d_aa = "/tmp/lr_d_aa.bin"; p_d_as = "/tmp/lr_d_as.bin"
    p_r2_aa_save = "/tmp/lr_r2aa_save.bin"; p_r2_as_save = "/tmp/lr_r2as_save.bin"
    p_gC = "/tmp/lr_gC.bin"; p_dh = "/tmp/lr_dh.bin"
    _w(p_a, active); _w(p_s, static_p)

    # density_aa: gives density from active-active SPH neighbors
    subprocess.run([BINARY, "density_aa_fwd",
                    str(n_a), str(h), str(mass),
                    p_a, p_d_aa, p_r2_aa_save], check=True)
    density_aa = np.fromfile(p_d_aa, dtype=np.float32)

    # density_as: gives density from active-static
    subprocess.run([BINARY, "density_as_fwd",
                    str(n_a), str(n_s), str(h), str(mass),
                    p_a, p_s, p_d_as, p_r2_as_save], check=True)
    density_as = np.fromfile(p_d_as, dtype=np.float32)

    density = density_aa + density_as

    # grad_C and denom_helper via M6.4
    subprocess.run([BINARY, "density_constraint_grad",
                    str(n_a), str(n_s), str(h), str(mass), str(rho_rest),
                    p_a, p_s, p_r2_aa_save, p_r2_as_save, p_gC, p_dh],
                   check=True)
    grad_C = np.fromfile(p_gC, dtype=np.float32).reshape(n_a, 3)
    denom_h = np.fromfile(p_dh, dtype=np.float32)

    return density, grad_C, denom_h, p_r2_aa_save, p_r2_as_save


def fwd_constraint(pos_pre, density, grad_C, denom_h, rho_rest, mass, alpha_inv_dt2):
    """Run solve_density_constraint."""
    n = len(pos_pre)
    p_pos = "/tmp/lr_pos.bin"; p_lam = "/tmp/lr_lam.bin"
    p_d = "/tmp/lr_d.bin"; p_gc = "/tmp/lr_gc.bin"; p_dh = "/tmp/lr_dh.bin"
    p_pos_o = "/tmp/lr_pos_o.bin"; p_lam_o = "/tmp/lr_lam_o.bin"
    _w(p_pos, pos_pre); _w(p_lam, np.zeros(n, dtype=np.float32))
    _w(p_d, density); _w(p_gc, grad_C); _w(p_dh, denom_h)
    subprocess.run([BINARY, "solve_density_constraint_fwd",
                    str(n), str(rho_rest), str(mass), str(alpha_inv_dt2),
                    p_pos, p_lam, p_d, p_gc, p_dh, p_pos_o, p_lam_o],
                   check=True)
    pos_post = np.fromfile(p_pos_o, dtype=np.float32).reshape(n, 3)
    return pos_post


def full_forward(active, static_p, h, mass, rho_rest, alpha_inv_dt2):
    density, grad_C, denom_h, _, _ = fwd_density_chain(
        active, static_p, h, mass, rho_rest)
    pos_post = fwd_constraint(active, density, grad_C, denom_h,
                              rho_rest, mass, alpha_inv_dt2)
    return pos_post, density, grad_C, denom_h


def grad_rho(active, static_p, h, mass, rho_rest, alpha_inv_dt2,
             grad_x_post):
    """Compute ∂L/∂ρ via M9.D + implicit grad_C dependency."""
    density, grad_C, denom_h, _, _ = fwd_density_chain(
        active, static_p, h, mass, rho_rest)
    n = len(active)
    p_d = "/tmp/lr_d.bin"; p_gc = "/tmp/lr_gc.bin"
    p_dh = "/tmp/lr_dh.bin"; p_lam = "/tmp/lr_lam.bin"
    p_gP = "/tmp/lr_gP.bin"; p_gL = "/tmp/lr_gL.bin"
    _w(p_d, density); _w(p_gc, grad_C); _w(p_dh, denom_h)
    _w(p_lam, np.zeros(n, dtype=np.float32))
    _w(p_gP, grad_x_post); _w(p_gL, np.zeros(n, dtype=np.float32))
    out_pos = "/tmp/lr_o_pos.bin"; out_lam = "/tmp/lr_o_lam.bin"
    out_d = "/tmp/lr_o_d.bin"; out_gc = "/tmp/lr_o_gc.bin"
    out_dh = "/tmp/lr_o_dh.bin"; out_r = "/tmp/lr_o_r.bin"
    subprocess.run([BINARY, "solve_density_constraint_bwd",
                    str(n), str(rho_rest), str(mass), str(alpha_inv_dt2),
                    p_d, p_gc, p_dh, p_lam, p_gP, p_gL,
                    out_pos, out_lam, out_d, out_gc, out_dh, out_r],
                   check=True)
    g_gc = np.fromfile(out_gc, dtype=np.float32).reshape(n, 3)
    g_r_per = np.fromfile(out_r, dtype=np.float32)
    g_r_partial = float(g_r_per.sum())
    # Implicit dependency: ∂L/∂ρ_via_grad_C = -(g_gc · grad_C) / ρ
    g_r_implicit = -float((g_gc * grad_C).sum()) / rho_rest
    return g_r_partial + g_r_implicit


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    np.random.seed(31)
    n_a, n_s = 4, 12
    h = 1.5
    mass = 1.0
    alpha_inv_dt2 = 0.5

    # Active particles cluster, static spread around
    active = np.random.randn(n_a, 3).astype(np.float32) * 0.4
    static_p = np.random.randn(n_s, 3).astype(np.float32) * 1.2

    # Generate observation at TRUE rho_rest
    true_rho = 1.5
    pos_obs, _, _, _ = full_forward(active, static_p, h, mass,
                                     true_rho, alpha_inv_dt2)

    print(f"=== M9.E: end-to-end SGD on rho_rest ===")
    print(f"  n_active={n_a}, n_static={n_s}, h={h}, α/dt²={alpha_inv_dt2}")
    print(f"  True ρ = {true_rho}")

    # Validate full ∂L/∂ρ vs finite-diff at a perturbed value
    bad_rho = 2.0
    pos_bad, dens_bad, gC_bad, dh_bad = full_forward(
        active, static_p, h, mass, bad_rho, alpha_inv_dt2)
    grad_x_post = (pos_bad - pos_obs).astype(np.float32)
    g_r_kernel = grad_rho(active, static_p, h, mass, bad_rho,
                           alpha_inv_dt2, grad_x_post)

    eps = 1e-3
    pos_p, _, _, _ = full_forward(active, static_p, h, mass,
                                   bad_rho + eps, alpha_inv_dt2)
    pos_m, _, _, _ = full_forward(active, static_p, h, mass,
                                   bad_rho - eps, alpha_inv_dt2)
    L_p = 0.5 * ((pos_p - pos_obs) ** 2).sum()
    L_m = 0.5 * ((pos_m - pos_obs) ** 2).sum()
    g_r_num = (L_p - L_m) / (2 * eps)

    print(f"  Gradient check at ρ={bad_rho:.3f}:")
    print(f"    ∂L/∂ρ kernel    = {g_r_kernel:+.6e}")
    print(f"    ∂L/∂ρ numerical = {g_r_num:+.6e}")
    rel = abs(g_r_kernel - g_r_num) / max(abs(g_r_num), 1e-9)
    print(f"    relative error  = {rel:.3e}")
    grad_pass = rel < 0.1
    print(f"    {'[PASS]' if grad_pass else '[FAIL]'} M9.E gradient")
    print()

    # SGD
    print(f"  SGD: learn ρ from observation (start {bad_rho}, target {true_rho})")
    learning_rate = 1.0
    cur = float(bad_rho)
    for it in range(40):
        pos_cur, _, _, _ = full_forward(active, static_p, h, mass,
                                         cur, alpha_inv_dt2)
        L = 0.5 * ((pos_cur - pos_obs) ** 2).sum()
        grad_x = (pos_cur - pos_obs).astype(np.float32)
        g = grad_rho(active, static_p, h, mass, cur, alpha_inv_dt2, grad_x)
        step = learning_rate * g
        # Cap step to ±20% of cur to avoid blowing up
        max_step = cur * 0.2
        step = max(min(step, max_step), -max_step)
        cur = max(cur - step, 0.01)
        if it % 5 == 0 or it == 39:
            print(f"    iter {it:3d}: ρ={cur:.4f}  L={L:.4e}  ∂L/∂ρ={g:+.3e}")

    pos_final, _, _, _ = full_forward(active, static_p, h, mass,
                                       cur, alpha_inv_dt2)
    L_final = 0.5 * ((pos_final - pos_obs) ** 2).sum()
    err = abs(cur - true_rho) / true_rho
    print(f"    learned ρ = {cur:.4f}  true = {true_rho}  rel err = {err:.3e}")
    print(f"    final L = {L_final:.3e}  (initial {0.5*((pos_bad-pos_obs)**2).sum():.3e})")
    sgd_pass = err < 0.05
    print(f"    {'[PASS]' if sgd_pass else '[FAIL]'} SGD recovered ρ to <5%")
    print()

    overall = grad_pass and sgd_pass
    if overall:
        print("[OVERALL PASS] M9.E: end-to-end SGD on rho_rest works")
    else:
        print("[OVERALL FAIL] see above")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
