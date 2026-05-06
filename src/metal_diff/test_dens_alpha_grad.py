"""FD-validation of the new ∂L/∂alpha_inv_dt2 from solve_density_constraint_backward.

Forward:
  D = |g|²/m + (m/ρ²)·helper + α
  Δλ = -(C + α·λ_pre) / D
  pos_post = pos_pre + g·Δλ/m
  λ_post = λ_pre + Δλ

We compute analytic ∂L/∂α (per-particle, summed) for a simple loss and
compare to a centered finite-difference ground truth.
"""
import os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")


def run_fwd(n, rho_rest, mass, alpha, density, grad_C, denom_h, lambda_pre,
             pos_pre, n_iters=1):
    """Single-iteration forward (matches kernel logic)."""
    pos_post = pos_pre.copy()
    lambda_post = lambda_pre.copy()
    for i in range(n):
        C = density[i] / rho_rest - 1.0
        if C <= 0: continue
        g = grad_C[i]
        g2 = float(np.dot(g, g))
        h = denom_h[i]
        D = g2 / mass + (mass / (rho_rest * rho_rest)) * h + alpha
        if D < 1e-9: continue
        dlambda = -(C + alpha * lambda_pre[i]) / D
        pos_post[i] = pos_pre[i] + g * dlambda / mass
        lambda_post[i] = lambda_pre[i] + dlambda
    return pos_post, lambda_post


def loss_fn(pos_post, lambda_post, w_pos, w_lam):
    """L = w_pos · pos_post + w_lam · λ_post (linear so easy to FD)."""
    return float(np.sum(pos_post * w_pos) + np.sum(lambda_post * w_lam))


def main():
    np.random.seed(7)
    n = 4
    rho_rest = 1000.0
    mass = 2e-12
    alpha = 1e-3 / (2e-5 ** 2)   # alpha_dens / dt² for typical Sibernetic params

    density = np.random.uniform(rho_rest * 1.05, rho_rest * 1.30, size=n).astype(np.float32)
    grad_C = np.random.randn(n, 3).astype(np.float32) * 0.01
    denom_h = np.random.uniform(0.5, 2.0, size=n).astype(np.float32)
    lambda_pre = np.random.randn(n).astype(np.float32) * 0.1
    pos_pre = np.random.randn(n, 3).astype(np.float32)

    w_pos = np.random.randn(n, 3).astype(np.float32)
    w_lam = np.random.randn(n).astype(np.float32)

    # === Analytic gradient via kernel ===
    # Compute grad_pos_post = w_pos, grad_lambda_post = w_lam
    paths = {}
    def w(name, arr):
        p = f"/tmp/test_dens_alpha_{name}.bin"
        arr.astype(np.float32).tofile(p)
        paths[name] = p
        return p
    w('density', density); w('grad_C', grad_C.flatten())
    w('denom_h', denom_h); w('lambda_pre', lambda_pre)
    w('grad_pos_post', w_pos.flatten()); w('grad_lambda_post', w_lam)

    out = {}
    for k in ['gp_out', 'gl_out', 'gd_out', 'ggC_out', 'gh_out', 'gr_out', 'galpha_out']:
        out[k] = f"/tmp/test_dens_alpha_{k}.bin"

    cmd = [BIN, "solve_density_constraint_bwd",
           str(n), f"{rho_rest:.6e}", f"{mass:.6e}", f"{alpha:.6e}",
           paths['density'], paths['grad_C'], paths['denom_h'], paths['lambda_pre'],
           paths['grad_pos_post'], paths['grad_lambda_post'],
           out['gp_out'], out['gl_out'], out['gd_out'],
           out['ggC_out'], out['gh_out'], out['gr_out'], out['galpha_out']]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("BIN failed:", r.stderr); return 1

    grad_alpha_per = np.fromfile(out['galpha_out'], dtype=np.float32)
    analytic = float(grad_alpha_per.sum())

    # === Numerical gradient via centered FD ===
    eps = max(abs(alpha) * 1e-4, 1e-8)
    pos_p, lam_p = run_fwd(n, rho_rest, mass, alpha + eps,
                            density, grad_C, denom_h, lambda_pre, pos_pre)
    pos_m, lam_m = run_fwd(n, rho_rest, mass, alpha - eps,
                            density, grad_C, denom_h, lambda_pre, pos_pre)
    L_p = loss_fn(pos_p, lam_p, w_pos, w_lam)
    L_m = loss_fn(pos_m, lam_m, w_pos, w_lam)
    numerical = (L_p - L_m) / (2 * eps)

    rel_err = abs(analytic - numerical) / max(abs(numerical), 1e-9)
    print(f"\n  ∂L/∂alpha_inv_dt2 (per α=alpha_inv_dt2={alpha:.3e}):")
    print(f"    analytic (per-particle sum): {analytic:+.6e}")
    print(f"    numerical (centered FD):     {numerical:+.6e}")
    print(f"    rel err: {rel_err:.3e}")

    print(f"\n  per-particle analytic: {grad_alpha_per}")

    if rel_err < 1e-3:
        print(f"\n  [PASS] ∂L/∂alpha_inv_dt2 within tolerance (1e-3)")
        return 0
    else:
        print(f"\n  [FAIL] ∂L/∂alpha_inv_dt2 rel err {rel_err:.3e} > 1e-3")
        return 1


if __name__ == "__main__":
    sys.exit(main())
