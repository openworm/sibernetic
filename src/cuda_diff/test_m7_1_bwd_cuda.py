"""FD-validate solve_density_constraint_bwd (M7.1_bwd / M9.D).

Forward (per active i):
    C = density/ρ_rest - 1
    if C ≤ 0: skip (one-sided projection, identity)
    else:
      D = |g|²/m + (m/ρ²)·helper + α/dt²
      Δλ = -(C + (α/dt²)·λ_pre) / D
      pos_post = pos_pre + g·Δλ/m
      λ_post   = λ_pre + Δλ

Loss: L = sum(ω · pos_post) + sum(ψ · λ_post)

Validates the kernel's 7 output gradients (grad_pos_pre, grad_lambda_pre,
grad_density, grad_grad_C, grad_denom_h, grad_rho_rest, grad_alpha) against
centered FD on each input. Densities are kept clearly above OR clearly
below ρ_rest by margin > 2·ε so FD perturbations don't flip the C-sign
branch.

Per-particle outputs (grad_rho_rest, grad_alpha) are summed before
comparing to the scalar-FD gradient.

Local-only; no Metal needed.
"""
import os
import platform
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY_NAME = "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda"
BINARY = os.path.join(HERE, BINARY_NAME)
TMP = tempfile.gettempdir()


def _run(args):
    subprocess.run([BINARY] + [str(a) for a in args], check=True)


def fwd(pos_pre, density, grad_C, denom_h, lambda_pre,
        rho_rest, mass, alpha_inv_dt2):
    """Python reference forward. Returns (pos_post, lambda_post)."""
    n = len(pos_pre)
    pos_post = pos_pre.copy()
    lambda_post = lambda_pre.copy()
    for i in range(n):
        C = density[i] / rho_rest - 1.0
        if C <= 0.0:
            continue
        g = grad_C[i]
        g2 = float(g @ g)
        helper = denom_h[i]
        D = g2 / mass + (mass / (rho_rest ** 2)) * helper + alpha_inv_dt2
        if D < 1e-9:
            continue
        dlambda = -(C + alpha_inv_dt2 * lambda_pre[i]) / D
        pos_post[i] = pos_pre[i] + g * (dlambda / mass)
        lambda_post[i] = lambda_pre[i] + dlambda
    return pos_post, lambda_post


def cuda_bwd(pos_pre, density, grad_C, denom_h, lambda_pre,
             omega, psi, rho_rest, mass, alpha_inv_dt2):
    n = len(pos_pre)
    p_dens  = os.path.join(TMP, "m71_density.bin")
    p_gC    = os.path.join(TMP, "m71_gradC.bin")
    p_dh    = os.path.join(TMP, "m71_denomh.bin")
    p_lp    = os.path.join(TMP, "m71_lampre.bin")
    p_gpp   = os.path.join(TMP, "m71_gpp.bin")
    p_glp   = os.path.join(TMP, "m71_glp.bin")
    p_gpre  = os.path.join(TMP, "m71_gpre.bin")
    p_glpre = os.path.join(TMP, "m71_glpre.bin")
    p_gdens = os.path.join(TMP, "m71_gdens.bin")
    p_ggC   = os.path.join(TMP, "m71_ggC.bin")
    p_gdh   = os.path.join(TMP, "m71_gdh.bin")
    p_grho  = os.path.join(TMP, "m71_grho.bin")
    p_galp  = os.path.join(TMP, "m71_galp.bin")

    density.astype(np.float32).tofile(p_dens)
    grad_C.astype(np.float32).tofile(p_gC)
    denom_h.astype(np.float32).tofile(p_dh)
    lambda_pre.astype(np.float32).tofile(p_lp)
    omega.astype(np.float32).tofile(p_gpp)
    psi.astype(np.float32).tofile(p_glp)
    np.zeros((n, 3), np.float32).tofile(p_gpre)

    _run(["solve_density_constraint_bwd",
          n, rho_rest, mass, alpha_inv_dt2,
          p_dens, p_gC, p_dh, p_lp, p_gpp, p_glp,
          p_gpre, p_glpre, p_gdens, p_ggC, p_gdh, p_grho, p_galp])

    return {
        "pos_pre":   np.fromfile(p_gpre,  dtype=np.float32).reshape(n, 3),
        "lambda_pre": np.fromfile(p_glpre, dtype=np.float32),
        "density":    np.fromfile(p_gdens, dtype=np.float32),
        "grad_C":     np.fromfile(p_ggC,   dtype=np.float32).reshape(n, 3),
        "denom_h":    np.fromfile(p_gdh,   dtype=np.float32),
        "rho_rest":   float(np.fromfile(p_grho,  dtype=np.float32).sum()),
        "alpha":      float(np.fromfile(p_galp,  dtype=np.float32).sum()),
    }


def fd_vector(loss_at, x, eps):
    """Centered FD over each component of x; loss_at(x_new) returns scalar L."""
    g = np.zeros_like(x)
    flat_x = x.ravel()
    flat_g = g.ravel()
    for k in range(flat_x.size):
        save = flat_x[k]
        flat_x[k] = save + eps; Lp = loss_at(x)
        flat_x[k] = save - eps; Lm = loss_at(x)
        flat_x[k] = save
        flat_g[k] = (Lp - Lm) / (2 * eps)
    return g


def fd_scalar(loss_at, x_ref, set_x, eps):
    """Centered FD over a scalar parameter."""
    set_x(x_ref + eps); Lp = loss_at()
    set_x(x_ref - eps); Lm = loss_at()
    set_x(x_ref)
    return (Lp - Lm) / (2 * eps)


def test_m7_1_bwd(eps=1e-3, tol=2e-3):
    np.random.seed(42)
    n = 5
    rho_rest = 1.0
    mass = 1.0
    alpha_inv_dt2 = 0.7

    pos_pre = np.random.randn(n, 3).astype(np.float64) * 0.3
    grad_C  = np.random.randn(n, 3).astype(np.float64) * 0.4
    denom_h = np.abs(np.random.randn(n).astype(np.float64)) * 0.5 + 0.1
    lambda_pre = np.random.randn(n).astype(np.float64) * 0.2

    # Densities: half clearly above ρ_rest (C > 0 branch), half clearly below
    # (C ≤ 0 branch). Margin >> ε so FD doesn't flip the mask.
    density = np.array([1.5, 1.4, 0.5, 0.6, 1.3], dtype=np.float64)

    omega = np.random.randn(n, 3).astype(np.float64) * 0.5
    psi   = np.random.randn(n).astype(np.float64) * 0.3

    state = {
        "pos_pre": pos_pre.astype(np.float32),
        "density": density.astype(np.float32),
        "grad_C":  grad_C.astype(np.float32),
        "denom_h": denom_h.astype(np.float32),
        "lambda_pre": lambda_pre.astype(np.float32),
        "rho_rest": rho_rest,
        "mass": mass,
        "alpha_inv_dt2": alpha_inv_dt2,
    }

    cuda = cuda_bwd(state["pos_pre"], state["density"], state["grad_C"],
                    state["denom_h"], state["lambda_pre"],
                    omega.astype(np.float32), psi.astype(np.float32),
                    rho_rest, mass, alpha_inv_dt2)

    def loss():
        pp, lp = fwd(state["pos_pre"], state["density"], state["grad_C"],
                     state["denom_h"], state["lambda_pre"],
                     state["rho_rest"], state["mass"], state["alpha_inv_dt2"])
        return float((omega * pp).sum() + (psi * lp).sum())

    def loss_at_arr(_unused):
        return loss()

    print("=== M7.1_bwd solve_density_constraint_bwd ===")
    print(f"  n={n}, C>0 particles: {int((density > rho_rest).sum())}, "
          f"C<=0 particles: {int((density <= rho_rest).sum())}")

    results = []
    # pos_pre FD
    fd_pp = fd_vector(loss_at_arr, state["pos_pre"], eps)
    results.append(("pos_pre",     cuda["pos_pre"],     fd_pp))

    # lambda_pre FD
    fd_lp = fd_vector(loss_at_arr, state["lambda_pre"], eps)
    results.append(("lambda_pre",  cuda["lambda_pre"],  fd_lp))

    # density FD
    fd_d = fd_vector(loss_at_arr, state["density"], eps)
    results.append(("density",     cuda["density"],     fd_d))

    # grad_C FD
    fd_g = fd_vector(loss_at_arr, state["grad_C"], eps)
    results.append(("grad_C",      cuda["grad_C"],      fd_g))

    # denom_h FD
    fd_dh = fd_vector(loss_at_arr, state["denom_h"], eps)
    results.append(("denom_h",     cuda["denom_h"],     fd_dh))

    # rho_rest FD (scalar)
    def set_rho(v):
        state["rho_rest"] = float(v)
    fd_rho = fd_scalar(loss, rho_rest, set_rho, eps)
    results.append(("rho_rest",    cuda["rho_rest"],    fd_rho))

    # alpha_inv_dt2 FD (scalar)
    def set_alpha(v):
        state["alpha_inv_dt2"] = float(v)
    fd_alpha = fd_scalar(loss, alpha_inv_dt2, set_alpha, eps)
    results.append(("alpha",       cuda["alpha"],       fd_alpha))

    all_ok = True
    print()
    print(f"  {'gradient':<13} {'analytic':>14} {'FD':>14} {'rel_err':>10}")
    print("  " + "-" * 55)
    for name, an, fd in results:
        an_arr = np.atleast_1d(np.asarray(an, dtype=np.float64))
        fd_arr = np.atleast_1d(np.asarray(fd, dtype=np.float64))
        err = float(np.abs(an_arr - fd_arr).max())
        mag = float(np.abs(fd_arr).max()) + 1e-9
        rel = err / mag
        ok = rel < tol
        all_ok &= ok
        if an_arr.size == 1:
            print(f"  {name:<13} {float(an_arr[0]):>14.6f} "
                  f"{float(fd_arr[0]):>14.6f} {rel:>10.2e} "
                  f"{'OK' if ok else 'FAIL'}")
        else:
            # Print max-magnitude element comparison.
            k = int(np.abs(fd_arr).argmax())
            print(f"  {name:<13} {float(an_arr.ravel()[k]):>14.6f} "
                  f"{float(fd_arr.ravel()[k]):>14.6f} {rel:>10.2e} "
                  f"{'OK' if ok else 'FAIL'}")
    print()
    print(f"  {'[PASS]' if all_ok else '[FAIL]'}")
    print()
    return all_ok


def main():
    if not os.path.exists(BINARY):
        print(f"sib_cuda binary missing at {BINARY}")
        return 1
    return 0 if test_m7_1_bwd() else 1


if __name__ == "__main__":
    sys.exit(main())
