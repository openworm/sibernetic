"""M9.C — validate density_constraint_grad_backward.

Forward (M6.4) at active particle i:
  grad_C[i]      = (mass/ρ_rest) · Σ_k ∇W_spiky(p_i - p_k)
  denom_helper[i] = Σ_k |∇W_spiky(p_i - p_k)|²

Where the sum spans active neighbors j ≠ i and static neighbors. Both
outputs depend on positions through ∇W_spiky's chain rule.

Test setup: small scene with active + static. Pick random ω = ∂L/∂grad_C
and ψ = ∂L/∂denom_helper. Define L = ω·grad_C + ψ·denom_helper. Compute
∂L/∂active via the new backward kernel and compare to finite-difference.
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def _write(path, arr):
    arr.astype(np.float32).tofile(path)


def fwd_chain(active, static_p, h, mass, rho_rest):
    """Compute grad_C and denom_helper from positions via existing ops."""
    n_a = len(active); n_s = len(static_p)
    p_a = "/tmp/cgb_a.bin"; p_s = "/tmp/cgb_s.bin"
    p_r2aa = "/tmp/cgb_r2aa.bin"; p_r2as = "/tmp/cgb_r2as.bin"
    p_gC = "/tmp/cgb_gC.bin"; p_dh = "/tmp/cgb_dh.bin"
    _write(p_a, active); _write(p_s, static_p)
    # 1) dist_active_active
    subprocess.run([BINARY, "dist_active_active",
                    str(n_a), p_a, p_r2aa], check=True)
    # 2) dist_active_static
    subprocess.run([BINARY, "dist_active_static",
                    str(n_a), str(n_s), p_a, p_s, p_r2as], check=True)
    # 3) density_constraint_grad
    subprocess.run([BINARY, "density_constraint_grad",
                    str(n_a), str(n_s), str(h), str(mass), str(rho_rest),
                    p_a, p_s, p_r2aa, p_r2as, p_gC, p_dh],
                   check=True)
    grad_C = np.fromfile(p_gC, dtype=np.float32).reshape(n_a, 3)
    denom_h = np.fromfile(p_dh, dtype=np.float32)
    return grad_C, denom_h


def bwd(active, static_p, h, mass, rho_rest, omega, psi):
    """Run the M9.C backward kernel; returns ∂L/∂active."""
    n_a = len(active); n_s = len(static_p)
    p_a = "/tmp/cgb_a.bin"; p_s = "/tmp/cgb_s.bin"
    p_r2aa = "/tmp/cgb_r2aa.bin"; p_r2as = "/tmp/cgb_r2as.bin"
    p_gomega = "/tmp/cgb_gomega.bin"; p_gpsi = "/tmp/cgb_gpsi.bin"
    p_out = "/tmp/cgb_out.bin"
    _write(p_a, active); _write(p_s, static_p)
    # Re-do the distance computations (cheap; we need the r² files
    # passed to backward kernel and they may have been overwritten).
    subprocess.run([BINARY, "dist_active_active",
                    str(n_a), p_a, p_r2aa], check=True)
    subprocess.run([BINARY, "dist_active_static",
                    str(n_a), str(n_s), p_a, p_s, p_r2as], check=True)
    _write(p_gomega, omega); _write(p_gpsi, psi)
    subprocess.run([BINARY, "density_constraint_grad_bwd",
                    str(n_a), str(n_s), str(h), str(mass), str(rho_rest),
                    p_a, p_s, p_r2aa, p_r2as,
                    p_gomega, p_gpsi, p_out],
                   check=True)
    return np.fromfile(p_out, dtype=np.float32).reshape(n_a, 3)


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    np.random.seed(13)
    n_a, n_s = 5, 20
    h = 1.5
    mass = 1.0
    rho_rest = 1.5

    active = np.random.randn(n_a, 3).astype(np.float32) * 0.6
    static_p = np.random.randn(n_s, 3).astype(np.float32) * 1.2

    print(f"=== M9.C: density_constraint_grad backward ===")
    print(f"  n_active={n_a}, n_static={n_s}, h={h}, ρ_rest={rho_rest}")

    grad_C, denom_h = fwd_chain(active, static_p, h, mass, rho_rest)
    print(f"  Forward outputs:")
    print(f"    |grad_C| mean={np.linalg.norm(grad_C, axis=1).mean():.4e}")
    print(f"    denom_helper mean={denom_h.mean():.4e}")

    omega = np.random.randn(n_a, 3).astype(np.float32)
    psi = np.random.randn(n_a).astype(np.float32)

    grad_kernel = bwd(active, static_p, h, mass, rho_rest, omega, psi)

    # Numerical: ∂L/∂active[i, ax] via central diff
    eps = 1e-3
    grad_num = np.zeros_like(active)
    for i in range(n_a):
        for ax in range(3):
            a_p = active.copy(); a_p[i, ax] += eps
            a_m = active.copy(); a_m[i, ax] -= eps
            gC_p, dh_p = fwd_chain(a_p, static_p, h, mass, rho_rest)
            gC_m, dh_m = fwd_chain(a_m, static_p, h, mass, rho_rest)
            L_p = (omega * gC_p).sum() + (psi * dh_p).sum()
            L_m = (omega * gC_m).sum() + (psi * dh_m).sum()
            grad_num[i, ax] = (L_p - L_m) / (2 * eps)

    err_abs = np.abs(grad_kernel - grad_num)
    err_max = float(err_abs.max())
    grad_mean = float(np.linalg.norm(grad_num, axis=1).mean())
    err_rel = err_max / max(grad_mean, 1e-9)

    print(f"\n  Sample (first 3 particles, all axes):")
    print(f"    {'p':>3} {'ax':>3}  {'kernel':>14}  {'numerical':>14}  {'abs err':>10}")
    for i in range(min(n_a, 3)):
        for ax in range(3):
            print(f"    {i:>3} {ax:>3}  {grad_kernel[i,ax]:>+14.6e}  "
                  f"{grad_num[i,ax]:>+14.6e}  {err_abs[i,ax]:>10.3e}")
    print()
    print(f"  Max abs error:  {err_max:.3e}")
    print(f"  Mean |∇|:       {grad_mean:.3e}")
    print(f"  Relative err:   {err_rel:.3e}")
    passed = err_rel < 5e-2
    print(f"  {'[PASS]' if passed else '[FAIL]'} M9.C: density_constraint_grad backward")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
