"""M9.A — validate density value backward chain (∂density/∂active_pos).

Forward: active positions → r²(active, static) → W(r²) → density(W)
Backward chain: ∂L/∂density → ∂L/∂W → ∂L/∂r² → ∂L/∂active

Test: small scene with active liquid + static boundary. Pick a loss
L = sum(w · density). Compute ∂L/∂active via the backward chain on
Metal, compare to finite-difference numerical gradient.

This is the foundation for M9 — once this passes, we can chain it
into the XPBD density constraint backward (next milestone).
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def fwd(active, static_p, h, mass):
    n_a = len(active); n_s = len(static_p)
    p_a = "/tmp/dens_a.bin"; p_s = "/tmp/dens_s.bin"
    p_d = "/tmp/dens_out.bin"; p_r2 = "/tmp/dens_r2.bin"
    active.astype(np.float32).tofile(p_a)
    static_p.astype(np.float32).tofile(p_s)
    subprocess.run(
        [BINARY, "density_as_fwd",
         str(n_a), str(n_s), str(h), str(mass), p_a, p_s, p_d, p_r2],
        check=True,
    )
    density = np.fromfile(p_d, dtype=np.float32)
    r2 = np.fromfile(p_r2, dtype=np.float32).reshape(n_a, n_s)
    return density, r2


def bwd(active, static_p, r2_saved, grad_density, h, mass):
    n_a = len(active); n_s = len(static_p)
    p_a = "/tmp/dens_a.bin"; p_s = "/tmp/dens_s.bin"
    p_r2 = "/tmp/dens_r2.bin"; p_gd = "/tmp/dens_gd.bin"
    active.astype(np.float32).tofile(p_a)
    static_p.astype(np.float32).tofile(p_s)
    r2_saved.astype(np.float32).tofile(p_r2)
    grad_density.astype(np.float32).tofile(p_gd)
    subprocess.run(
        [BINARY, "density_as_bwd",
         str(n_a), str(n_s), str(h), str(mass),
         p_a, p_s, p_r2, p_gd, "/tmp/dens_grad_a.bin"],
        check=True,
    )
    grad_active = np.fromfile("/tmp/density_as_grad_active.bin",
                              dtype=np.float32).reshape(n_a, 3)
    return grad_active


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    np.random.seed(7)
    n_a, n_s = 6, 30
    h = 1.5
    mass = 1.0

    # Active particles clustered, static spread around them so distances
    # vary across the kernel-active range.
    active = np.random.randn(n_a, 3).astype(np.float32) * 0.5
    static_p = np.random.randn(n_s, 3).astype(np.float32) * 1.5

    print(f"=== M9.A: density backward chain ===")
    print(f"  n_active={n_a}, n_static={n_s}, h={h}, mass={mass}")

    # Forward
    density, r2 = fwd(active, static_p, h, mass)
    print(f"  Density: min={density.min():.4f}  mean={density.mean():.4f}  "
          f"max={density.max():.4f}")

    # Random loss weights
    w = np.random.randn(n_a).astype(np.float32)
    L = float((w * density).sum())
    print(f"  L = w·density = {L:.4e}")
    print()

    # ── Kernel backward ──
    grad_active_kernel = bwd(active, static_p, r2, w, h, mass)

    # ── Numerical gradient (finite-difference per active particle, per axis) ──
    eps = 1e-3
    grad_active_numerical = np.zeros_like(active)
    for i in range(n_a):
        for ax in range(3):
            a_p = active.copy(); a_p[i, ax] += eps
            a_m = active.copy(); a_m[i, ax] -= eps
            d_p, _ = fwd(a_p, static_p, h, mass)
            d_m, _ = fwd(a_m, static_p, h, mass)
            L_p = float((w * d_p).sum())
            L_m = float((w * d_m).sum())
            grad_active_numerical[i, ax] = (L_p - L_m) / (2 * eps)

    # Compare
    err_abs = np.abs(grad_active_kernel - grad_active_numerical)
    err_max = float(err_abs.max())
    grad_mag_mean = float(np.linalg.norm(grad_active_numerical, axis=1).mean())
    err_rel = err_max / max(grad_mag_mean, 1e-9)

    print(f"  Gradient check (one element each row):")
    print(f"    {'particle':>8} {'axis':>4}   {'kernel':>14}   {'numerical':>14}   {'|abs|':>12}")
    for i in range(min(n_a, 4)):
        for ax in range(3):
            print(f"    {i:>8} {ax:>4}   "
                  f"{grad_active_kernel[i, ax]:>+14.6e}   "
                  f"{grad_active_numerical[i, ax]:>+14.6e}   "
                  f"{err_abs[i, ax]:>12.3e}")
    print()
    print(f"  Max abs error:           {err_max:.3e}")
    print(f"  Mean |∇| (numerical):    {grad_mag_mean:.3e}")
    print(f"  Max-error / mean-grad:   {err_rel:.3e}")

    passed = err_rel < 5e-2
    print(f"  {'[PASS]' if passed else '[FAIL]'} M9.A: density backward agrees with finite-diff")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
