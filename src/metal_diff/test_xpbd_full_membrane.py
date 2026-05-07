"""End-to-end FD validation of xpbd_full_bwd through M10 membranes.

Synthetic 4-particle scene: 3 elastic verts forming one membrane
triangle + 1 liquid above. Run K small XPBD steps with membranes ON,
compute scalar loss on final liquid y, validate ∂L/∂x_init from
xpbd_full_bwd against centered finite differences on each particle.

PASS threshold: rel_err < 5e-2 per particle (looser than the
single-kernel tests because the multi-step chain accumulates fp32
roundoff and chaos).
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "sib_metal")
TMP = "/tmp"


def write_floats(path, arr):
    np.asarray(arr, dtype=np.float32).tofile(path)


def write_ints(path, arr):
    np.asarray(arr, dtype=np.int32).tofile(path)


def setup():
    """Write scenario binaries to /tmp/xfm_*. Return params dict."""
    n_active = 4
    n_elastic = 3
    n_membranes = 1
    n_static = 1   # one boundary particle far away (kept for signature compat)
    r0 = 2.0

    # Triangle in xz-plane at y=0; liquid hovers ~0.6 above centroid.
    pos = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.05, 0.0],
        [0.5, 0.0, 1.0],
        [0.5, 0.7, 0.3],
    ], dtype=np.float32)
    vel = np.zeros_like(pos)
    pos_static = np.array([[100.0, 100.0, 100.0]], dtype=np.float32)
    membranes = np.array([[0, 1, 2]], dtype=np.int32)
    pmem_idx = np.full((n_elastic, 7), -1, dtype=np.int32)
    for v in range(n_elastic):
        pmem_idx[v, 0] = 0

    write_floats(f"{TMP}/xfm_pos.bin", pos)
    write_floats(f"{TMP}/xfm_vel.bin", vel)
    write_floats(f"{TMP}/xfm_pos_static.bin", pos_static)
    write_ints  (f"{TMP}/xfm_mem.bin", membranes)
    write_ints  (f"{TMP}/xfm_pmi.bin", pmem_idx)

    return dict(n_active=n_active, n_elastic=n_elastic,
                n_static=n_static, n_membranes=n_membranes,
                r0=r0, pos=pos, vel=vel)


def run_fwd(p, K, h=2.0, mass=1.0, rho=0.3, dt=0.001, gravity=-2.0,
            alpha_dens=1e-3, sim_scale=1.0, use_floor=False, floor_y=-10.0):
    """Run xpbd_full_fwd with membranes; return final pos [n_active, 3]."""
    state_path = f"{TMP}/xfm_state.bin"
    cmd = [BIN, "xpbd_full_fwd",
           str(p['n_active']), str(p['n_static']), str(K),
           str(h), str(mass), f"{rho:.15e}",
           str(dt), str(gravity), f"{alpha_dens:.15e}",
           f"{TMP}/xfm_pos.bin", f"{TMP}/xfm_vel.bin",
           f"{TMP}/xfm_pos_static.bin",
           state_path,
           f"{sim_scale:.15e}",
           "0.0",   # visc_pair_coef = 0
           "0.0",   # spring_K = 0
           "/tmp/xfm_dummy_bonds.bin"]   # not used (springs off)
    # Floor args only if enabled (presence of arg toggles use_floor in C++).
    if use_floor:
        cmd.append(f"{floor_y:.6f}")
        cmd.append("0.0")   # restitution
    # Membrane args (positional after restitution; need restitution if floor off
    # to keep ordinal positions consistent — ops_xpbd_full.mm expects them
    # after argv[20]=restitution).
    if not use_floor:
        # Pad: floor_y + restitution slots so membrane args land at the right argv index.
        # But if use_floor=false, the C++ won't read argv[19/20]. So we shouldn't pad.
        # Actually we MUST set use_floor=true to use membranes (the ordinals require it).
        raise RuntimeError("test currently requires use_floor=True due to argv ordering")
    cmd.extend([
        str(p['n_membranes']), str(p['n_elastic']),
        f"{p['r0']:.15e}",
        f"{TMP}/xfm_mem.bin", f"{TMP}/xfm_pmi.bin"])
    # Create dummy bonds.bin since the path is required when use_floor.
    with open("/tmp/xfm_dummy_bonds.bin", "wb") as f:
        pass
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"fwd: {r.stderr}")

    # Decode state file: K * per_step + (K+1)*n*3 + n*3
    n = p['n_active']
    extra = (1 if use_floor else 0) + 3
    per_step = n * (3 + 3 + 1 + 3 + 1 + extra)
    raw = np.fromfile(state_path, dtype=np.float32)
    traj_off = K * per_step
    traj = raw[traj_off:traj_off + (K + 1) * n * 3].reshape(K + 1, n, 3)
    return traj[-1]


def run_bwd(p, K, grad_x_final, h=2.0, mass=1.0, rho=0.3, dt=0.001,
            gravity=-2.0, alpha_dens=1e-3, sim_scale=1.0,
            use_floor=True, floor_y=-10.0):
    """Run xpbd_full_bwd with membranes; return ∂L/∂x_init [n_active, 3]."""
    write_floats(f"{TMP}/xfm_grad_xf.bin", grad_x_final)
    cmd = [BIN, "xpbd_full_bwd",
           str(p['n_active']), str(p['n_static']), str(K),
           str(h), str(mass), f"{rho:.15e}",
           str(dt), str(gravity), f"{alpha_dens:.15e}",
           f"{TMP}/xfm_state.bin",
           f"{TMP}/xfm_pos_static.bin",
           f"{TMP}/xfm_grad_xf.bin",
           f"{TMP}/xfm_grad_x0.bin",
           f"{TMP}/xfm_grad_v0.bin",
           f"{TMP}/xfm_grad_rho.bin",
           f"{sim_scale:.15e}",
           "0.0",
           "0.0",
           "/tmp/xfm_dummy_bonds.bin",
           f"{TMP}/xfm_grad_K.bin",
           f"{TMP}/xfm_grad_vK.bin",
           f"{floor_y:.6f}",
           f"{TMP}/xfm_grad_alpha.bin",
           "0.0",   # restitution
           str(p['n_membranes']), str(p['n_elastic']),
           f"{p['r0']:.15e}",
           f"{TMP}/xfm_mem.bin", f"{TMP}/xfm_pmi.bin"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bwd: {r.stderr}")
    return np.fromfile(f"{TMP}/xfm_grad_x0.bin", dtype=np.float32).reshape(-1, 3)


def main():
    p = setup()
    K = int(os.environ.get('K_STEPS', 1))
    rng = np.random.default_rng(11)

    # Loss = <random_seed, x_final>: linear in x_final, so ∂L/∂x_final = seed.
    grad_seed = rng.standard_normal((p['n_active'], 3)).astype(np.float32)

    use_floor = True
    floor_y = -10.0   # well below all particles → no clamping

    print(f"=== Multi-step xpbd_full bwd FD validation (K={K}, membranes ON) ===")
    x_final = run_fwd(p, K, use_floor=use_floor, floor_y=floor_y)
    L = float(np.dot(grad_seed.flatten(), x_final.flatten()))
    print(f"  x_final[liquid]: {x_final[3]}")
    print(f"  L = <seed, x_final> = {L:.6f}")

    g_an = run_bwd(p, K, grad_seed, use_floor=use_floor, floor_y=floor_y)
    print(f"  analytic ∂L/∂x_init:\n{g_an}")

    # FD: perturb each pos[k, axis], rerun fwd, recompute L.
    eps = 1e-3
    g_fd = np.zeros_like(g_an)
    pos_orig = p['pos'].copy()
    for k in range(p['n_active']):
        for ax in range(3):
            pp = pos_orig.copy(); pp[k, ax] += eps
            write_floats(f"{TMP}/xfm_pos.bin", pp)
            xfp = run_fwd(p, K, use_floor=use_floor, floor_y=floor_y)
            Lp = float(np.dot(grad_seed.flatten(), xfp.flatten()))
            pp = pos_orig.copy(); pp[k, ax] -= eps
            write_floats(f"{TMP}/xfm_pos.bin", pp)
            xfn = run_fwd(p, K, use_floor=use_floor, floor_y=floor_y)
            Ln = float(np.dot(grad_seed.flatten(), xfn.flatten()))
            g_fd[k, ax] = (Lp - Ln) / (2 * eps)
    write_floats(f"{TMP}/xfm_pos.bin", pos_orig)   # restore

    print(f"  FD       ∂L/∂x_init:\n{g_fd}")
    print()

    rel_thresh = 5e-2
    eps_floor = 1e-4
    all_pass = True
    for k in range(p['n_active']):
        denom = max(float(np.linalg.norm(g_an[k])),
                    float(np.linalg.norm(g_fd[k])),
                    eps_floor)
        diff = float(np.linalg.norm(g_an[k] - g_fd[k]))
        rel  = diff / denom
        ok = rel < rel_thresh
        all_pass = all_pass and ok
        tag = "OK" if ok else "FAIL"
        print(f"    [{tag}] k={k}: rel={rel:.3e}  "
              f"||an||={np.linalg.norm(g_an[k]):.3e}  "
              f"||fd||={np.linalg.norm(g_fd[k]):.3e}")

    if all_pass:
        print("\n  [PASS] multi-step xpbd_full backward through membranes matches FD")
        return 0
    else:
        print("\n  [FAIL] multi-step xpbd_full backward through membranes diverges")
        return 1


if __name__ == "__main__":
    sys.exit(main())
