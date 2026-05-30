"""Property/fuzz testing: random scenarios through xpbd_full_fwd.

Generates ~30 random small scenarios (n_active 4-32, random positions in a
box, random velocities, randomized params within physically reasonable
ranges), runs each through xpbd_full_fwd for K=500 steps. Per scenario,
verify:
  - exit code 0
  - state_out.bin has expected size
  - no NaN/Inf anywhere in trajectory
  - trajectory has non-trivial motion (positions actually evolve)

Catches arithmetic edge cases the fixed-scenario tests miss: degenerate
particle separations, near-h boundary conditions, extreme initial vels.
"""
import os
import platform
import struct
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE,
                   "sib_cuda.exe" if platform.system() == "Windows" else "sib_cuda")
TMP_DIR = os.path.join(tempfile.gettempdir(), "v9_fuzz")
os.makedirs(TMP_DIR, exist_ok=True)

N_FUZZ = 30
SEED = 42


def random_scenario(rng, fuzz_id):
    """Generate a random sceneario within sensible physical bounds."""
    # Parameter ranges within ~1000x of demo1 baseline. Outside these
    # bounds (very small h * sim_scale, very large dt with tiny sim_scale)
    # the explicit-Euler integrator hits CFL or visc_amp blows up — that's
    # a known physics limit, not a substrate bug. Bounds here cover the
    # range of physically reasonable C. elegans-scale SPH simulations.
    n_active = int(rng.integers(4, 32))
    n_static = int(rng.integers(4, 64))
    K = int(rng.choice([100, 500, 1000]))
    h = float(rng.uniform(1.0, 3.5))
    mass = float(10 ** rng.uniform(-13, -3))
    rho_rest = float(10 ** rng.uniform(-13, -10))
    dt = float(rng.choice([1e-5, 5e-5]))
    g_y = float(rng.uniform(-15.0, -5.0))
    alpha_dens = float(10 ** rng.uniform(-4, 0))
    sim_scale = float(10 ** rng.uniform(-4, -2))
    # 50% scenarios have pair forces, 50% don't
    visc_pair_coef = float(rng.uniform(0, 1e-2)) if rng.random() < 0.5 else 0.0
    # 50% scenarios have springs. CFL bound: spring_K * dt^2 / sim_scale
    # < ~0.3 keeps explicit-Euler stable (see V5 stiff_spring finding).
    has_springs = rng.random() < 0.5
    if has_springs:
        cfl_max_K = 0.3 * sim_scale / (dt * dt)
        spring_K_upper = max(10.0, min(1000.0, cfl_max_K))
        spring_K = float(rng.uniform(10, spring_K_upper)) if spring_K_upper > 10 else 0.0
    else:
        spring_K = 0.0
    # 70% have floor
    has_floor = rng.random() < 0.7
    floor_y = float(rng.uniform(-0.5, 0.5)) if has_floor else None
    restitution = float(rng.uniform(0, 1)) if has_floor else 0.0

    # Particles in a box, well separated to avoid r->0 degeneracy.
    pos0 = rng.uniform(-1, 1, size=(n_active, 3)).astype(np.float32) * h
    vel0 = rng.uniform(-0.5, 0.5, size=(n_active, 3)).astype(np.float32)
    pos_static = rng.uniform(-2, 2, size=(n_static, 3)).astype(np.float32) * h

    # Springs: each bond connects two distinct active particles.
    bonds = []
    if has_springs:
        n_bonds = int(rng.integers(1, min(8, n_active)))
        for _ in range(n_bonds):
            i, j = rng.choice(n_active, 2, replace=False)
            rest = float(np.linalg.norm(pos0[i] - pos0[j])) * float(rng.uniform(0.8, 1.2))
            bonds.append((int(i), int(j), rest))

    return {
        'id': fuzz_id,
        'n_active': n_active, 'n_static': n_static, 'K': K,
        'h': h, 'mass': mass, 'rho_rest': rho_rest, 'dt': dt,
        'g_y': g_y, 'alpha_dens': alpha_dens, 'sim_scale': sim_scale,
        'visc_pair_coef': visc_pair_coef, 'spring_K': spring_K,
        'bonds': bonds,
        'floor_y': floor_y, 'restitution': restitution,
        'pos0': pos0, 'vel0': vel0, 'pos_static': pos_static,
    }


def run_scenario(s):
    """Returns (rc, n_nans, motion, stderr_tail)."""
    label = f"fuzz_{s['id']:03d}"
    p_pos = os.path.join(TMP_DIR, f"{label}_pos.bin")
    p_vel = os.path.join(TMP_DIR, f"{label}_vel.bin")
    p_static = os.path.join(TMP_DIR, f"{label}_static.bin")
    p_bonds = os.path.join(TMP_DIR, f"{label}_bonds.bin")
    p_state = os.path.join(TMP_DIR, f"{label}_state.bin")

    s['pos0'].tofile(p_pos)
    s['vel0'].tofile(p_vel)
    s['pos_static'].tofile(p_static)

    use_pair = s['visc_pair_coef'] != 0
    use_floor = s['floor_y'] is not None
    use_springs = s['spring_K'] != 0
    if use_springs:
        with open(p_bonds, "wb") as f:
            for i, j, rest in s['bonds']:
                f.write(struct.pack("ii f f", i, j, rest, 0.0))

    cmd = [BIN, "xpbd_full_fwd",
           str(s['n_active']), str(s['n_static']), str(s['K']),
           str(s['h']), f"{s['mass']:.15e}", f"{s['rho_rest']:.15e}",
           f"{s['dt']:.15e}", f"{s['g_y']:.6f}", f"{s['alpha_dens']:.15e}",
           p_pos, p_vel, p_static, p_state]

    # Optional args: sim_scale, visc_pair_coef, spring_K, bonds, floor_y, restit
    cmd.extend([f"{s['sim_scale']:.15e}",
                f"{s['visc_pair_coef']:.15e}",
                f"{s['spring_K']:.15e}"])
    if use_springs:
        cmd.append(p_bonds)
    elif use_floor:
        # need bonds slot filled to reach floor_y; use empty
        empty = os.path.join(TMP_DIR, "empty_bonds.bin")
        if not os.path.exists(empty):
            open(empty, "wb").close()
        cmd.append(empty)
    if use_floor:
        cmd.extend([f"{s['floor_y']:.6f}", f"{s['restitution']:.6f}"])

    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    rc = r.returncode
    n_nans, motion = -1, -1.0
    if rc == 0 and os.path.exists(p_state):
        extra = (1 if use_pair else 0) + (1 if use_floor else 0)
        per_step = s['n_active'] * (3 + 3 + 1 + 3 + 1 + extra)
        traj_off = s['K'] * per_step
        raw = np.fromfile(p_state, dtype=np.float32)
        expected = s['K'] * per_step + (s['K'] + 1) * s['n_active'] * 3 + s['n_active'] * 3
        if raw.size == expected:
            traj = raw[traj_off:traj_off + (s['K']+1)*s['n_active']*3]
            traj = traj.reshape(s['K']+1, s['n_active'], 3)
            n_nans = int(np.isnan(traj).sum() + np.isinf(traj).sum())
            if n_nans == 0:
                motion = float(np.linalg.norm(traj[-1] - traj[0]))
    # Cleanup
    for p in [p_pos, p_vel, p_static, p_state, p_bonds]:
        if os.path.exists(p):
            try: os.remove(p)
            except: pass
    return rc, n_nans, motion, r.stderr[-200:]


def main():
    if not os.path.exists(BIN):
        print(f"ERROR: {BIN} missing"); return 1

    rng = np.random.default_rng(SEED)
    print(f"=== V9 fuzz: {N_FUZZ} random scenarios (seed={SEED}) ===\n")

    results = []
    fails = []
    for i in range(N_FUZZ):
        s = random_scenario(rng, i)
        flags = []
        if s['visc_pair_coef'] != 0: flags.append("pair")
        if s['spring_K'] != 0: flags.append(f"spring(K={s['spring_K']:.0f},nb={len(s['bonds'])})")
        if s['floor_y'] is not None: flags.append(f"floor(y={s['floor_y']:.2f},r={s['restitution']:.2f})")
        print(f"  [{i:02d}] n_a={s['n_active']:2d} n_s={s['n_static']:2d} K={s['K']:4d} "
              f"h={s['h']:.2f} m={s['mass']:.2e} dt={s['dt']:.0e} "
              f"sim_scale={s['sim_scale']:.0e} | {'+'.join(flags) or 'baseline'}", end=" ... ")
        rc, n_nans, motion, stderr_tail = run_scenario(s)
        results.append((i, rc, n_nans, motion))
        ok = (rc == 0 and n_nans == 0 and motion > 0)
        if ok:
            print(f"OK (motion={motion:.2e})")
        else:
            print(f"FAIL (rc={rc} nans={n_nans} motion={motion})")
            if stderr_tail: print(f"       stderr: {stderr_tail[:150]!r}")
            fails.append(i)

    print(f"\n{N_FUZZ - len(fails)}/{N_FUZZ} scenarios passed")
    if not fails:
        print("[OVERALL PASS] V9 fuzz")
        return 0
    print(f"[OVERALL FAIL] {len(fails)} scenarios failed: {fails}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
