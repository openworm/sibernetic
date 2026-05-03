"""End-to-end XPBD step test on Metal.

Drops a small cube of "liquid" particles onto a floor with a few static
"boundary" particles. Validates the orchestration: predict → 3-iter
project → velocity update.

Sanity checks:
  - Particles fall under gravity (mean Y decreases over steps).
  - Particles settle on/above the floor (min Y stays >= floor_y).
  - System doesn't explode (no NaN, no escape to infinity).
  - Particle cluster doesn't pancake (extent retention reasonable).

Run:
    .venv/bin/python src/metal_diff/test_xpbd.py
"""
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "sib_metal")


def make_cube(n_side: float, spacing: float, center_y: float):
    """Return [n_side³, 3] positions of a small cube cluster."""
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


def make_bonds_for_cube(n_side: int, spacing: float):
    """Build axis-aligned neighbor bonds for an n_side³ cube grid.

    Returns array of dtype matching the C struct (i:int32, j:int32,
    rest:float32, pad:float32) per bond. For n_side=4: 144 bonds.
    """
    bonds_dtype = np.dtype([
        ('i', np.int32), ('j', np.int32),
        ('rest', np.float32), ('pad', np.float32),
    ])
    bonds = []

    def idx(ix, iy, iz):
        return ix * n_side * n_side + iy * n_side + iz

    for ix in range(n_side):
        for iy in range(n_side):
            for iz in range(n_side):
                if ix < n_side - 1:
                    bonds.append((idx(ix, iy, iz),
                                  idx(ix + 1, iy, iz), spacing, 0.0))
                if iy < n_side - 1:
                    bonds.append((idx(ix, iy, iz),
                                  idx(ix, iy + 1, iz), spacing, 0.0))
                if iz < n_side - 1:
                    bonds.append((idx(ix, iy, iz),
                                  idx(ix, iy, iz + 1), spacing, 0.0))
    return np.array(bonds, dtype=bonds_dtype)


def run_step(n_active, n_static, h, mass, rho_rest, dt, gravity_y,
             floor_y, alpha_dens, n_iters,
             pos_active, vel_active, pos_static,
             bonds=None, alpha_dist=0.0, bench_steps=1):
    p_a = "/tmp/xpbd_in_pos_a.bin"
    p_v = "/tmp/xpbd_in_vel_a.bin"
    p_s = "/tmp/xpbd_in_pos_s.bin"
    p_b = "/tmp/xpbd_in_bonds.bin"
    pos_active.astype(np.float32).tofile(p_a)
    vel_active.astype(np.float32).tofile(p_v)
    pos_static.astype(np.float32).tofile(p_s)
    n_bonds = len(bonds) if bonds is not None else 0
    if n_bonds > 0:
        bonds.tofile(p_b)
    else:
        # Write a single dummy bond so the file exists; n_bonds=0 means
        # the C code skips reading it anyway.
        np.zeros(1, dtype=[('i', np.int32), ('j', np.int32),
                           ('rest', np.float32), ('pad', np.float32)]).tofile(p_b)
    subprocess.run(
        [BINARY, "xpbd_step",
         str(n_active), str(n_static),
         str(h), str(mass), str(rho_rest), str(dt),
         str(gravity_y), str(floor_y), str(alpha_dens),
         str(n_iters),
         p_a, p_v, p_s,
         str(n_bonds), p_b, str(alpha_dist),
         str(bench_steps)],
        check=True,
    )
    pos_out = np.fromfile("/tmp/xpbd_pos_out.bin", dtype=np.float32) \
        .reshape(n_active, 3)
    vel_out = np.fromfile("/tmp/xpbd_vel_out.bin", dtype=np.float32) \
        .reshape(n_active, 3)
    return pos_out, vel_out


def main():
    if not os.path.exists(BINARY):
        print(f"sib_metal binary missing — run {HERE}/build.sh first")
        return 1

    # Small test scene:
    #   - 4×4×4 = 64 active particles forming a cube
    #   - A "floor plate" of 100 static particles spread out at y=0
    #   - Gravity pulls cube down, floor stops it
    n_side = 4
    n_active = n_side ** 3
    h = 1.5
    spacing = 0.6  # particles overlap kernel support
    cube_init_y = 5.0

    active = make_cube(n_side, spacing, cube_init_y)
    vel_active = np.zeros_like(active)

    # Floor plate: 10×10 grid of static particles at y=-0.1
    n_floor_side = 10
    n_static = n_floor_side ** 2
    floor = np.zeros((n_static, 3), dtype=np.float32)
    plate_spacing = spacing
    for ix in range(n_floor_side):
        for iz in range(n_floor_side):
            p = ix * n_floor_side + iz
            floor[p, 0] = (ix - n_floor_side / 2) * plate_spacing
            floor[p, 1] = -0.1
            floor[p, 2] = (iz - n_floor_side / 2) * plate_spacing

    bonds = make_bonds_for_cube(n_side, spacing)

    print(f"=== M7 XPBD substrate test ===")
    print(f"  Active: {n_active} particles in a cube at y={cube_init_y}")
    print(f"  Static: {n_static} floor-plate particles at y=-0.1")
    print(f"  Bonds:  {len(bonds)} elastic axis-aligned bonds (M7.B)")
    print(f"  h={h}, spacing={spacing}")

    init_extent_y = active[:, 1].max() - active[:, 1].min()
    init_mean_y = active[:, 1].mean()
    init_min_y = active[:, 1].min()
    print(f"  Initial: mean_y={init_mean_y:.2f}  min_y={init_min_y:.2f}  "
          f"extent_y={init_extent_y:.2f}")

    # Run N steps and track stats every 50 steps.
    dt = 0.005
    gravity_y = -9.8
    floor_y = 0.0
    alpha_dens = 1e-4   # mild incompressibility
    alpha_dist = 1e-6   # stiff bonds (~1 MPa equivalent compliance)
    n_iters = 3
    rho_rest = 5.0      # tuned for the layout
    mass = 1.0

    print(f"  Stepping: dt={dt}, gravity={gravity_y}, "
          f"alpha_dens={alpha_dens}, alpha_dist={alpha_dist}, "
          f"n_iters={n_iters}, rho_rest={rho_rest}")

    n_steps_per_chunk = 100
    n_chunks = 6  # 600 steps total = 3s sim, well past impact

    cur_pos = active.copy()
    cur_vel = vel_active.copy()
    init_vol = ((cur_pos[:,0].max() - cur_pos[:,0].min()) *
                (cur_pos[:,1].max() - cur_pos[:,1].min()) *
                (cur_pos[:,2].max() - cur_pos[:,2].min()))

    print()
    print(f"  step  mean_y  min_y  ext_x  ext_y  ext_z  vol  vol/init  max|v|")

    for chunk in range(n_chunks):
        pos, vel = run_step(
            n_active, n_static, h, mass, rho_rest, dt, gravity_y,
            floor_y, alpha_dens, n_iters,
            cur_pos, cur_vel, floor,
            bonds=bonds, alpha_dist=alpha_dist,
            bench_steps=n_steps_per_chunk,
        )
        cur_step = (chunk + 1) * n_steps_per_chunk
        has_nan = bool(np.isnan(pos).any() or np.isnan(vel).any())
        ex = pos[:,0].max() - pos[:,0].min()
        ey = pos[:,1].max() - pos[:,1].min()
        ez = pos[:,2].max() - pos[:,2].min()
        vol = ex * ey * ez
        print(f"  {cur_step:4d}  {pos[:,1].mean():6.2f} {pos[:,1].min():6.2f}  "
              f"{ex:5.2f}  {ey:5.2f}  {ez:5.2f}  {vol:5.2f}  "
              f"{vol/init_vol:7.2f}x  {np.abs(vel).max():.2f}")
        cur_pos, cur_vel = pos, vel
        if has_nan:
            print("  STOPPED: NaN detected")
            break

    print()
    final_ex = cur_pos[:, 0].max() - cur_pos[:, 0].min()
    final_ey = cur_pos[:, 1].max() - cur_pos[:, 1].min()
    final_ez = cur_pos[:, 2].max() - cur_pos[:, 2].min()
    final_vol = final_ex * final_ey * final_ez
    final_mean_y = cur_pos[:, 1].mean()
    final_min_y = cur_pos[:, 1].min()

    has_nan = bool(np.isnan(cur_pos).any() or np.isnan(cur_vel).any())
    no_explosion = bool(np.abs(cur_pos).max() < 100.0)
    fell_significantly = bool(final_mean_y < init_mean_y - 1.0)
    above_floor = bool(final_min_y > floor_y - 0.5)
    # Cube integrity: final volume should be within 0.5x to 2.5x of initial
    # (some elastic deformation is expected, but no pancaking or explosion).
    integrity = bool(0.5 < final_vol / init_vol < 2.5)

    print(f"  Final: mean_y={final_mean_y:.2f}  min_y={final_min_y:.2f}")
    print(f"         extent x={final_ex:.2f} y={final_ey:.2f} z={final_ez:.2f}  "
          f"vol={final_vol:.2f} ({final_vol/init_vol:.2f}x initial)")
    print(f"  Checks:")
    print(f"    no NaN:                  {'[OK]' if not has_nan else '[FAIL]'}")
    print(f"    no explosion (<100m):    {'[OK]' if no_explosion else '[FAIL]'}")
    print(f"    fell significantly:      {'[OK]' if fell_significantly else '[FAIL]'}")
    print(f"    stayed above floor:      {'[OK]' if above_floor else '[FAIL]'}")
    print(f"    cube integrity (0.5–2.5x vol):  "
          f"{'[OK]' if integrity else '[FAIL]'}")

    all_ok = (not has_nan and no_explosion and fell_significantly and
              above_floor and integrity)
    print()
    print(f"  {'[PASS] M7.A+B: XPBD with elastic bonds keeps cube intact post-impact' if all_ok else '[FAIL] M7.A+B: see failures above'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
