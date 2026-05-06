"""Side-by-side OpenCL vs Metal cube-drop MP4 on unmodified demo1.

Reads two Sibernetic-format position_buffer.txt files (with optional
.times.txt sidecar for variable chunk sizes), samples N matched sim
times with TIME-WARPED sampling — dense in the first descent window,
sparse in the settle phase — so the cube actually appears to fall under
gravity instead of teleporting from top to floor in one frame.

Renderer uses pyvista shaded spheres (not matplotlib scatter) for
properly lit particles, like Sibernetic's official `render_movie.py`.

Usage:
    python3 scripts/render_demo1_parity.py OUT.mp4 OPENCL.txt METAL.txt \\
        [--samples 80] [--descent-fraction 0.5] [--descent-window 0.15]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import imageio.v2 as imageio  # noqa: E402
import pyvista as pv  # noqa: E402


def parse_position_buffer(path: Path):
    """Return (n_active, n_boundary, dt, log_step, boundary_xyz, frames, times)."""
    with open(path) as f:
        header = [f.readline().rstrip('\n') for _ in range(11)]
        n_elastic = int(header[6])
        n_liquid = int(header[7])
        n_boundary = int(header[8])
        dt = float(header[9])
        log_step = int(header[10])
        n_active = n_elastic + n_liquid

        boundary_xyz, active_xyz_first, active_t_first = [], [], []
        for _ in range(n_active + n_boundary):
            ws = f.readline().split()
            x, y, z, t = float(ws[0]), float(ws[1]), float(ws[2]), float(ws[3])
            if t >= 3.0:
                boundary_xyz.append([x, y, z])
            else:
                active_xyz_first.append([x, y, z])
                active_t_first.append(t)
        boundary_xyz = np.array(boundary_xyz)
        frames = [(np.array(active_xyz_first), np.array(active_t_first))]

        rest = [ln for ln in f if ln.strip()]
        if rest:
            chunk_active = n_active
            chunk_full = n_active + n_boundary
            if len(rest) % chunk_active == 0:
                chunk = chunk_active
            elif len(rest) % chunk_full == 0:
                chunk = chunk_full
            else:
                chunk = chunk_active
            for i in range(0, len(rest), chunk):
                xyz, ts = [], []
                for ln in rest[i:i + chunk]:
                    ws = ln.split()
                    t = float(ws[3])
                    if t < 3.0:
                        xyz.append([float(ws[0]), float(ws[1]), float(ws[2])])
                        ts.append(t)
                if len(xyz) == n_active:
                    frames.append((np.array(xyz), np.array(ts)))

    sidecar = path.with_suffix(path.suffix + '.times.txt')
    times = None
    if sidecar.exists():
        with open(sidecar) as f:
            times = [float(ln.strip()) for ln in f if ln.strip()]
        if len(times) != len(frames):
            print(f'WARN: sidecar {sidecar} has {len(times)} times '
                  f'vs {len(frames)} frames — falling back to uniform.')
            times = None
    if times is None:
        frame_dt = dt * log_step
        times = [i * frame_dt for i in range(len(frames))]
    return n_active, n_boundary, dt, log_step, boundary_xyz, frames, times


def time_warped_targets(t_max, n_total, descent_window, descent_fraction,
                          descent_curve='quad'):
    """Return n_total target sim times.

    `descent_fraction` of frames placed in [0, descent_window], remaining in
    [descent_window, t_max]. The descent samples are non-uniformly placed
    according to `descent_curve`:
      - 'linear': evenly spaced in time
      - 'quad'  : t = descent_window * (i/(n-1))^2 — starts dense at t=0,
                  spreads out toward descent_window. Good for capturing the
                  very first physics steps where things move fastest.
    """
    n_descent = max(2, int(round(n_total * descent_fraction)))
    n_settle = max(2, n_total - n_descent)
    if descent_curve == 'quad':
        u = np.linspace(0, 1, n_descent)
        descent_t = descent_window * (u ** 2)
    else:
        descent_t = np.linspace(0, descent_window, n_descent)
    settle_t = np.linspace(descent_window, t_max, n_settle + 1)[1:]
    return np.concatenate([descent_t, settle_t])


def find_index_at_time(times_arr, t):
    return int(np.argmin(np.abs(times_arr - t)))


def setup_plotter(plotter, title, bbox, view_up=(0, 0, 1)):
    """Configure one panel of the side-by-side plot."""
    plotter.set_background('white')
    plotter.add_text(title, position='upper_edge', font_size=12,
                      color='black')
    # Camera position: looking at center of the bottom region from above and to the side
    cx = (bbox[0] + bbox[1]) / 2
    cz = (bbox[4] + bbox[5]) / 2
    cy_center = (bbox[2] + bbox[3]) / 2
    cam_dist = (bbox[1] - bbox[0]) * 1.4
    plotter.camera_position = [
        (cx + cam_dist * 0.7, cy_center + cam_dist * 0.5, cz + cam_dist * 0.7),
        (cx, cy_center, cz),
        (0, 1, 0),
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out', type=Path)
    ap.add_argument('opencl', type=Path)
    ap.add_argument('metal', type=Path)
    ap.add_argument('--samples', type=int, default=80,
                    help="total frames in the MP4")
    ap.add_argument('--descent-window', type=float, default=0.15,
                    help="seconds: dense sampling in [0, this]")
    ap.add_argument('--descent-fraction', type=float, default=0.55,
                    help="fraction of frames in the descent window")
    ap.add_argument('--descent-curve', choices=['linear', 'quad'], default='quad',
                    help="sample placement inside descent window")
    ap.add_argument('--t-max', type=float, default=None,
                    help="cap simulation time shown (seconds). Default = full data.")
    ap.add_argument('--uniform', action='store_true',
                    help="disable time-warping; sample uniformly across [0, t-max]")
    ap.add_argument('--fps', type=int, default=15)
    ap.add_argument('--width', type=int, default=900)
    ap.add_argument('--height', type=int, default=600)
    ap.add_argument('--particle-radius', type=float, default=0.6,
                    help="sphere radius in sim units")
    ap.add_argument('--floor-band', type=float, default=3.0,
                    help="show boundary particles with y < this")
    ap.add_argument('--hide-liquid', action='store_true',
                    help="hide type-1 (liquid) particles, show only elastic")
    args = ap.parse_args()

    print(f"Loading OpenCL: {args.opencl}")
    o_n, o_nb, o_dt, o_ls, o_boundary, o_frames, o_times = parse_position_buffer(args.opencl)
    print(f"  {len(o_frames)} frames over 0..{o_times[-1]:.3f}s")

    print(f"Loading Metal:  {args.metal}")
    m_n, m_nb, m_dt, m_ls, m_boundary, m_frames, m_times = parse_position_buffer(args.metal)
    print(f"  {len(m_frames)} frames over 0..{m_times[-1]:.3f}s")

    t_max = min(o_times[-1], m_times[-1])
    if args.t_max is not None:
        t_max = min(t_max, args.t_max)

    if args.uniform:
        target_times = np.linspace(0, t_max, args.samples)
        print(f"Sampling {len(target_times)} frames uniformly across [0, {t_max*1000:.1f} ms]"
              f" (frame_dt = {t_max/(args.samples-1)*1000:.2f} ms)")
    else:
        target_times = time_warped_targets(t_max, args.samples,
                                            args.descent_window,
                                            args.descent_fraction,
                                            args.descent_curve)
        print(f"Sampling {len(target_times)} frames "
              f"({int(args.samples * args.descent_fraction)} in [0,{args.descent_window}s], "
              f"rest in [{args.descent_window}s, {t_max:.2f}s])")

    o_t = np.array(o_times)
    m_t = np.array(m_times)
    o_indices = [find_index_at_time(o_t, t) for t in target_times]
    m_indices = [find_index_at_time(m_t, t) for t in target_times]

    # Tight bounding box: cube initial region + floor.
    bbox = [
        float(o_boundary[:, 0].min()),
        float(o_boundary[:, 0].max()),
        0.0,
        50.0,           # cube starts at y~44, slight head-room
        float(o_boundary[:, 2].min()),
        float(o_boundary[:, 2].max()),
    ]

    floor_o = o_boundary[o_boundary[:, 1] < args.floor_band]
    floor_m = m_boundary[m_boundary[:, 1] < args.floor_band]
    print(f"  bbox (rendering region): x={bbox[:2]}, y={bbox[2:4]}, z={bbox[4:6]}")
    print(f"  floor-band particles: opencl={len(floor_o)}, metal={len(floor_m)}")

    # Plotter with two viewports.
    plotter = pv.Plotter(off_screen=True, window_size=(args.width, args.height),
                         shape=(1, 2), border=False)

    plotter.subplot(0, 0)
    setup_plotter(plotter, "OpenCL (gold standard)", bbox)
    floor_o_cloud = pv.PolyData(floor_o)
    plotter.add_mesh(floor_o_cloud, color='lightgray', point_size=4,
                      render_points_as_spheres=True, opacity=0.5)

    plotter.subplot(0, 1)
    setup_plotter(plotter, "Metal Native Substrate", bbox)
    floor_m_cloud = pv.PolyData(floor_m)
    plotter.add_mesh(floor_m_cloud, color='lightgray', point_size=4,
                      render_points_as_spheres=True, opacity=0.5)

    # Dummy first-frame particle clouds, will be updated each frame.
    sphere_template = pv.Sphere(radius=args.particle_radius,
                                 theta_resolution=10, phi_resolution=10)

    def filter_xyz(xyz, ts):
        if args.hide_liquid:
            mask = (ts >= 2.0) & (ts < 3.0)  # elastic only
            return xyz[mask]
        return xyz

    o_xyz_f = filter_xyz(*o_frames[0])
    m_xyz_f = filter_xyz(*m_frames[0])

    plotter.subplot(0, 0)
    o_glyph = pv.PolyData(o_xyz_f).glyph(geom=sphere_template, scale=False, factor=1.0)
    o_actor = plotter.add_mesh(o_glyph, color='red', smooth_shading=True)

    plotter.subplot(0, 1)
    m_glyph = pv.PolyData(m_xyz_f).glyph(geom=sphere_template, scale=False, factor=1.0)
    m_actor = plotter.add_mesh(m_glyph, color='red', smooth_shading=True)

    time_text_o = plotter.add_text("", position='lower_left', font_size=10,
                                    color='black', viewport=True)

    print(f"Rendering {len(target_times)} frames to {args.out}")
    writer = imageio.get_writer(str(args.out), fps=args.fps,
                                  codec='libx264', quality=8,
                                  macro_block_size=1)

    for fi, (oi, mi, sim_t) in enumerate(zip(o_indices, m_indices, target_times)):
        o_xyz_f = filter_xyz(*o_frames[oi])
        m_xyz_f = filter_xyz(*m_frames[mi])

        plotter.subplot(0, 0)
        plotter.remove_actor(o_actor)
        o_glyph = pv.PolyData(o_xyz_f).glyph(geom=sphere_template, scale=False, factor=1.0)
        o_actor = plotter.add_mesh(o_glyph, color='red', smooth_shading=True)

        plotter.subplot(0, 1)
        plotter.remove_actor(m_actor)
        m_glyph = pv.PolyData(m_xyz_f).glyph(geom=sphere_template, scale=False, factor=1.0)
        m_actor = plotter.add_mesh(m_glyph, color='red', smooth_shading=True)

        plotter.remove_actor(time_text_o)
        time_text_o = plotter.add_text(f"sim t = {sim_t * 1000:.1f} ms",
                                         position='lower_edge', font_size=12,
                                         color='black', viewport=True)

        img = plotter.screenshot(return_img=True)
        writer.append_data(img)
        if (fi + 1) % 10 == 0:
            print(f"  frame {fi+1}/{len(target_times)} (t={sim_t*1000:.0f}ms)")

    writer.close()
    plotter.close()
    print(f"Saved {args.out}")


if __name__ == '__main__':
    main()
