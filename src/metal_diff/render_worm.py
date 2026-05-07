"""Render a worm-config trajectory to MP4.

Designed for `worm_alone_half_resolution` and similar small worm configs
(2 290 elastic + 388 liquid particles in this case). Renders:

- Worm body as a closed mesh from the [membranes] triangle list.
- Liquid particles as a small blue point cloud.
- Box outline (lightgray).

Auto-detects OpenCL vs Metal trajectory format using the same heuristic
as `render_one_sprig.py.parse()`.
"""
import argparse
import os
from pathlib import Path
import numpy as np

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
import imageio.v2 as imageio
import pyvista as pv


def _count_lines(path):
    n = 0
    with open(path) as f:
        for _ in f: n += 1
    return n


def parse_trajectory(path):
    """Return (active_traj [n_frames, n_active, 4], boundary_pts, frame_dt_ms,
    metadata dict {n_e, n_l, n_b, dt, ls})."""
    total = _count_lines(path)
    with open(path) as f:
        for _ in range(6): f.readline()
        n_e = int(f.readline()); n_l = int(f.readline()); n_b = int(f.readline())
        dt = float(f.readline()); ls = float(f.readline())
        n_a = n_e + n_l
        frame_dt_ms = dt * ls * 1000

        remaining = total - 11 - n_a - n_b
        if remaining > 0 and (n_a + n_b) > 1 and remaining % (n_a + n_b) == 0:
            includes_boundary_per_frame = True
            n_extra = remaining // (n_a + n_b)
        else:
            includes_boundary_per_frame = False
            n_extra = remaining // max(n_a, 1)

        traj = []
        # Frame 0: n_a active + n_b boundary
        f0 = []
        for _ in range(n_a):
            line = f.readline()
            if not line: break
            w = line.split()
            f0.append([float(w[0]), float(w[1]), float(w[2]), float(w[3])])
        traj.append(f0)

        boundary = []
        for _ in range(n_b):
            line = f.readline()
            if not line: break
            w = line.split()
            if len(w) >= 3:
                boundary.append([float(w[0]), float(w[1]), float(w[2])])
        boundary = np.array(boundary) if boundary else np.zeros((0, 3))

        # Subsequent frames
        for _ in range(n_extra):
            fk = []
            ok = True
            for _ in range(n_a):
                line = f.readline()
                if not line: ok = False; break
                w = line.split()
                if len(w) < 4: ok = False; break
                fk.append([float(w[0]), float(w[1]), float(w[2]), float(w[3])])
            if not ok or len(fk) != n_a: break
            traj.append(fk)
            if includes_boundary_per_frame:
                for _ in range(n_b):
                    line = f.readline()
                    if not line: break

    arr = np.array(traj)
    fmt = "Metal" if includes_boundary_per_frame else "OpenCL"
    print(f'parsed: {n_a} active ({n_e} elastic, {n_l} liquid), {n_b} boundary, '
          f'{arr.shape[0]} frames, frame_dt={frame_dt_ms:.4f} ms, format={fmt}')
    return arr, boundary, frame_dt_ms, {'n_e': n_e, 'n_l': n_l, 'n_b': n_b,
                                          'dt': dt, 'ls': ls}


def load_membranes(scenario, repo_root='/Users/slarson/Documents/sibernetic'):
    """Read the (n_membranes, 3) int32 array from the .bin emitted by load_config."""
    path = Path(repo_root) / 'src' / 'metal_diff' / 'cache' / f'{scenario}_membranes.bin'
    if not path.exists():
        # Fall back to parsing the original config
        cfg = Path(repo_root) / 'configuration' / scenario
        tris = []
        in_mem = False
        with open(cfg) as f:
            for line in f:
                s = line.strip()
                if s == '[membranes]': in_mem = True; continue
                if s.startswith('[') and s.endswith(']'): in_mem = False; continue
                if in_mem and s:
                    parts = s.split()
                    if len(parts) == 3:
                        tris.append([int(parts[0]), int(parts[1]), int(parts[2])])
        return np.array(tris, dtype=np.int32) if tris else np.zeros((0, 3), dtype=np.int32)
    return np.fromfile(str(path), dtype=np.int32).reshape(-1, 3)


def render(traj, boundary, frame_dt_ms, membranes, n_e, n_l, out_path, *,
           fps=30, width=900, height=900, title='worm', box_min=(0,0,0),
           box_max=(100.2, 50.1, 267.2), max_frames=0, zoom=1.0,
           focus_on_worm=True, view='iso', z_focus_frac=0.5,
           z_extent_frac=1.0, show_edges=False):
    """Render worm trajectory.

    view options:
      'iso'      — 3/4 view of the entire worm length (default)
      'closeup'  — closeup on a fraction of the worm centered at z_focus_frac;
                   width controlled by z_extent_frac (1.0 = whole worm).
      'cross'    — view perpendicular to the worm's z axis, showing the
                   cross-section (cylindrical diameter under gravity).
    """
    if max_frames > 0 and traj.shape[0] > max_frames:
        traj = traj[:max_frames]
    n_frames = traj.shape[0]

    box_min = np.array(box_min); box_max = np.array(box_max)
    cx, cy, cz = (box_min + box_max) / 2

    plotter = pv.Plotter(off_screen=True, window_size=(width, height))
    plotter.set_background('white')
    plotter.add_text(title, position='upper_edge', font_size=10, color='black',
                     viewport=True)

    # Stronger lighting for 3D depth perception.
    # Pyvista's default scene has ambient lighting that flattens out the
    # surface. We replace it with a key + fill + back-light setup.
    plotter.remove_all_lights()
    plotter.add_light(pv.Light(position=(1.0, 1.0, 0.5), focal_point=(0, 0, 0),
                                color='white', intensity=0.85,
                                light_type='scenelight'))
    plotter.add_light(pv.Light(position=(-0.8, 0.4, -0.6), focal_point=(0, 0, 0),
                                color='white', intensity=0.45,
                                light_type='scenelight'))
    plotter.add_light(pv.Light(position=(0.0, -0.8, 0.3), focal_point=(0, 0, 0),
                                color='white', intensity=0.30,
                                light_type='scenelight'))

    # Camera target: when focus_on_worm, frame the elastic body's centroid
    # at frame 0, and use the worm's actual extent to size the camera distance
    # rather than the (mostly-empty) sim box.
    if focus_on_worm:
        elastic0 = traj[0, :n_e, :3]
        worm_min = elastic0.min(axis=0)
        worm_max = elastic0.max(axis=0)
        cx, cy, cz = (worm_min + worm_max) / 2
        worm_size = worm_max - worm_min
        long_dim = worm_size.max()
        # Cross-section diameter (xy plane, ignoring z length)
        cross_dim = max(worm_size[0], worm_size[1])
    else:
        long_dim = (box_max - box_min).max()
        cross_dim = long_dim
        worm_min = box_min; worm_max = box_max

    if view == 'closeup' or view == 'cross':
        # Center on a fraction of the worm length
        z_min, z_max = worm_min[2], worm_max[2]
        cz = z_min + (z_max - z_min) * z_focus_frac
        # View extent along z = z_extent_frac × full worm length
        view_long = (z_max - z_min) * z_extent_frac
        if view == 'cross':
            # 3/4 view from near one end of the visible section, looking
            # along the worm's z axis. Reveals cylindrical cross-section
            # plus a length of body for context. Distance scaled by the
            # whichever-is-larger of cross-section diameter and the
            # length-fraction in view.
            view_extent = max(cross_dim * 4.0, view_long * 1.5)
            cam_dist = view_extent / max(zoom, 1e-3)
            plotter.camera_position = [
                (cx + cam_dist * 0.30, cy + cam_dist * 0.20, cz + cam_dist * 0.85),
                (cx, cy, cz),
                (0, 1, 0),
            ]
        else:
            # Closeup: side-on perspective, slightly elevated, framing the
            # selected section. Cam distance scaled so the long axis
            # (z-extent) just fits at zoom=1.
            cam_dist = view_long * 2.4 / max(zoom, 1e-3)
            plotter.camera_position = [
                (cx + cam_dist * 0.50, cy + cam_dist * 0.20, cz + cam_dist * 0.40),
                (cx, cy, cz),
                (0, 1, 0),
            ]
    else:
        # iso (default): 3/4 view of the whole worm
        cam_dist = long_dim * 1.6 / max(zoom, 1e-3)
        plotter.camera_position = [
            (cx + cam_dist * 0.55, cy + cam_dist * 0.35, cz + cam_dist * 0.55),
            (cx, cy, cz),
            (0, 1, 0),
        ]

    # Box outline
    box = pv.Box(bounds=[box_min[0], box_max[0], box_min[1], box_max[1],
                          box_min[2], box_max[2]]).extract_feature_edges()
    plotter.add_mesh(box, color='lightgray', line_width=1, opacity=0.3)

    # Pre-build a faces array for the membrane mesh.
    # PyVista format: [3, i, j, k, 3, i, j, k, ...] flattened.
    if len(membranes):
        faces_flat = np.empty((membranes.shape[0], 4), dtype=np.int64)
        faces_flat[:, 0] = 3
        faces_flat[:, 1:] = membranes
        faces_flat = faces_flat.flatten()
    else:
        faces_flat = None

    actor_mem = None
    actor_liq = None
    time_text = plotter.add_text('t = 0', position='lower_edge', font_size=12,
                                  color='black', viewport=True)

    writer = imageio.get_writer(str(out_path), fps=fps, codec='libx264',
                                 quality=8, macro_block_size=1)

    for fi in range(n_frames):
        active = traj[fi, :, :3]   # (n_e + n_l, 3)
        elastic_pts = active[:n_e]
        liquid_pts = active[n_e:n_e + n_l]

        # Re-build membrane mesh each frame (positions change)
        if actor_mem is not None: plotter.remove_actor(actor_mem)
        if actor_liq is not None: plotter.remove_actor(actor_liq)
        if faces_flat is not None:
            mem = pv.PolyData(elastic_pts, faces_flat)
            actor_mem = plotter.add_mesh(mem, color='lightcoral', opacity=0.85,
                                          smooth_shading=True,
                                          show_edges=show_edges,
                                          edge_color='maroon',
                                          line_width=0.4,
                                          specular=0.6, specular_power=20,
                                          ambient=0.20, diffuse=0.85)

        if len(liquid_pts):
            liq = pv.PolyData(liquid_pts)
            actor_liq = plotter.add_mesh(liq, color='royalblue', point_size=4,
                                          render_points_as_spheres=True,
                                          opacity=0.85)

        plotter.remove_actor(time_text)
        t_ms = fi * frame_dt_ms
        # Quick parity-friendly summary
        ey_mean = float(elastic_pts[:, 1].mean())
        time_text = plotter.add_text(
            f't = {t_ms:.3f} ms   <y_elastic> = {ey_mean:.3f}',
            position='lower_edge', font_size=12, color='black', viewport=True)

        img = plotter.screenshot(return_img=True)
        writer.append_data(img)
        if (fi + 1) % 25 == 0:
            print(f'    frame {fi+1}/{n_frames}')

    writer.close()
    plotter.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--scenario', default='worm_alone_half_resolution',
                    help='used to locate the membrane bin')
    ap.add_argument('--fps', type=int, default=30)
    ap.add_argument('--width', type=int, default=900)
    ap.add_argument('--height', type=int, default=900)
    ap.add_argument('--title', default='worm')
    ap.add_argument('--max-frames', type=int, default=0)
    ap.add_argument('--box-min', nargs=3, type=float, default=[0, 0, 0])
    ap.add_argument('--box-max', nargs=3, type=float, default=[100.2, 50.1, 267.2])
    ap.add_argument('--zoom', type=float, default=1.0,
                    help='>1 = closer in, <1 = farther out')
    ap.add_argument('--no-focus-worm', action='store_true',
                    help='Frame the full sim box instead of just the worm')
    ap.add_argument('--view', choices=['iso', 'closeup', 'cross'], default='iso',
                    help='iso=3/4 of whole worm, closeup=section centered at '
                         'z_focus_frac, cross=perpendicular to worm length')
    ap.add_argument('--z-focus-frac', type=float, default=0.5,
                    help='where along worm length to center closeup/cross view '
                         '(0=tail, 1=head)')
    ap.add_argument('--z-extent-frac', type=float, default=1.0,
                    help='fraction of worm length visible in closeup view '
                         '(smaller = tighter)')
    ap.add_argument('--show-edges', action='store_true',
                    help='draw membrane triangle edges (helps depth perception)')
    args = ap.parse_args()

    traj, boundary, frame_dt, meta = parse_trajectory(args.input)
    membranes = load_membranes(args.scenario)
    print(f'  membranes: {membranes.shape[0]} triangles')

    render(traj, boundary, frame_dt, membranes, meta['n_e'], meta['n_l'],
           args.output, fps=args.fps, width=args.width, height=args.height,
           title=args.title, box_min=tuple(args.box_min),
           box_max=tuple(args.box_max), max_frames=args.max_frames,
           zoom=args.zoom, focus_on_worm=not args.no_focus_worm,
           view=args.view, z_focus_frac=args.z_focus_frac,
           z_extent_frac=args.z_extent_frac, show_edges=args.show_edges)
    print(f'Saved {args.output}')


if __name__ == '__main__':
    main()
