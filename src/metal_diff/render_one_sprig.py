"""Render a one_sprig_test position trajectory to MP4.

Handles both trajectory formats:
- OpenCL Sibernetic: header + frame0(active+boundary) + frame_k(active only)
- Native Metal (dump_metal_trajectory.py): header + frame_k(active+boundary)
  for ALL frames including frame 0

Format is auto-detected from file size after parsing the header. The
elastic particle is the first active particle; the anchor is at
(elastic_x_init, anchor_y, elastic_z_init).
"""
import argparse
import os
import sys
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


def parse(path):
    """Parse position_buffer; return (elastic_traj, boundary_pts, frame_dt_ms)."""
    total = _count_lines(path)
    with open(path) as f:
        # Header: 6 box-bound floats, n_e, n_l, n_b, dt, logstep
        for _ in range(6): f.readline()
        n_e = int(f.readline()); n_l = int(f.readline()); n_b = int(f.readline())
        dt = float(f.readline()); ls = float(f.readline())
        n_a = n_e + n_l
        frame_dt_ms = dt * ls * 1000

        # Auto-detect: does every frame include boundary, or only frame 0?
        # OpenCL: total = 11 + n_a + n_b + n_extra * n_a
        # Metal:  total = 11 + (n_a + n_b) * n_total_frames
        # Discriminator: divisibility of "remaining" by (n_a + n_b).
        # If n_b >> n_a, this rarely matches by accident.
        remaining = total - 11 - n_a - n_b
        if remaining > 0 and (n_a + n_b) > 1 and remaining % (n_a + n_b) == 0:
            includes_boundary_per_frame = True
            n_extra_frames = remaining // (n_a + n_b)
        else:
            includes_boundary_per_frame = False
            n_extra_frames = remaining // max(n_a, 1)

        # Frame 0
        active_traj = []  # list of (n_a, 4) per frame
        f0 = []
        for _ in range(n_a):
            w = f.readline().split()
            f0.append([float(w[0]), float(w[1]), float(w[2]), float(w[3])])
        active_traj.append(f0)

        boundary_pts = []
        for _ in range(n_b):
            w = f.readline().split()
            if len(w) < 3: continue
            boundary_pts.append([float(w[0]), float(w[1]), float(w[2])])
        boundary_pts = np.array(boundary_pts) if boundary_pts else np.zeros((0, 3))

        # Subsequent frames
        for fi in range(n_extra_frames):
            fk = []
            for _ in range(n_a):
                line = f.readline()
                if not line: break
                w = line.split()
                if len(w) < 4: break
                fk.append([float(w[0]), float(w[1]), float(w[2]), float(w[3])])
            if len(fk) != n_a: break
            active_traj.append(fk)
            if includes_boundary_per_frame:
                # skip n_b boundary lines (we already have them from frame 0)
                for _ in range(n_b):
                    line = f.readline()
                    if not line: break

    traj = np.array(active_traj)  # (n_frames, n_a, 4)
    print(f'parsed: {n_a} active, {n_b} boundary, {traj.shape[0]} frames, '
          f'frame_dt={frame_dt_ms:.4f} ms, '
          f'format={"Metal" if includes_boundary_per_frame else "OpenCL"}')
    return traj, boundary_pts, frame_dt_ms


def render(traj, boundary, frame_dt_ms, out_path, *, fps=30,
           width=900, height=900, anchor_y=32.565, title='one_sprig_test',
           box_min=(0, 0, 0), box_max=(33.4, 33.4, 33.4),
           y_label_anchor='auto'):
    """Render the trajectory of the first active particle as the elastic.

    traj shape: (n_frames, n_a, 4) — only traj[:, 0] is rendered as elastic.
    """
    elastic = traj[:, 0, :3]   # (n_frames, 3)
    print(f'  elastic y range: [{elastic[:,1].min():.3f}, {elastic[:,1].max():.3f}]')
    print(f'  elastic x range: [{elastic[:,0].min():.4f}, {elastic[:,0].max():.4f}]')
    print(f'  elastic z range: [{elastic[:,2].min():.4f}, {elastic[:,2].max():.4f}]')

    box_min = np.array(box_min); box_max = np.array(box_max)
    cx, cy, cz = (box_min + box_max) / 2

    plotter = pv.Plotter(off_screen=True, window_size=(width, height))
    plotter.set_background('white')
    # Title rendered with smaller font so long titles fit
    plotter.add_text(title, position='upper_edge', font_size=10, color='black',
                     viewport=True)

    cam_dist = (box_max - box_min).max() * 1.6
    plotter.camera_position = [
        (cx + cam_dist * 0.7, cy + cam_dist * 0.3, cz + cam_dist * 1.0),
        (cx, cy, cz),
        (0, 1, 0),
    ]

    # Box outline
    box_lines = pv.Box(bounds=[box_min[0], box_max[0], box_min[1], box_max[1],
                                box_min[2], box_max[2]]).extract_feature_edges()
    plotter.add_mesh(box_lines, color='lightgray', line_width=1, opacity=0.4)

    # Anchor (gray sphere) — position in (x_init, anchor_y, z_init)
    anchor_pos = np.array([elastic[0, 0], anchor_y, elastic[0, 2]])
    anchor_sphere = pv.Sphere(center=anchor_pos, radius=0.5,
                              theta_resolution=12, phi_resolution=12)
    plotter.add_mesh(anchor_sphere, color='dimgray')

    # Floor scatter (just the bottom-row boundary, faint)
    if len(boundary):
        floor = boundary[boundary[:, 1] < 1.0]
        if len(floor):
            plotter.add_mesh(pv.PolyData(floor), color='lightgray',
                             point_size=2, render_points_as_spheres=True, opacity=0.3)

    sphere_template = pv.Sphere(radius=1.5, theta_resolution=16, phi_resolution=16)
    actor_elastic = None
    actor_spring = None
    time_text = plotter.add_text('t = 0.000 ms', position='lower_edge',
                                  font_size=12, color='black', viewport=True)

    writer = imageio.get_writer(str(out_path), fps=fps, codec='libx264',
                                 quality=8, macro_block_size=1)

    for fi in range(elastic.shape[0]):
        ep = elastic[fi]
        if actor_elastic is not None: plotter.remove_actor(actor_elastic)
        if actor_spring is not None: plotter.remove_actor(actor_spring)
        elastic_sphere = pv.Sphere(center=ep, radius=1.5,
                                    theta_resolution=16, phi_resolution=16)
        actor_elastic = plotter.add_mesh(elastic_sphere, color='royalblue',
                                          smooth_shading=True)
        spring_line = pv.Line(anchor_pos, ep)
        actor_spring = plotter.add_mesh(spring_line, color='black', line_width=2)

        plotter.remove_actor(time_text)
        t_ms = fi * frame_dt_ms
        time_text = plotter.add_text(
            f't = {t_ms:.3f} ms   y = {ep[1]:.3f}',
            position='lower_edge', font_size=12, color='black', viewport=True)

        img = plotter.screenshot(return_img=True)
        writer.append_data(img)
        if (fi + 1) % 100 == 0:
            print(f'    frame {fi+1}/{elastic.shape[0]}')

    writer.close()
    plotter.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--fps', type=int, default=30)
    ap.add_argument('--width', type=int, default=900)
    ap.add_argument('--height', type=int, default=900)
    ap.add_argument('--anchor-y', type=float, default=32.565)
    ap.add_argument('--title', default='one_sprig_test')
    ap.add_argument('--max-frames', type=int, default=0,
                    help='Truncate to first N frames (0 = all)')
    args = ap.parse_args()

    traj, boundary, frame_dt = parse(args.input)
    if args.max_frames > 0 and traj.shape[0] > args.max_frames:
        traj = traj[:args.max_frames]
        print(f'  truncated to {traj.shape[0]} frames')
    render(traj, boundary, frame_dt, args.output,
           fps=args.fps, width=args.width, height=args.height,
           anchor_y=args.anchor_y, title=args.title)
    print(f'Saved {args.output}')


if __name__ == '__main__':
    main()
