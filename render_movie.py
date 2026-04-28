#!/usr/bin/env python3
"""Render simulation buffers to MP4 video using pyvista."""

import pyvista as pv
import sys
import os
import numpy as np

def render_to_movie(position_file, output_file="worm_simulation.mp4", fps=30,
                    frame_skip=10, max_frames=500, include_boundary=False,
                    camera="side"):
    """Load position buffer and render to MP4.

    Args:
        frame_skip: Only keep every Nth frame (default 10)
        max_frames: Maximum frames to render (default 500)
        include_boundary: If True, render boundary (type 3) particles too
        camera: "side" (default, original behaviour) or "iso" (3/4 view from
            above-front-side so the floor plane is visible).
    """

    colours = {1.1: "blue", 2.1: "green", 2.2: "turquoise", 3: "#cccccc"}

    points_per_frame = []
    types_per_frame = []

    line_count = 0
    pcount = 0
    points = []
    types = []

    numOfElasticP = 0
    numOfLiquidP = 0
    numOfBoundaryP = 0
    time_step = None
    log_step = None
    frame_count = 0

    print(f"Loading {position_file} (skip={frame_skip}, max={max_frames})...")

    for line in open(position_file):
        ws = line.split()

        if line_count == 6:
            numOfElasticP = int(ws[0])
        if line_count == 7:
            numOfLiquidP = int(ws[0])
        if line_count == 8:
            numOfBoundaryP = int(ws[0])
        if line_count == 9:
            time_step = float(ws[0])
        if line_count == 10:
            log_step = int(ws[0])

        if len(ws) == 4:
            ptype = float(ws[3])

            if include_boundary or ptype != 3:
                points.append([float(ws[0]), float(ws[1]), float(ws[2])])
                types.append(ptype)

        if log_step is not None:
            pcount += 1

            total_particles = numOfBoundaryP + numOfElasticP + numOfLiquidP
            if pcount == total_particles:
                # Only keep every Nth frame
                if frame_count % frame_skip == 0 and len(points_per_frame) < max_frames:
                    points_per_frame.append(np.array(points))
                    types_per_frame.append(np.array(types))
                    if len(points_per_frame) % 50 == 0:
                        print(f"  Loaded {len(points_per_frame)} frames...")
                points = []
                types = []
                pcount = 0
                frame_count += 1

                # Early exit if we have enough frames
                if len(points_per_frame) >= max_frames:
                    print(f"  Reached max_frames={max_frames}, stopping load.")
                    break

        line_count += 1

    num_frames = len(points_per_frame)
    print(f"Loaded {num_frames} frames")
    print(f"Time step: {time_step}, Log step: {log_step}")

    if num_frames == 0:
        print("No frames to render!")
        return

    # Setup plotter for offscreen rendering
    plotter = pv.Plotter(off_screen=True, window_size=[1920, 1080])
    plotter.set_background("white")

    # Get bounds across the whole simulation so the camera frames the
    # eventual position of the cube (including impact with the floor),
    # not just its initial location.
    all_points = np.concatenate(points_per_frame, axis=0)
    mins = all_points.min(axis=0)
    maxs = all_points.max(axis=0)
    center = (mins + maxs) / 2.0
    extent = maxs - mins
    max_extent = max(extent)

    if camera == "iso":
        # 3/4 view: camera above and to one side, looking down so the floor
        # plane (low Y) and the descending cube are both visible.
        plotter.camera_position = [
            (center[0] + max_extent * 1.4,
             maxs[1] + max_extent * 0.6,
             center[2] + max_extent * 1.4),
            center,
            (0, 1, 0),
        ]
    else:
        # Original side view, framed on the simulation bounds.
        plotter.camera_position = [
            (center[0] + max_extent * 1.5,
             center[1] + max_extent * 0.5,
             center[2]),
            center,
            (0, 1, 0),
        ]

    # Open movie file
    print(f"Rendering to {output_file}...")
    plotter.open_movie(output_file, framerate=fps, quality=8)

    # When boundary particles are included they outnumber the body 50:1
    # (~17.5K vs ~340 for demo1), and rendering them as spheres tanks
    # offscreen rendering on macOS. Render boundary as plain points, and
    # keep only the lower portion (floor + bottom of walls) so the impact
    # surface is visible without obscuring the camera path.
    is_boundary_per_frame = [t >= 3 for t in types_per_frame]
    keep_floor_per_frame = []
    if include_boundary:
        floor_cutoff = mins[1] + 0.15 * (maxs[1] - mins[1])
        for pts, mask in zip(points_per_frame, is_boundary_per_frame):
            keep_floor_per_frame.append(mask & (pts[:, 1] <= floor_cutoff))

        boundary_mesh = pv.PolyData(points_per_frame[0][keep_floor_per_frame[0]])
        plotter.add_mesh(
            boundary_mesh,
            render_points_as_spheres=False,  # spheres × 17K = hang
            color="lightgray",
            opacity=0.6,
            point_size=2,
            show_scalar_bar=False,
        )
    else:
        boundary_mesh = None

    body_pts0 = points_per_frame[0][~is_boundary_per_frame[0]] if include_boundary else points_per_frame[0]
    body_types0 = types_per_frame[0][~is_boundary_per_frame[0]] if include_boundary else types_per_frame[0]
    mesh = pv.PolyData(body_pts0)
    mesh["types"] = body_types0

    actor = plotter.add_mesh(
        mesh,
        render_points_as_spheres=True,
        scalars="types",
        cmap=["blue", "green", "turquoise", "lightgray"],
        point_size=8,
        show_scalar_bar=False
    )

    # Add time annotation
    time_text = plotter.add_text(
        "Time: 0.000 s",
        position="upper_left",
        font_size=14,
        color="black"
    )

    # Render each frame
    for i, (pts, tps) in enumerate(zip(points_per_frame, types_per_frame)):
        if i % 10 == 0:
            pct = 100 * i / num_frames
            print(f"  Frame {i}/{num_frames} ({pct:.1f}%)")

        if include_boundary:
            mask = is_boundary_per_frame[i]
            boundary_mesh.points = pts[keep_floor_per_frame[i]]
            mesh.points = pts[~mask]
            mesh["types"] = tps[~mask]
        else:
            mesh.points = pts
            mesh["types"] = tps

        # Update time text
        sim_time = i * log_step * time_step
        plotter.remove_actor(time_text)
        time_text = plotter.add_text(
            f"Time: {sim_time:.3f} s",
            position="upper_left",
            font_size=14,
            color="black"
        )

        plotter.write_frame()

    plotter.close()
    print(f"Done! Saved to {output_file}")


if __name__ == "__main__":
    position_file = "buffers/position_buffer.txt"
    output_file = "worm_simulation.mp4"

    if len(sys.argv) >= 2:
        position_file = sys.argv[1]
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]

    render_to_movie(position_file, output_file)
