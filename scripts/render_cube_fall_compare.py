"""Side-by-side OpenCL vs Metal cube-drop comparison MP4.

Reads both position_buffer.txt files, samples N evenly-spaced frames at
matched sim times, renders each frame as 3D scatter plot, combines into
a single MP4 with both panels visible.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter


def parse_position_buffer(path):
    """Return: (n_active, n_boundary, dt, log_step, frames),
    where frames is list of (active_xyz [n_active, 3], type_codes)."""
    with open(path) as f:
        header = [f.readline().rstrip('\n') for _ in range(11)]
        n_elastic = int(header[6])
        n_liquid = int(header[7])
        n_boundary = int(header[8])
        dt = float(header[9])
        log_step = int(header[10])
        n_active = n_elastic + n_liquid

        # First frame: all particles
        boundary_xyz = []
        active_xyz_first = []
        active_types_first = []
        for _ in range(n_active + n_boundary):
            ws = f.readline().split()
            x, y, z, t = float(ws[0]), float(ws[1]), float(ws[2]), float(ws[3])
            if t >= 3.0:
                boundary_xyz.append([x, y, z])
            else:
                active_xyz_first.append([x, y, z])
                active_types_first.append(t)
        boundary_xyz = np.array(boundary_xyz)
        active_xyz_first = np.array(active_xyz_first)
        active_types_first = np.array(active_types_first)

        frames = [(active_xyz_first, active_types_first)]

        # Subsequent frames: each may have just active particles (Sibernetic)
        # or active+boundary (our normalized Metal dump).
        # We detect which by reading lines and looking for boundary types.
        # Simpler: read line at a time, batch into N=n_active or N=n_active+n_boundary frames.
        line_buf = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            line_buf.append(line)

        # Heuristic: try interpreting as N=n_active per frame first.
        if len(line_buf) % n_active == 0:
            chunk = n_active
        elif len(line_buf) % (n_active + n_boundary) == 0:
            chunk = n_active + n_boundary
        else:
            chunk = n_active  # default

        for i in range(0, len(line_buf), chunk):
            chunk_lines = line_buf[i:i+chunk]
            xyz = []
            ts = []
            for ln in chunk_lines:
                ws = ln.split()
                t = float(ws[3])
                if t < 3.0:
                    xyz.append([float(ws[0]), float(ws[1]), float(ws[2])])
                    ts.append(t)
            if len(xyz) == n_active:
                frames.append((np.array(xyz), np.array(ts)))

    return n_active, n_boundary, dt, log_step, boundary_xyz, frames


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else "/Users/slarson/Downloads/cube_drop_compare.mp4"
    left_path = sys.argv[2] if len(sys.argv) > 2 else "/tmp/sib_centered_dump/position_buffer.txt"
    right_path = sys.argv[3] if len(sys.argv) > 3 else "/tmp/metal_centered_nopair.txt"
    left_label = sys.argv[4] if len(sys.argv) > 4 else "Reference"
    right_label = sys.argv[5] if len(sys.argv) > 5 else "Metal Native Substrate (no pair forces)"

    print(f"Loading LEFT: {left_path}")
    o_na, o_nb, o_dt, o_ls, o_boundary, o_frames = parse_position_buffer(left_path)
    print(f"  {len(o_frames)} frames, dt={o_dt}, log_step={o_ls}")
    print(f"  sim_t range: 0 to {(len(o_frames)-1)*o_ls*o_dt:.3f}s")

    print(f"Loading RIGHT: {right_path}")
    m_na, m_nb, m_dt, m_ls, m_boundary, m_frames = parse_position_buffer(right_path)
    print(f"  {len(m_frames)} frames, dt={m_dt}, log_step={m_ls}")
    print(f"  sim_t range: 0 to {(len(m_frames)-1)*m_ls*m_dt:.3f}s")

    # Sample N matched sim times
    N = 60  # frames in output
    target_t_max = min((len(o_frames)-1)*o_ls*o_dt, (len(m_frames)-1)*m_ls*m_dt)
    print(f"Sampling {N} frames over 0..{target_t_max:.3f}s sim time")

    o_dt_per_frame = o_ls * o_dt
    m_dt_per_frame = m_ls * m_dt
    target_times = np.linspace(0, target_t_max, N)
    o_indices = [min(int(round(t / o_dt_per_frame)), len(o_frames)-1) for t in target_times]
    m_indices = [min(int(round(t / m_dt_per_frame)), len(m_frames)-1) for t in target_times]

    # Setup figure
    fig = plt.figure(figsize=(14, 7))
    ax_o = fig.add_subplot(1, 2, 1, projection='3d')
    ax_m = fig.add_subplot(1, 2, 2, projection='3d')

    bbox = [0, 88.84]  # demo1 box

    def setup_axis(ax, title):
        ax.set_xlim(bbox); ax.set_ylim(bbox); ax.set_zlim(bbox)
        ax.set_box_aspect([1, 1, 1])
        ax.set_title(title, fontsize=14)
        ax.view_init(elev=15, azim=-60)
        ax.set_xlabel('x'); ax.set_ylabel('z'); ax.set_zlabel('y')
        ax.grid(False)
        # Boundary as thin scatter at floor
        # Hide box outline for cleaner view
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

    setup_axis(ax_o, left_label)
    setup_axis(ax_m, right_label)

    # Plot boundary floor (subset of boundary at low y) for spatial reference
    floor_o = o_boundary[o_boundary[:, 1] < 5]  # bottom 5 units
    floor_m = m_boundary[m_boundary[:, 1] < 5]
    ax_o.scatter(floor_o[:, 0], floor_o[:, 2], floor_o[:, 1],
                  s=0.3, c='lightgray', alpha=0.4)
    ax_m.scatter(floor_m[:, 0], floor_m[:, 2], floor_m[:, 1],
                  s=0.3, c='lightgray', alpha=0.4)

    # Active particle scatters (initially empty)
    o_xyz, o_ts = o_frames[0]
    m_xyz, m_ts = m_frames[0]
    o_colors = ['blue' if t < 2.0 else 'green' for t in o_ts]
    m_colors = ['blue' if t < 2.0 else 'green' for t in m_ts]

    o_scatter = ax_o.scatter(o_xyz[:, 0], o_xyz[:, 2], o_xyz[:, 1],
                              s=8, c=o_colors)
    m_scatter = ax_m.scatter(m_xyz[:, 0], m_xyz[:, 2], m_xyz[:, 1],
                              s=8, c=m_colors)
    time_text = fig.suptitle("t = 0.000 s", fontsize=16)

    def update(frame_idx):
        oi = o_indices[frame_idx]
        mi = m_indices[frame_idx]
        sim_t = target_times[frame_idx]
        o_xyz, _ = o_frames[oi]
        m_xyz, _ = m_frames[mi]
        o_scatter._offsets3d = (o_xyz[:, 0], o_xyz[:, 2], o_xyz[:, 1])
        m_scatter._offsets3d = (m_xyz[:, 0], m_xyz[:, 2], m_xyz[:, 1])
        time_text.set_text(f"t = {sim_t:.3f} s")
        return o_scatter, m_scatter, time_text

    print(f"Rendering {N} frames to {out_path}")
    writer = FFMpegWriter(fps=10, codec='libx264', bitrate=2000)
    anim = FuncAnimation(fig, update, frames=N, blit=False)
    anim.save(out_path, writer=writer, dpi=120)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
