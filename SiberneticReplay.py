"""
A PyVista based viewer/player for saved Sibernetic simulations

Loads in the generated position_buffer.txt file

"""

import pyvista as pv
import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt

from SibSimulation import SibSimulation

from enum import Enum

last_meshes = {}
last_actors = {}

show_liquid = True
show_elastic = True
show_boundary_particles = True
show_footprint_outline = True

replay_speed = 0.02  # seconds between frames
replaying = False

sim_positions = None

plotter = None
offset3d_ = (0, 0, 0)
slider = None

show_boundary = False
max_time = None

verbose = False

musc_chart = None
curv_chart = None

downsample = 1  # only load every nth time point of 3d positions


def print_(msg):
    prefix = "SibReplay: "
    print(prefix + str(msg).replace("\n", "\n" + prefix))


class State(Enum):
    PAUSED = "Paused"
    RUNNING = "Running"


class ReplayController:
    slider_view = None

    def __init__(self, times=None, show_vtk_mesh=False):
        self.times = list(times)
        self.state = State.PAUSED
        self.current_time_index = 0
        self.show_vtk_mesh = show_vtk_mesh
        self.vtk_actor = None

    def play(self, should_play, step=1):
        if should_play:
            print_(" > Starting replay playback.")

            if self.current_time_index == len(self.times) - 1:
                print_(" > Replay at end of time, resetting to start.")
                self.set_to_time(0)

            self.state = State.RUNNING
            while self.state == State.RUNNING and self.current_time_index + 1 < len(
                self.times
            ):
                self.current_time_index = min(
                    self.current_time_index + step, len(self.times) - 1
                )
                print_(f" > Advancing to time index: {self.current_time_index}")
                self.render_all()
                time.sleep(replay_speed)
            print_(" > Replay playback finished or paused.")
            self.state = State.PAUSED

        else:
            print_(" > Pausing replay playback.")
            self.state = State.PAUSED

    def step_forward(self):
        if self.state == State.RUNNING:
            self.state = State.PAUSED

        if self.current_time_index + 1 >= len(self.times):
            print_(" > Replay at end of time, cannot step forward.")
            return

        self.current_time_index += 1
        print_(f" > Stepping forward one time step to index: {self.current_time_index}")
        self.render_all()

    def step_backward(self):
        if self.state == State.RUNNING:
            self.state = State.PAUSED

        if self.current_time_index == 0:
            print_(" > Replay at start of time, cannot step backward.")
            return

        self.current_time_index -= 1
        print_(
            f" > Stepping backward one time step to index: {self.current_time_index}"
        )
        self.render_all()

    def set_to_time(self, time_value):
        if time_value == 0:
            closest_index = 0
            closest_time = 0
        else:
            if time_value in self.times:
                closest_index = self.times.index(time_value)
                closest_time = self.times[closest_index]
                print_(
                    f" > .Finding closest time to {time_value}, got index: {closest_index}"
                )
            else:
                closest = min(self.times, key=lambda x: abs(x - time_value))
                closest_index = self.times.index(closest)
                print_(
                    f" > Finding closest time to {time_value}, got index: {closest_index} of {len(self.times)} times {self.times[0]}-{self.times[-1]}"
                )
                closest_time = self.times[closest_index]

        self.current_time_index = closest_index

        print_(
            " > Replay requested to be set to: %s; being set to time value %s (index: %d)."
            % (time_value, self.current_time_index, closest_time)
        )
        self.state = State.PAUSED
        self.render_all()

    def render_all(self):
        print_(
            "\n > Rendering replay at time index: %d (time: %s)"
            % (self.current_time_index, self.times[self.current_time_index])
        )
        if self.slider_view is not None:
            self.slider_view.GetSliderRepresentation().SetValue(
                self.times[self.current_time_index]
            )
        create_mesh(self.current_time_index)

        if len(sim_positions.vtp_files) > 0 and self.show_vtk_mesh:
            if self.vtk_actor is not None:
                plotter.remove_actor(self.vtk_actor)

            vtk_mesh = pv.read(sim_positions.vtp_files[self.current_time_index])

            self.vtk_actor = plotter.add_mesh(vtk_mesh, color="grey", style="wireframe")

        plotter.render()
        try:
            plotter.update()
        except Exception:
            pass  # plotter may not be initialised yet

    def get_state(self):
        return f" > Replay state: {self.state}, current time index: {self.current_time_index}, so time is {self.times[self.current_time_index]} of max time {self.times[-1]} ({len(self.times)} time points)."


replay_controller = None


def get_color_info_for_type(type_):
    """
    Get color, info string and point size for a given point type
    returns: color, info, size
    """

    if type_ == 1.1:
        return "cyan", "liquid 1", 2
    elif type_ == 1.2:
        return "darkturquoise", "liquid 2", 3

    elif type_ == 2.1:
        return "pink", "elastic 1", 5
    elif type_ == 2.2:
        return "#FF0000", "elastic 2", 5
    elif type_ >= 2.3 and type_ < 2.4:
        return "lightyellow", "elastic variable", 5
    elif type_ > 2 and type_ < 3:
        return "#00cc00", "elastic variable", 7

    elif type_ == 3:
        return "grey", "boundary 0", 3
    elif type_ == 3.1:
        return "black", "boundary 1", 7
    else:
        return "orange", "unknown", 5


def add_muscle_activation_chart(sim_dir, pl, duration=None):
    global musc_chart

    muscle_activation_file = os.path.join(sim_dir, "muscles_activity_buffer.txt")
    print_(f"Loading muscle activation file from: {muscle_activation_file}")
    musc_dat = np.loadtxt(muscle_activation_file, delimiter="\t").T

    if np.sum(musc_dat) == 0:
        print_("No muscle data, not plotting")
        return

    # print(musc_dat.shape)
    # plt.imshow(musc_dat, interpolation="none", aspect="auto", cmap="YlOrRd")

    f_musc, ax_musc = plt.subplots(tight_layout=True)

    im = ax_musc.imshow(musc_dat, interpolation="none", aspect="auto", cmap="YlOrRd")
    f_musc.canvas.manager.set_window_title("Muscle Activation Heatmap")

    f_musc.colorbar(im)

    if duration is not None:
        num_ticks = 5
        ax_musc.set_xticks(np.linspace(0, musc_dat.shape[1], num_ticks))
        ax_musc.set_xticklabels(np.linspace(0, duration, num_ticks))

        ax_musc.set_xlabel("Time (ms)")
    else:
        ax_musc.set_xlabel("Time point")

    _ = ax_musc.set_ylabel("Muscle")

    musc_chart = pv.ChartMPL(f_musc, size=(0.35, 0.35), loc=(0.02, 0.06))
    musc_chart.title = None
    musc_chart.border_color = "white"
    musc_chart.background_color = (1.0, 1.0, 1.0, 0.4)

    pl.add_chart(
        musc_chart,
    )


def add_body_curvature_chart(sim_dir, pl, duration=None):
    global curv_chart

    from wcon.generate_wcon import generate_wcon

    sib_position_file = os.path.join(sim_dir, "worm_motion_log.txt")
    wcon_output_file = "/tmp/worm_motion_log.wcon"

    x, y, z, ts, body_curv_data = generate_wcon(
        sib_position_file,
        wcon_file_name=wcon_output_file,
        rate_to_plot=1,
        plot=False,
        verbose=False,
    )

    if np.sum(x[ts[-1]]) + np.sum(y[ts[-1]]) == 0:
        print_("No worm body data, not plotting")
        return

    print_(f"Temporary WCON file (re)generated at: {wcon_output_file}")

    f_curv, ax_curv = plt.subplots(tight_layout=True)
    im = ax_curv.imshow(
        body_curv_data.transpose(),
        interpolation="none",
        aspect="auto",
        cmap="bwr",
        vmin=170,
        vmax=190,
    )
    f_curv.colorbar(im)
    f_curv.canvas.manager.set_window_title("Body Curvature")

    if duration is not None:
        ax_curv.set_xlabel("Time (ms)")

        num_ticks = 5
        ax_curv.set_xticks(np.linspace(0, body_curv_data.shape[0], num_ticks))
        ax_curv.set_xticklabels(np.linspace(0, duration, num_ticks))
    else:
        ax_curv.set_xlabel("Time point")

    _ = ax_curv.set_ylabel("Body curv.")
    ax_curv.set_in_layout(False)

    curv_chart = pv.ChartMPL(f_curv, size=(0.35, 0.35), loc=(0.62, 0.06))
    curv_chart.title = None
    curv_chart.border_color = "white"
    curv_chart.background_color = (1.0, 1.0, 1.0, 0.4)
    pl.add_chart(
        curv_chart,
    )


def add_sibernetic_model(
    pl,
    position_file="Sibernetic/position_buffer.txt",
    report_file=None,
    swap_y_z=False,
    offset3d=(0, 0, 0),
    include_boundary=False,
    show_footprint=True,
):
    global \
        sim_positions, \
        last_meshes, \
        plotter, \
        offset3d_, \
        slider, \
        show_boundary, \
        show_boundary_particles, \
        max_time, \
        replay_controller, \
        show_footprint_outline

    print_(
        f"Adding Sibernetic model from position file: {position_file}, report file: {report_file}"
    )
    print_(
        f"Swap Y/Z: {swap_y_z}, offset3d: {offset3d}, include_boundary: {include_boundary}  "
    )

    offset3d_ = offset3d
    plotter = pl
    show_boundary = include_boundary
    show_boundary_particles = include_boundary
    show_footprint_outline = show_footprint

    sim_positions = SibSimulation(
        position_file=position_file,
        report_file=report_file,
        downsample=downsample,
        swap_y_z=swap_y_z,
    )

    if sim_positions.has_worm_data():
        add_muscle_activation_chart(sim_positions.sim_dir, pl, sim_positions.duration)
        add_body_curvature_chart(sim_positions.sim_dir, pl, sim_positions.duration)
    elif (
        sim_positions.report_data is None and sim_positions.has_muscle_activation_data()
    ):
        print_(
            f"Found muscle activation file at: {sim_positions.sim_dir}/muscles_activity_buffer.txt, adding muscle activation chart."
        )
        add_muscle_activation_chart(sim_positions.sim_dir, pl)
        add_body_curvature_chart(sim_positions.sim_dir, pl)
    else:
        print_(
            "No report file provided and no muscle activation file found, skipping muscle activation chart."
        )

    if replay_controller is None:
        replay_controller = ReplayController(times=sim_positions.loaded_time_points)

    create_mesh(0)

    slider_text = "Time point"

    if sim_positions.duration is not None:
        max_time = sim_positions.duration
        slider_text = "Time (ms)"
    else:
        max_time = sim_positions.num_time_points() - 1

    slider = pl.add_slider_widget(
        slider_updated, rng=[0, max_time], value=0, title=slider_text, style="modern"
    )

    replay_controller.slider_view = slider

    button_height = 10
    button_separation = 70
    txt_offset = 8
    txt_voffset = 12

    b0 = 10
    pl.add_checkbox_button_widget(
        info_checkbox_pressed,
        value=False,
        position=(b0, button_height),
        color_on="lightgrey",
        color_off="lightgrey",
    )
    pl.add_text(
        "  i",
        position=(b0 + txt_offset, button_height + txt_voffset),
        font_size=12,
        color="black",
    )

    b1 = b0 + button_separation
    pl.add_checkbox_button_widget(
        back_checkbox_pressed,
        value=False,
        position=(b1, button_height),
        color_on="lightgrey",
        color_off="lightgrey",
    )
    pl.add_text(
        "<|",
        position=(b1 + txt_offset, button_height + txt_voffset),
        font_size=12,
        color="black",
    )

    b2 = b1 + button_separation
    pl.add_checkbox_button_widget(
        play_checkbox_pressed,
        value=False,
        position=(b2, button_height),
        color_on="lightgrey",
        color_off="darkgrey",
    )
    pl.add_text(
        " >",
        position=(b2 + txt_offset, button_height + txt_voffset),
        font_size=12,
        color="black",
    )

    b3 = b2 + button_separation
    pl.add_checkbox_button_widget(
        ff_checkbox_pressed,
        value=False,
        position=(b3, button_height),
        color_on="lightgrey",
        color_off="darkgrey",
    )
    pl.add_text(
        ">>",
        position=(b3 + txt_offset, button_height + txt_voffset),
        font_size=12,
        color="black",
    )

    b4 = b3 + button_separation
    pl.add_checkbox_button_widget(
        fwd_checkbox_pressed,
        value=False,
        position=(b4, button_height),
        color_on="lightgrey",
        color_off="lightgrey",
    )
    pl.add_text(
        "|>",
        position=(b4 + txt_offset, button_height + txt_voffset),
        font_size=12,
        color="black",
    )

    if len(sim_positions.vtp_files) > 0:
        pl.add_checkbox_button_widget(
            toggle_vtk_mesh,
            value=False,
            position=(1180, button_height),
            color_on="lightblue",
            color_off="darkgrey",
        )

    pl.add_checkbox_button_widget(
        toggle_boundary_particles,
        value=show_boundary,
        position=(1250, button_height),
        color_on="grey",
        color_off="darkgrey",
    )

    pl.add_checkbox_button_widget(
        toggle_elastic,
        value=True,
        position=(1320, button_height),
        color_on="pink",
        color_off="darkgrey",
    )

    pl.add_checkbox_button_widget(
        toggle_liquid,
        value=True,
        position=(1390, button_height),
        color_on="cyan",
        color_off="darkgrey",
    )

    if musc_chart is not None:
        pl.add_checkbox_button_widget(
            toggle_musc_chart,
            value=True,
            position=(1470, button_height),
            color_on="green",
            color_off="darkgrey",
        )

    if curv_chart is not None:
        pl.add_checkbox_button_widget(
            toggle_curv_chart,
            value=True,
            position=(1540, button_height),
            color_on="green",
            color_off="darkgrey",
        )


def slider_updated(value):
    global replay_controller
    print_(
        f" > Slider updated to value: {value}, replay: {replay_controller.get_state()}"
    )

    replay_controller.set_to_time(value)


def info_checkbox_pressed(value):
    global replay_controller
    print_(f" > Info checkbox pressed, value: {value}")
    if value:
        print_(" > Showing sim info:")

        fig, ax = plt.subplots(figsize=(8, 6))

        fig.suptitle("Sibernetic Replay Info")
        info_lines = [
            f"Total time points loaded: {sim_positions.num_time_points()}",
            f"Total points per time point: {len(sim_positions.all_point_types[0])}",
        ]
        if sim_positions.report_data is not None:
            for key, val in sim_positions.report_data.items():
                info_lines.append(f"{key}:  {val}")

        # Remove the axes
        ax.axis("off")

        # Add text to the center of the figure
        ax.text(
            0.5,
            0.5,
            "\n".join(info_lines),
            fontsize=8,
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.8),
        )

        # Display the window
        plt.show()


def fwd_checkbox_pressed(value):
    global replay_controller
    print_(f" > Fwd checkbox pressed, value: {value}")
    replay_controller.step_forward()


def play_checkbox_pressed(value):
    global replay_controller
    print_(f" > Play checkbox pressed, value: {value}")
    replay_controller.play(value, 1)


def ff_checkbox_pressed(value):
    global replay_controller
    print_(f" > FF checkbox pressed, value: {value}")
    replay_controller.play(value, 3)


def back_checkbox_pressed(value):
    global replay_controller
    print_(f" > Back checkbox pressed, value: {value}")
    replay_controller.step_backward()


def toggle_musc_chart(value):
    global musc_chart, plotter
    print_(f" > Muscle activation chart toggle, value: {value}")
    if musc_chart is not None:
        musc_chart.visible = value
        plotter.render()


def toggle_curv_chart(value):
    global curv_chart, plotter
    print_(f" > Body curvature chart toggle, value: {value}")
    if curv_chart is not None:
        curv_chart.visible = value
        plotter.render()


def toggle_liquid(value):
    global show_liquid, last_actors, plotter
    show_liquid = value
    for type_, actor in last_actors.items():
        _, info, _ = get_color_info_for_type(type_)
        if "liquid" in info:
            actor.visibility = value
    plotter.render()


def toggle_elastic(value):
    global show_elastic, last_actors, plotter
    show_elastic = value
    for type_, actor in last_actors.items():
        _, info, _ = get_color_info_for_type(type_)
        if "elastic" in info:
            actor.visibility = value
    plotter.render()


def toggle_boundary_particles(value):
    global show_boundary_particles, last_actors, plotter
    print_(f" > Boundary particles toggle, value: {value}")
    show_boundary_particles = value
    for type_, actor in last_actors.items():
        _, info, _ = get_color_info_for_type(type_)
        if "boundary" in info:
            actor.visibility = value
    plotter.render()


def toggle_vtk_mesh(value):
    global replay_controller, plotter
    print_(f" > VTK mesh toggle, value: {value}")
    replay_controller.show_vtk_mesh = value
    if not value and replay_controller.vtk_actor is not None:
        plotter.remove_actor(replay_controller.vtk_actor)
        replay_controller.vtk_actor = None
        plotter.render()
    elif value:
        replay_controller.render_all()


def create_mesh(time_index):
    global \
        sim_positions, \
        last_meshes, \
        last_actors, \
        plotter, \
        offset3d_, \
        show_boundary, \
        show_liquid, \
        show_elastic, \
        show_boundary_particles

    if time_index >= sim_positions.num_time_points():
        print_(
            "Index %i out of bounds for loaded time points with length %i"
            % (time_index, sim_positions.num_time_points())
        )
        return

    print_(
        "   -- Creating new mesh at time point index: %s/%s"
        % (time_index, sim_positions.num_time_points())
    )
    curr_points_dict = sim_positions.get_points_at(time_index)

    print_("      Plotting %i point types" % (len(curr_points_dict)))

    for type_, curr_points in curr_points_dict.items():
        color, info, size = get_color_info_for_type(type_)

        is_boundary = "boundary" in info

        if is_boundary and time_index == 0 and not show_boundary_particles:
            # draw bounding box outline as a visual placeholder when boundary particles are hidden
            mx = np.max(curr_points, axis=0)
            mn = np.min(curr_points, axis=0)
            swap = False
            if swap:
                a = [mn[0], mn[1], mn[2]]
                b = [mn[0], mx[1], mn[2]]
                c = [mn[0], mx[1], mx[2]]
                d = [mn[0], mn[1], mx[2]]
            else:
                a = [mn[0], mn[1], mn[2]]
                b = [mx[0], mn[1], mn[2]]
                c = [mx[0], mn[1], mx[2]]
                d = [mn[0], mn[1], mx[2]]

            print_(f"        >>>>>>>>>>   Boundary box points: {a}, {b}, {c}, {d}")
            points = np.array([a, b, b, c, c, d, d, a])
            if show_footprint_outline:
                plotter.add_lines(points, color="grey", width=2)
                # fall through so boundary particles are still added to last_actors (hidden)

        if verbose:
            print_(
                "       - Plotting %i points of type '%s' (%s), color: %s, size: %i"
                % (len(curr_points), type_, info, color, size)
            )

        if len(curr_points) == 0:
            continue
        if type_ not in last_meshes:
            last_meshes[type_] = pv.PolyData(curr_points)
            last_meshes[type_].translate(offset3d_, inplace=True)

            actor = plotter.add_mesh(
                last_meshes[type_],
                render_points_as_spheres=True,
                point_size=size,
                color=color,
            )
            last_actors[type_] = actor
            if "liquid" in info and not show_liquid:
                actor.visibility = False
            elif "elastic" in info and not show_elastic:
                actor.visibility = False
            elif "boundary" in info and not show_boundary_particles:
                actor.visibility = False
        else:
            if not is_boundary:
                last_meshes[type_].points = curr_points
                last_meshes[type_].translate(
                    (offset3d_[0], offset3d_[1], offset3d_[2]), inplace=True
                )
            else:
                print_("Boundary points not translated")

    plotter.render()
    # time.sleep(0.1)

    return


if __name__ == "__main__":
    plotter = pv.Plotter()

    default_position_file = "buffers/position_buffer.txt"  # can be overwritten by arg
    report_file = None

    if not os.path.isfile(default_position_file):
        default_position_file = (
            "Sibernetic/position_buffer.txt"  # example location in Worm3DViewer repo
        )

    include_boundary = False

    if "-b" in sys.argv:
        include_boundary = True
    else:
        print_("Run with -b to display boundary box")

    if len(sys.argv) > 1:
        if os.path.isdir(sys.argv[1]):
            if os.path.isfile(os.path.join(sys.argv[1], "report.json")):
                position_file = None
                report_file = os.path.join(sys.argv[1], "report.json")
                dir_name = os.path.dirname(report_file)

            elif os.path.isfile(os.path.join(sys.argv[1], "position_buffer.txt")):
                position_file = os.path.join(sys.argv[1], "position_buffer.txt")
                dir_name = os.path.dirname(position_file)
            else:
                raise ValueError(
                    f"Provided argument is a directory but no report.json or position_buffer.txt file found in it: {sys.argv[1]}"
                )

        elif "json" in sys.argv[1] and os.path.isfile(sys.argv[1]):
            position_file = None
            report_file = sys.argv[1]
            dir_name = os.path.dirname(report_file)

        else:
            if not os.path.isfile(sys.argv[1]):
                raise ValueError(
                    f"Provided argument is not a valid file or directory: {sys.argv[1]}"
                )
            position_file = sys.argv[1]
            dir_name = os.path.dirname(position_file)

    else:
        position_file = default_position_file

    swap_y_z = False

    add_sibernetic_model(
        plotter,
        position_file,
        report_file,
        swap_y_z=swap_y_z,
        include_boundary=include_boundary,
    )
    plotter.window_size = [1600, 800]

    plotter.set_background("white")
    plotter.add_axes()

    if swap_y_z:
        plotter.camera_position = "zx"
        plotter.camera.roll = 90
        plotter.camera.elevation = 45
    else:
        plotter.camera_position = "yz"
        plotter.camera.roll = 0
        plotter.camera.elevation = 25

    # print(plotter.camera_position)

    def on_close_callback(plotter):
        global replay_controller
        print_(
            f"Plotter window is closing. Performing actions now (replay: {replay_controller.get_state()})."
        )
        replay_controller.state = State.PAUSED

    if "-nogui" not in sys.argv:
        plotter.show(before_close_callback=on_close_callback, auto_close=True)
        print_("Done showing")
