from nicegui import ui
import sys
import os
from SibSimulation import SibSimulation
from W2DSimulation import W2DSimulation
from SiberneticReplay import add_sibernetic_model

import pyvista as pv

simulations = []


def load_simulation(sim_dir, full_load=False):

    sibernetic_report = os.path.join(sim_dir, "report.json")

    if os.path.isdir(sim_dir) and os.path.isfile(sibernetic_report):
        print("--- Adding sim dir: " + sim_dir)
        ss = SibSimulation(
            report_file=sibernetic_report, load_positions=False, load_traces=full_load
        )
        return ss

    worm2d_report = os.path.join(sim_dir, "simulation_pars.json")

    if os.path.isdir(sim_dir) and os.path.isfile(worm2d_report):
        print("--- Adding sim dir: " + sim_dir)
        ss = W2DSimulation(sim_dir=sim_dir, verbose=True, load_traces=full_load)
        return ss

    print(f"  - Directory {sim_dir} does not contain valid simulation files.")
    return None


def main(sims_dir=None):

    for f in sorted(os.listdir(sims_dir)):
        sim_dir = os.path.join(sims_dir, f)
        if os.path.isdir(sim_dir):
            sim = load_simulation(sim_dir)
            if sim is not None:
                simulations.append(sim)

    if len(simulations) == 0:
        print("No simulations found in directory: " + sims_dir)
        sys.exit(1)

    rows = []
    for s in simulations:
        print(s)
        duration = s.get_duration_ms()
        run_time = s.get_run_time_s()
        timestamp = s.get_timestamp()
        sim_dir_name = os.path.basename(s.sim_dir) if hasattr(s, "sim_dir") else "???"
        rows.append(
            {
                "name": sim_dir_name,
                "timestamp": timestamp,
                "duration": duration,
                "run_time": run_time,
                "sim_dir": s.sim_dir,
            }
        )

    # set the title of the nicegui window to "Sibernetic Simulation Browser"
    ui.label("OpenWorm Simulation Browser").style(
        "font-size: 18px; font-weight: bold; margin-bottom: 20px;"
    )

    table = ui.table(
        columns=[
            {"name": "name", "label": "Name", "field": "name", "sortable": True},
            {
                "name": "timestamp",
                "label": "Timestamp",
                "field": "timestamp",
                "sortable": True,
                # sort on the raw epoch value, but display it as a readable date/time
                ":format": "val => new Date(val * 1000).toLocaleString()",
            },
            {
                "name": "duration",
                "label": "Duration (ms)",
                "field": "duration",
                "sortable": True,
            },
            {
                "name": "run_time",
                "label": "Run Time",
                "field": "run_time",
                "sortable": True,
                # if the time in seconds is greater than 60, display it in minutes, otherwise in seconds, or if greater than 3600, display in hours
                ":format": "val => val > 3600 ? (val / 3600).toFixed(2) + ' hr' : val > 60 ? (val / 60).toFixed(2) + ' min' : val.toFixed(2) + ' s'",
            },
            {"name": "action", "label": "Replay 3D", "align": "center"},
            {"name": "action2", "label": "Replay 2D", "align": "center"},
        ],
        rows=rows,
        row_key="name",
        # newest first, so the most recently run simulation is in the top row
        pagination={"sortBy": "timestamp", "descending": True, "rowsPerPage": 25},
    )

    table.add_slot(
        "body-cell-action",
        """
        <q-td :props="props">
            <q-btn flat label="Load" @click="() => $parent.$emit('load', props.row.sim_dir)" />
        </q-td>
    """,
    )
    table.add_slot(
        "body-cell-action2",
        """
        <q-td :props="props">
            <q-btn flat label="Load" @click="() => $parent.$emit('load2', props.row.sim_dir)" />
        </q-td>
    """,
    )
    table.on("load", lambda e: load_sim(e))
    table.on("load2", lambda e: load_sim_2d(e))

    ui.run(
        native=True, window_size=(900, 600)
    )  # remove native=True to serve it as a web app


def load_sim(e):
    sim_name = e.args
    report_file = os.path.join(sim_name, "report.json")
    print(f"Loading simulation: {sim_name}, with report file: {report_file}")

    plotter = pv.Plotter()
    swap_y_z = False
    add_sibernetic_model(
        plotter,
        position_file=None,
        report_file=report_file,
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
        print("Closing...")
        """
        global replay_controller
        print(
            f"Plotter window is closing. Performing actions now (replay: {replay_controller.get_state()})."
        )
        replay_controller.state = State.PAUSED"""

    if "-nogui" not in sys.argv:
        plotter.show(before_close_callback=on_close_callback, auto_close=True)
        print("Done showing")

    # Here you can add code to load the simulation based on sim_name


def load_sim_2d(e):

    from wconviewer.WormView import show_worm_view

    sim_dir = e.args
    wcon_file = os.path.join(sim_dir, "worm_motion_log.wcon")
    if not os.path.isfile(wcon_file):
        wcon_file = os.path.join(sim_dir, "output.wcon")

    sim = load_simulation(sim_dir, full_load=True)

    sim.display_traces(show=False)

    print(f"Loading 2D simulation: {sim_dir} with WCON file: {wcon_file}")

    show_worm_view(
        wcon_file,
        show_head=False,
        zoom_to_worm=False,
        show_grid=False,
        nogui=False,
        start_play=False,
    )


if __name__ in {"__main__", "__mp_main__"}:
    if len(sys.argv) < 2:
        print("Usage: python browse.py <sims_dir>")
        sys.exit(1)

    sims_dir = sys.argv[1]
    main(sims_dir)
