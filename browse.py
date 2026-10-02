from nicegui import ui
import sys
import os
from SibSimulation import SibSimulation
from SiberneticReplay import add_sibernetic_model

import pyvista as pv


def main(sims_dir=None):

    simulations = []
    for f in sorted(os.listdir(sims_dir)):
        ff = os.path.join(sims_dir, f)
        report = os.path.join(ff, "report.json")
        if os.path.isdir(ff) and os.path.isfile(report):
            print("--- Adding sim dir: " + ff)
            ss = SibSimulation(report_file=report, load_positions=False)
            simulations.append(ss)

    rows = []
    for s in simulations:
        print(s)
        duration = float(s.report_data.get("duration").split()[0])
        run_time = round(float(s.report_data.get("run_time").split()[0]), 5)
        rows.append(
            {
                "name": s.report_data.get("sim_ref", "???"),
                "duration": duration,
                "run_time": run_time,
                "report_file": s.report_file,
            }
        )

    table = ui.table(
        columns=[
            {"name": "name", "label": "Name", "field": "name", "sortable": True},
            {
                "name": "duration",
                "label": "Duration (ms)",
                "field": "duration",
                "sortable": True,
            },
            {
                "name": "run_time",
                "label": "Run Time (s)",
                "field": "run_time",
                "sortable": True,
            },
            {"name": "action", "label": "Replay 3D", "align": "center"},
            {"name": "action2", "label": "Replay 2D", "align": "center"},
        ],
        rows=rows,
        row_key="name",
    )

    table.add_slot(
        "body-cell-action",
        """
        <q-td :props="props">
            <q-btn flat label="Load" @click="() => $parent.$emit('load', props.row.report_file)" />
        </q-td>
    """,
    )
    table.add_slot(
        "body-cell-action2",
        """
        <q-td :props="props">
            <q-btn flat label="Load" @click="() => $parent.$emit('load2', props.row.report_file)" />
        </q-td>
    """,
    )
    table.on("load", lambda e: load_sim(e))
    table.on("load2", lambda e: load_sim_2d(e))

    ui.run(native=True)  # remove native=True to serve it as a web app


def load_sim(e):
    sim_name = e.args
    print(f"Loading simulation: {sim_name}")

    plotter = pv.Plotter()
    swap_y_z = False
    add_sibernetic_model(
        plotter,
        position_file=None,
        report_file=sim_name,
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
    sim_name = e.args
    print(f"Loading 2D simulation: {sim_name}")

    # Here you can add code to load the 2D simulation based on sim_name


if __name__ in {"__main__", "__mp_main__"}:
    if len(sys.argv) < 2:
        print("Usage: python browse.py <sims_dir>")
        sys.exit(1)

    sims_dir = sys.argv[1]
    main(sims_dir)
