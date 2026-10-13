"""
Manages all generated output files from a Worm2D simulation run.
"""

import os
import json
import numpy as np
import sys

from SibSimulation import plot_muscle_activity, plot_neuron_activity

import matplotlib.pyplot as plt


def print_(msg, print_it=True):
    prefix = "W2DSim: "
    if print_it:
        print(prefix + str(msg).replace("\n", "\n" + prefix))


class W2DSimulation:
    """
    Loads and manages all output files from a Worm2D simulation run.

    """

    def __init__(
        self,
        sim_dir,
        load_traces=True,
        verbose=True,
    ):
        """
        Parameters
        ----------
        ...
        """

        self.sim_dir = sim_dir
        self.verbose = verbose

        self.sim_json_file = os.path.join(self.sim_dir, "simulation_pars.json")

        self._load_json_files()

        if load_traces:
            self._load_traces()

    def __repr__(self):
        return f"W2DSimulation(sim_dir={self.sim_dir}, duration={self.duration}, transient={self.transient})"

    def get_duration_ms(self):
        return self.duration

    def get_run_time_s(self):
        return -1

    def get_timestamp(self):
        # when the simulation was run based on last modified date of sim_json_file
        timestamp = os.path.getmtime(self.sim_json_file)
        return timestamp

    def _load_traces(self):

        self.muscle_states = {}
        self.neuron_states = {}

        act_file = os.path.join(self.sim_dir, "act.dat")

        if os.path.isfile(act_file):
            ma = np.loadtxt(act_file)
            print_(f"Loaded activity from {act_file}: shape {ma.shape}")

            self.neuron_states["t"] = ma[:, 0]
            self.muscle_states["t"] = ma[:, 0]

        cell_ids_file = os.path.join(self.sim_dir, "cell_Ids.json")

        if os.path.isfile(cell_ids_file):
            cell_ids = json.load(open(cell_ids_file))
            infos = cell_ids["Nervous System"]
            print_(f"Loaded {len(infos)} cell ID infos from {cell_ids_file}")

            for i, cell_id in enumerate(infos):
                cell_index = cell_id["Ind"]
                cell_type = cell_id["Cell"]
                cell_name = f"{cell_type}{cell_index}"
                print(f"Handling cell {i}: {cell_name}")
                self.neuron_states[cell_name] = ma[:, i + 1]

            for i in range(24):
                for dv in ["D", "V"]:
                    muscle_name = f"M{dv}L{i + 1}"
                    offset = 0 if dv == "D" else 1
                    self.muscle_states[muscle_name] = ma[
                        :, len(infos) + i * 2 + 1 + offset
                    ]
        else:
            print_(
                f"Cell IDs file not found: {cell_ids_file}. Cannot load neuron and muscle states."
            )
            self.neuron_states = None
            self.muscle_states = None

    def _load_json_files(self):

        self.sim_data = json.load(open(self.sim_json_file))
        print_(f"Loaded simulation data:\n{self.sim_data}")
        self.duration = float(self.sim_data["Duration"])
        self.transient = float(self.sim_data["Transient"])

    def display_traces(self, show=True):

        if self.neuron_states != None:
            plot_neuron_activity(
                self.neuron_states,
                label="",
                title="Neuronal activation (Worm2D)",
            )

        if self.muscle_states != None:
            plot_muscle_activity(
                self.muscle_states,
                label="",
                title="Muscle activation (Worm2D)",
            )

        if show:
            plt.show()


if __name__ == "__main__":
    filename = sys.argv[1] if len(sys.argv) > 1 else "."

    ss = W2DSimulation(sim_dir=filename, verbose=True)

    print("Loaded simulation: %s" % (ss))

    ss.display_traces()
