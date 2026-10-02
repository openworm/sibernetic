"""
Manages all generated output files from a Sibernetic simulation run.
Handles reading position_buffer.txt, report.json, VTK files, and associated
simulation data, independent of any GUI or visualisation layer.
"""

import os
import json
import numpy as np


def print_(msg, print_it=True):
    prefix = "SibSim: "
    if print_it:
        print(prefix + str(msg).replace("\n", "\n" + prefix))


class SibSimulation:
    """
    Loads and manages all output files from a Sibernetic simulation run.

    Can be initialised from a position_buffer.txt file directly, or from a
    report.json file (which also provides timing metadata and locates the
    position file automatically).

    Attributes
    ----------
    all_3D_points : list of dict
        One entry per loaded time point. Each dict maps particle type (float)
        to a list of [x, y, z] positions.
    all_point_types : list of list
        One entry per loaded time point. Each inner list contains the type
        float for every particle in that frame (same order as the file).
    loaded_time_points : list
        The time value (ms if dt is known, otherwise frame index) for each
        loaded frame.
    vtp_files : list of str
        Sorted paths to any .vtp files found alongside the position file.
    report_data : dict or None
        Parsed report.json content, if a report file was provided.
    sim_dir : str
        Directory containing the simulation output files.
    """

    def __init__(
        self,
        position_file=None,
        report_file=None,
        downsample=1,
        swap_y_z=False,
        verbose=False,
    ):
        """
        Parameters
        ----------
        position_file : str, optional
            Path to position_buffer.txt. Inferred from report_file directory
            if not given.
        report_file : str, optional
            Path to report.json. Provides dt/duration/logstep metadata and
            locates position_file automatically.
        downsample : int
            Load every nth time point (1 = load all).
        swap_y_z : bool
            If True, swap the x and y coordinates on load (matches the
            swap_y_z behaviour in SiberneticReplay).
        """
        self.downsample = downsample
        self.swap_y_z = swap_y_z

        self.report_data = None
        self.dt = None
        self.duration = None
        self.log_step = None
        self.sim_dir = None
        self.time_step = None

        self.num_elastic = 0
        self.num_liquid = 0
        self.num_boundary = 0

        self.all_3D_points = []
        self.all_point_types = []
        self.loaded_time_points = []
        self.vtp_files = []
        self.verbose = verbose

        if report_file is not None:
            self._load_report(report_file)
            self.sim_dir = os.path.dirname(os.path.abspath(report_file))
            if position_file is None:
                position_file = os.path.join(self.sim_dir, "position_buffer.txt")
        else:
            self.sim_dir = os.path.dirname(os.path.abspath(position_file))

        self._load_positions(position_file)
        self._find_vtp_files()

    def _load_report(self, report_file):
        self.report_data = json.load(open(report_file))
        print_(f"Loaded report_data:\n{self.report_data}")
        self.dt = float(self.report_data["dt"].split(" ")[0])
        self.duration = float(self.report_data["duration"].split(" ")[0])
        self.log_step = int(self.report_data["logstep"])

        neuron_time_points = np.linspace(
            0, self.duration, int(self.duration / self.dt) + 1
        )
        sibernetic_time_points = np.linspace(
            0, self.duration, int((self.duration / self.dt) / self.log_step)
        )
        print_(
            "Simulation dt: %s ms, duration: %s ms, times simulated (%i): %s; "
            "sibernetic logged times (%i): %s"
            % (
                self.dt,
                self.duration,
                len(neuron_time_points),
                neuron_time_points,
                len(sibernetic_time_points),
                sibernetic_time_points,
            )
        )

    def _load_positions(self, position_file):
        points = {}
        types = []
        line_count = 0
        pcount = 0
        time_count = 0
        log_step = None
        sampled = 1e6  # force first sample to be included
        first_pass_complete = False
        count_point_types = {}

        num_elastic = 0
        num_liquid = 0
        num_boundary = 0

        for line in open(position_file):
            ws = line.split()

            if line_count == 6:
                num_elastic = int(ws[0])
            elif line_count == 7:
                num_liquid = int(ws[0])
            elif line_count == 8:
                num_boundary = int(ws[0])
            elif line_count == 9:
                self.time_step = float(ws[0])
            elif line_count == 10:
                log_step = int(ws[0])

            if len(ws) == 4:
                type_ = float(ws[3])
                if type_ not in points:
                    points[type_] = []

                if not first_pass_complete:
                    count_point_types[type_] = count_point_types.get(type_, 0) + 1

                if self.swap_y_z:
                    points[type_].append([float(ws[1]), float(ws[0]), float(ws[2])])
                else:
                    points[type_].append([float(ws[0]), float(ws[1]), float(ws[2])])

                types.append(type_)

            if log_step is not None:
                pcount += 1

                if pcount == num_boundary + num_elastic + num_liquid:
                    first_pass_complete = True
                    sampled += 1
                    if self.verbose:
                        print_(
                            "End of one batch of %i total points (%i types), at line %i, time point: %i"
                            % (pcount, len(points), line_count, time_count)
                        )

                    if sampled < self.downsample:
                        print_(
                            "  -- Skipping sample %i due to downsampling factor %i"
                            % (sampled, self.downsample),
                            self.verbose,
                        )
                    else:
                        print_(
                            "  -- Including sample %i, downsampling factor %i"
                            % (sampled, self.downsample),
                            self.verbose,
                        )
                        self.all_3D_points.append(points)
                        self.all_point_types.append(types)
                        sampled = 0

                        if self.dt is not None:
                            time_calculated = time_count * log_step * self.dt
                            self.loaded_time_points.append(time_calculated)
                            print_(
                                f"Time calculated as: {time_calculated}", self.verbose
                            )
                        else:
                            self.loaded_time_points.append(time_count)

                    points = {}
                    types = []
                    num_boundary = 0
                    pcount = 0
                    time_count += 1

            line_count += 1

        self.num_elastic = num_elastic
        self.num_liquid = num_liquid
        self.num_boundary = num_boundary
        if self.log_step is None:
            self.log_step = log_step

        print_(
            "\nLoaded positions with %i elastic, %i liquid and %i boundary points "
            "(%i total), over %i lines"
            % (
                num_elastic,
                num_liquid,
                num_boundary,
                num_elastic + num_liquid + num_boundary,
                line_count,
            )
        )
        print_(
            f"Num of time points loaded: {len(self.all_3D_points)} (total: {time_count})"
        )
        print_(f"Loaded time points: {self.loaded_time_points}")
        print_(f"Count of point types found: {dict(sorted(count_point_types.items()))}")

    def _find_vtp_files(self):
        self.vtp_files = sorted(
            [
                os.path.join(self.sim_dir, f)
                for f in os.listdir(self.sim_dir)
                if f.endswith(".vtp")
            ]
        )
        if self.vtp_files:
            print_(
                f"Found {len(self.vtp_files)} VTK *.vtp files: "
                f"[{self.vtp_files[0]},..., {self.vtp_files[-1]}]"
            )

    def has_worm_data(self):
        """True if report.json indicates a worm-body configuration."""
        return self.report_data is not None and "worm" in self.report_data.get(
            "configuration", ""
        )

    def has_muscle_activation_data(self):
        """True if a muscles_activity_buffer.txt file exists in sim_dir."""
        return os.path.isfile(os.path.join(self.sim_dir, "muscles_activity_buffer.txt"))

    def num_time_points(self):
        """Number of loaded time frames."""
        return len(self.all_3D_points)

    def get_points_at(self, time_index):
        """Return the points dict {type_: [[x,y,z], ...]} at time_index."""
        return self.all_3D_points[time_index]

    # Method to get the t, x, y, z positions of a single particle of a given type
    def get_particle_position(self, particle_type, particle_index=0):
        """
        Get the time series of positions for a single particle of a given type.

        Parameters
        ----------
        particle_type : float
            The type of particle to track (e.g., 2.1 for elastic).
        particle_index : int
            The index of the particle of the given type to track (default 0).
        """
        positions = []
        times = []
        for time_idx, points in enumerate(self.all_3D_points):
            if particle_type in points and len(points[particle_type]) > particle_index:
                positions.append(points[particle_type][particle_index])
                times.append(self.loaded_time_points[time_idx])
        return times, np.array(positions)
