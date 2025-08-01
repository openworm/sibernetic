"""
Utility module used by the C++ simulator to generate
synthetic muscle activation patterns. Can also be run
standalone to visualize these waves.
"""

import math
import sys

import numpy as np


muscle_row_count = 24

default_time_per_step = 0.000001  #  s
time_per_step = default_time_per_step  #  s

quadrant0 = "MDR"
quadrant1 = "MVR"
quadrant2 = "MVL"
quadrant3 = "MDL"

colours = {}
colours[quadrant0] = "#000000"
colours[quadrant1] = "#00ff00"
colours[quadrant2] = "#0000ff"
colours[quadrant3] = "#ff0000"


def print_(msg):
    pre = "Python >> "
    print("%s %s" % (pre, msg.replace("\n", "\n" + pre)))


"""

Get list of muscle names in same order as waves generated below.
Based on info here:
https://github.com/openworm/Smoothed-Particle-Hydrodynamics/blob/3da1edc3b018c2e5c7c1a25e2f8d44b54b1a1c47/src/owWorldSimulation.cpp#L475

"""


def get_muscle_names():
    names = []
    for i in range(muscle_row_count):
        names.append(get_muscle_name(quadrant0, i))
        names.append(get_muscle_name(quadrant1, i))
        names.append(get_muscle_name(quadrant2, i))
        names.append(get_muscle_name(quadrant3, i))

    return names


def get_muscle_name(quadrant, index):
    return "%s%s" % (quadrant, index + 1 if index > 8 else ("0%i" % (index + 1)))


def parallel_waves(
    n=muscle_row_count,  # 24 for our first test?
    step=0,
    phi=math.pi,
    amplitude=1,
    # velocity=0.000008):
    velocity_s=0.000015 * 3.7 * 1.94 / 1.76,  # swimming
    velocity_c=0.000015 * 0.72,
):  # crawling // 0.9 // 0.65
    """
    Array of two travelling waves, second one starts
    half way through the array
    """
    j = n / 2
    max_muscle_force_coeff = 1.0

    # "<" = first 6 seconds crawling, then swimming
    # ">" = first 6 seconds swimming, then crawling
    if step > 1200000:
        velocity = 4 * 0.000015 * 0.72  # crawling
        max_muscle_force_coeff = 1.0
        row_positions = np.linspace(0, 2.97 * math.pi, int(j))
        wave_m = np.linspace(1, 0.6, int(j))
    else:
        velocity = 4 * 0.000015 * 3.7  # swimming
        max_muscle_force_coeff = 0.575
        row_positions = np.linspace(0, 0.81 * math.pi, int(j))
        # wave_m = [0.8,0.7,0.8,0.93,1.0,1.0,1.0,1.0,0.93,0.8,0.6,0.4]
        # wave_m = [0.8,0.7,0.8,0.93,1.0,1.0,1.0,1.0,0.93,0.8,0.65,0.5]
        # wave_m = [0.8,0.7,0.8,0.93,1.0,1.0,1.0,0.93,0.8,0.6,0.4,0.2]
        # wave_m = [0.81,0.90,0.97,1.00,0.99,0.95,0.88,0.78,0.65,0.50,0.33,0.15]
        # wave_m = np.linspace(1,1,j)
        wave_m = [
            0.81,
            0.90,
            0.97,
            1.00,
            0.99,
            0.95,
            0.88,
            0.78,
            0.65,
            0.53,
            0.40,
            0.25,
        ]  # 6
        # wave_m = [0.81,0.90,0.97,1.00,0.99,0.95,0.90,0.83,0.75,0.65,0.55,0.45] 	#7

    if n % 2 != 0:
        raise NotImplementedError("Currently only supports even number of muscles!")

    wave_1 = map(math.sin, (row_positions - velocity * step))
    wave_2 = map(math.sin, (row_positions - velocity * step + (math.pi)))

    normalize_sine = lambda x: abs(x * (x > 0))  # (x + 1)/2   # noqa: E731
    wave_1 = map(normalize_sine, wave_1)
    wave_2 = map(normalize_sine, wave_2)

    ###### sinusoidal signal correction ##################################
    # normalize_sine = lambda x : x*x
    # wave_1 = map(normalize_sine, wave_1)
    # wave_2 = map(normalize_sine, wave_2)

    wave_1 = map(lambda x, y: max_muscle_force_coeff * x * y, wave_1, wave_m)
    wave_2 = map(lambda x, y: max_muscle_force_coeff * x * y, wave_2, wave_m)
    ###### smooth start###################################################
    if step < (10000 / 4):
        normalize_sine = lambda x: x * step / (10000 / 4)  # noqa: E731
        wave_1 = map(normalize_sine, wave_1)
        wave_2 = map(normalize_sine, wave_2)
    ######################################################################

    double_wave_1 = []
    double_wave_2 = []

    for i in wave_1:
        double_wave_1.append(i)
        double_wave_1.append(i)
    for i in wave_2:
        double_wave_2.append(i)
        double_wave_2.append(i)

    return (double_wave_1, double_wave_2)


class MuscleSimulation:
    def __init__(self, increment=1.0):
        self.increment = increment
        self.step = 0

    def set_timestep(self, dt):
        pass

    def run(self, skip_to_time=0, do_plot=True):
        self.contraction_array = parallel_waves(step=self.step)
        self.step += self.increment
        # if (self.step>1000000):
        #    self.increment = -1.0
        # else:
        #    self.increment = self.increment

        return list(
            np.concatenate(
                [
                    self.contraction_array[0],
                    self.contraction_array[1],
                    self.contraction_array[1],
                    self.contraction_array[0],
                ]
            )
        )

    def save_results(self):
        print_("MuscleSimulation does NOT save results")


class C302NRNSimulation:
    max_ca = 4e-7
    max_ca_found = -1

    def __init__(self, tstop=100, dt=0.005, activity_file=None, verbose=True):
        # from LEMS_c302_C1_Full_nrn import NeuronSimulation
        # from LEMS_c302_nrn import NeuronSimulation

        print("Initialising C302NRNSimulation...")

        self.tstop = tstop
        self.verbose = verbose

        # self.ns = NeuronSimulation(tstop, dt)
        # print_("Initialised C302NRNSimulation of length %s ms and dt = %s ms..."%(tstop,dt))

    def set_timestep(self, dt):
        print("Setting timestep to %s..." % dt)

        dt = float("{:0.1e}".format(dt)) * 1000.0  # memory issue fix

        try:
            from LEMS_c302_nrn import NeuronSimulation
            import neuron

            self.h = neuron.h
        except Exception as e:
            raise Exception("Python import error in C302NRNSimulation: %s.." % e)

        self.ns = NeuronSimulation(self.tstop, dt)
        print_(
            "Initialised C302NRNSimulation of length %s ms and dt = %s ms..."
            % (self.tstop, dt)
        )

    def save_results(self):
        print_("> Saving results at time: %s" % self.h.t)

        self.ns.save_results()

    def run(self, skip_to_time=-1):
        print_("> Current NEURON time: %s ms" % self.h.t)

        self.ns.advance()

        print_("< Current NEURON time: %s ms" % self.h.t)

        values = []
        vars_read = []
        for i in range(24):
            var = "a_MDR%s" % (i + 1 if i > 8 else ("0%i" % (i + 1)))
            try:
                val = getattr(self.h, var)[0].soma.cai
            except AttributeError as e:
                print(
                    "Problem passing neuronal output of %s to muscle in Sibernetic: %s"
                    % (var, e)
                )
                continue
                val = 0
            scaled_val = self._scale(val)
            values.append(scaled_val)
            vars_read.append(var)
        for i in range(24):
            var = "a_MVR%s" % (i + 1 if i > 8 else ("0%i" % (i + 1)))
            if i == 23:
                var = "a_MVR23"
            try:
                val = getattr(self.h, var)[0].soma.cai
            except AttributeError as e:
                print(
                    "Problem passing neuronal output of %s to muscle in Sibernetic: %s"
                    % (var, e)
                )
                val = 0
                continue
            scaled_val = self._scale(val)
            values.append(scaled_val)
            vars_read.append(var)
        for i in range(24):
            var = "a_MVL%s" % (i + 1 if i > 8 else ("0%i" % (i + 1)))
            try:
                val = getattr(self.h, var)[0].soma.cai
            except AttributeError as e:
                if var == "a_MVL24":
                    extra = (
                        "Note: not an issue as no muscle MVL24 in the real C. elegans"
                    )
                else:
                    extra = ""
                print(
                    "Problem passing output of %s to muscle in Sibernetic: %s %s"
                    % (var, e, extra)
                )

                val = 0
                continue
            scaled_val = self._scale(val)
            values.append(scaled_val)
            vars_read.append(var)
        for i in range(24):
            var = "a_MDL%s" % (i + 1 if i > 8 else ("0%i" % (i + 1)))
            try:
                val = getattr(self.h, var)[0].soma.cai
            except AttributeError as e:
                print(
                    "Problem passing neuronal output of %s to muscle in Sibernetic: %s"
                    % (var, e)
                )
                val = 0
                continue
            scaled_val = self._scale(val)
            values.append(scaled_val)
            vars_read.append(var)

        if self.verbose:
            print_("Returning %s values: %s; %s" % (len(values), values, vars_read))
        return values

    def _scale(self, ca, print_it=False):
        self.max_ca_found = max(ca, self.max_ca_found)
        scaled = min(1, (ca / self.max_ca))
        if print_it:
            print_(
                "- Scaling %s to %s (max found: %s)" % (ca, scaled, self.max_ca_found)
            )
        return scaled


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    print_("This script is used by the Sibernetic C++ application")
    print_(
        "Running it directly in Python will only plot the waves being generated for sending to the muscle cells..."
    )

    try_c302 = "-c302dat" in sys.argv
    try_c302_nrn = "-c302nrn" in sys.argv
    testnrn = "-testnrn" in sys.argv

    skip_to_time = 0

    max_time = 2.0  # s

    time_per_step = 0.001  #  s
    increment = time_per_step / default_time_per_step
    num_plots = 2

    ms = MuscleSimulation(increment=increment)

    if try_c302:
        # ms = C302Simulation('configuration/test/c302/c302_B_Muscles.muscles.activity.dat', scale_to_max=True)
        ms = C302Simulation(  # noqa
            "configuration/test/c302/c302_C1_Muscles.muscles.activity.dat",  # noqa
            scale_to_max=True,  # noqa
        )  # noqa
        # ms = C302Simulation('../../../neuroConstruct/osb/invertebrate/celegans/CElegansNeuroML/CElegans/pythonScripts/c302/TestMuscles.activity.dat')
        # ms = C302Simulation('../../neuroConstruct/osb/invertebrate/celegans/CElegansNeuroML/CElegans/pythonScripts/c302/c302_B_Oscillator.muscles.activity.dat')
        skip_to_time = 0.05
        max_time = 0.2

    elif try_c302_nrn or testnrn:
        dt = 0.1  # ms
        max_time = 0.5  # s
        maxt = max_time * 1000

        time_per_step = dt / 1000  #  s
        increment = time_per_step / default_time_per_step

        ms = C302NRNSimulation(tstop=maxt, dt=dt, verbose=False)

    activation = {}
    row = "11"
    row_int = int(row)
    m0 = "%s%s" % (quadrant0, row)
    m1 = "%s%s" % (quadrant1, row)
    m2 = "%s%s" % (quadrant2, row)
    m3 = "%s%s" % (quadrant3, row)

    for m in get_muscle_names():
        activation[m] = []
    times = []

    num_steps = int(max_time / time_per_step)
    steps_between_plots = int(num_steps / num_plots)

    show_all = True

    if testnrn:
        for step in range(num_steps):
            ms.run()

        ms.save_results()

        quit()

    for step in range(num_steps):
        t = step * time_per_step

        output = ms.run(skip_to_time=skip_to_time)

        for i in range(muscle_row_count):
            mq0 = get_muscle_name(quadrant0, i)
            activation[mq0].append(output[i])
            mq1 = get_muscle_name(quadrant1, i)
            activation[mq1].append(output[i + muscle_row_count])
            mq2 = get_muscle_name(quadrant2, i)
            activation[mq2].append(output[i + muscle_row_count * 2])
            mq3 = get_muscle_name(quadrant3, i)
            activation[mq3].append(output[i + muscle_row_count * 3])

        times.append(t)

        if step == 0 or step % steps_between_plots == 0:
            print_("At step %s (%s s)" % (step, t))
            if show_all:
                figV = plt.figure()
                figV.suptitle("Muscle activation waves at step %s (%s s)" % (step, t))
                plV = figV.add_subplot(111, autoscale_on=True)

                plV.plot(
                    output[0:muscle_row_count],
                    label="%s*" % quadrant0,
                    color=colours[quadrant0],
                    linestyle="-",
                    marker="o",
                )
                plV.plot(
                    output[muscle_row_count : 2 * muscle_row_count],
                    label="%s*" % quadrant1,
                    color=colours[quadrant1],
                    linestyle="-",
                    marker="o",
                )
                plV.plot(
                    output[2 * muscle_row_count : 3 * muscle_row_count],
                    label="%s*" % quadrant2,
                    color=colours[quadrant2],
                    linestyle="-",
                )
                plV.plot(
                    output[3 * muscle_row_count : 4 * muscle_row_count],
                    label="%s*" % quadrant3,
                    color=colours[quadrant3],
                    linestyle="-",
                )

                plV.legend()

    if show_all:
        fig0 = plt.figure()
        fig0.suptitle(
            "Muscle activation waves of [%s, %s, %s, %s] vs time" % (m0, m1, m2, m3)
        )
        pl0 = fig0.add_subplot(111, autoscale_on=True)
        pl0.plot(
            times, activation[m0], label=m0, color=colours[quadrant0], linestyle="-"
        )
        pl0.plot(
            times, activation[m1], label=m1, color=colours[quadrant1], linestyle="-"
        )
        pl0.plot(
            times, activation[m2], label=m2, color=colours[quadrant2], linestyle="--"
        )
        pl0.plot(
            times, activation[m3], label=m3, color=colours[quadrant3], linestyle="--"
        )

        pl0.legend()

    f, a = plt.subplots(4, sharex=True, sharey=True)

    a[0].set_title(quadrant0)
    a[0].set_ylabel("muscle #")
    a[1].set_title(quadrant3)
    a[1].set_ylabel("muscle #")
    a[2].set_title(quadrant1)
    a[2].set_ylabel("muscle #")
    a[3].set_title(quadrant2)
    a[3].set_ylabel("muscle #")
    a[3].set_xlabel("time step (%s s -> %s s)" % (skip_to_time, max_time))

    arr0 = []
    arr1 = []
    arr2 = []
    arr3 = []
    for i in range(muscle_row_count):
        arr0.append(activation[get_muscle_name(quadrant0, i)])
        arr1.append(activation[get_muscle_name(quadrant3, i)])
        arr2.append(activation[get_muscle_name(quadrant1, i)])
        arr3.append(activation[get_muscle_name(quadrant2, i)])

    a[0].imshow(arr0, interpolation="none", aspect="auto")
    a[1].imshow(arr1, interpolation="none", aspect="auto")
    a[2].imshow(arr2, interpolation="none", aspect="auto")
    a[3].imshow(arr3, interpolation="none", aspect="auto")

    # plt.colorbar()

    plt.show()
