"""
Helper script to visualise particle positions and muscle
activations logged by the simulator. Used by `sibernetic_c302.py`
but can be invoked directly.
"""

import sys


def print_(msg):
    print("c302 plot >>> %s" % (msg))


"""
    Plot the positions from a saved Sibernetic position_buffer.txt file
"""


def plot_positions(pos_file_name, rate_to_plot=100, save_figure=True, show_plot=True):
    postions_file = open(pos_file_name)

    rate_to_plot = max(1, int(rate_to_plot))

    index = 0

    xmin = -1
    xmax = -1
    ymin = -1
    ymax = -1
    zmin = -1
    zmax = -1

    import matplotlib.pyplot as plt

    print_("Loading: %s" % pos_file_name)

    fig = plt.figure()

    ax = fig.add_subplot(111)
    # plt.xlim([-350, 450]) # based on standard configuration
    # plt.ylim([-50, 750])

    max_lines = 2e121
    max_time_ms = 10000

    # colors = {2.1:'blue', 2.2:'red',1.1:'black',3.0:'green'}
    points_plotted = 0

    a = b = c = dt = steps_per_frame = 0

    frame = 0
    in_frame = 0
    xs = []
    ys = []

    num_plotted_frames = 0

    for line in postions_file:
        if index > max_lines:
            print_("Finished parsing file, as max number of lines reached!")
            break

        if index == 0:
            xmin = float(line)
        elif index == 1:
            xmax = float(line)
        elif index == 2:
            ymin = float(line)  # noqa: F841
        elif index == 3:
            ymax = float(line)  # noqa: F841
        elif index == 4:
            zmin = float(line)
        elif index == 5:
            zmax = float(line)
        elif index == 6:
            a = float(line)
        elif index == 7:
            b = float(line)
        elif index == 8:
            c = float(line)
        elif index == 9:
            dt = float(line)
        elif index == 10:
            steps_per_frame = float(line)

        elif index == 11:
            x = [xmin, xmin, xmax, xmax, xmin]
            y = [zmin, zmax, zmax, zmin, zmin]
            # plt.plot(x,y,'-', color='grey')

        elif index >= 12:
            t_ms = frame * steps_per_frame * dt * 1000

            if t_ms > max_time_ms:
                print_(
                    "Finished parsing file, as max time (%s ms) reached!" % max_time_ms
                )

                break

            plot_frame = frame % rate_to_plot == 0

            if plot_frame:
                w = line.split()
                m = float(w[3])
                if m > 2 and m < 3:
                    x = float(w[0])
                    y = float(w[1])  # noqa: F841
                    z = float(w[2])
                    # plt.plot(x,z,'.',color=colors[m])
                    xs.append(x + num_plotted_frames * 30)
                    ys.append(z)
                    points_plotted += 1

            in_frame += 1
            if in_frame == a + b + c - 1:
                if plot_frame:
                    print_(
                        " >> Plotting frame %i at %s ms; line %i: %s...\n"
                        % (num_plotted_frames, t_ms, index, line)
                    )
                    ax.plot(xs, ys, ".", markersize=1)
                    ax.axis("equal")
                    num_plotted_frames += 1
                    if num_plotted_frames % 3 == 1:
                        time = (
                            "%sms" % t_ms
                            if not t_ms == int(t_ms)
                            else "%sms" % int(t_ms)
                        )
                        ax.text(
                            50 + ((num_plotted_frames - 1) * 30), 510, time, fontsize=12
                        )

                frame += 1
                in_frame = 0
                print_(
                    "New positions (#%i) at time %s ms; line %i: %s"
                    % (frame, t_ms, index, line)
                )
                xs = []
                ys = []

        index += 1

    print_(
        "Loaded: %s points from %s, showing %s points in %i plots"
        % (index, pos_file_name, points_plotted, num_plotted_frames)
    )

    if save_figure:
        plt.savefig("%s.png" % pos_file_name, bbox_inches="tight")

    if show_plot:
        plt.show()


def plot_muscle_activity(
    muscle_file_name, dt, logstep, save_figure=True, show_plot=True
):
    muscle_file = open(muscle_file_name)

    count = 0
    a = []
    times = []
    muscle_names = []

    quadrant0 = "MDR"
    quadrant1 = "MVR"
    quadrant2 = "MVL"
    quadrant3 = "MDL"
    qs = [quadrant0, quadrant1, quadrant2, quadrant3]

    activations = {}

    for q in qs:
        activations[q] = {}
        for i in range(24):
            activations[q][i] = []
            muscle_names.append("%s%s" % (q, i))
    max_act = 0
    for line in muscle_file:
        acts = [float(w) for w in line.split()]

        print_(
            "Found %s activation values (%s,...,%s) at line %s"
            % (len(acts), acts[0], acts[-1], count)
        )
        a.append(acts)

        for qi in range(len(qs)):
            q = qs[qi]
            for mi in range(24):
                index = qi * 24 + mi
                activations[q][mi].append(acts[index])
                max_act = max(max_act, acts[index])

        times.append(count * dt * logstep / 1000)

        count += 1

    import matplotlib.pyplot as plt

    fig, ax0 = plt.subplots(4, sharex=True, sharey=True)

    info = "Muscle activation values per quadrant"
    fig.canvas.manager.set_window_title(info)

    ax0[0].set_title(
        "Muscle activation values per quadrant - " + quadrant0, size="small"
    )
    ax0[0].set_ylabel("muscle #")
    ax0[1].set_title(quadrant3, size="small")
    ax0[1].set_ylabel("muscle #")
    ax0[2].set_title(quadrant1, size="small")
    ax0[2].set_ylabel("muscle #")
    ax0[3].set_title(quadrant2, size="small")
    ax0[3].set_ylabel("muscle #")
    ax0[3].set_xlabel("time (s)" % ())

    arr0 = []
    arr1 = []
    arr2 = []
    arr3 = []
    for i in range(24):
        arr0.append(activations[quadrant0][i])
        arr1.append(activations[quadrant3][i])
        arr2.append(activations[quadrant1][i])
        arr3.append(activations[quadrant2][i])

    ax0[0].imshow(arr0, interpolation="none", aspect="auto", vmin=0, vmax=max_act)
    ax0[1].imshow(arr1, interpolation="none", aspect="auto", vmin=0, vmax=max_act)
    ax0[2].imshow(arr2, interpolation="none", aspect="auto", vmin=0, vmax=max_act)
    ax0[3].imshow(arr3, interpolation="none", aspect="auto", vmin=0, vmax=max_act)
    # fig.colorbar(im0)

    ax1 = plt.gca()
    xt = ax1.get_xticks()
    time_ticks = [times[int(ti)] if (ti >= 0 and ti < len(times)) else 0 for ti in xt]
    ax1.set_xticklabels(time_ticks)
    # print ax1.get_yticks()
    # ax1.set_yticks([i-1 for i in ax1.get_yticks()])
    # print ax1.get_yticks()

    if save_figure:
        plt.savefig("%s.png" % muscle_file_name, bbox_inches="tight")
    """
    aa = np.array(a).transpose()
    fig, ax = plt.subplots()

    #plot0 = ax.pcolormesh(aa)

    plot0 = plt.imshow(aa, interpolation='nearest', aspect='auto')
    ax = plt.gca();


    ax.set_yticks(np.arange(aa.shape[0]) + 0.5, minor=False)
    ax.set_yticklabels(muscle_names, size='x-small')

    ax1 = plt.gca();
    xt = ax1.get_xticks()
    time_ticks = [times[int(ti)] if (ti>=0 and ti<len(times)) else 0 for ti in xt]
    ax1.set_xticklabels(time_ticks)

    fig.colorbar(plot0)

    if save_figure:
        plt.savefig('%s0.png'%muscle_file_name,bbox_inches='tight')"""

    if show_plot:
        plt.show()

    return activations, times


if __name__ == "__main__":
    if len(sys.argv) == 2:
        pos_file_name = sys.argv[1]
    else:
        pos_file_name = "buffers/position_buffer.txt"

    plot_positions(pos_file_name, rate_to_plot=50, show_plot=True, save_figure=True)
